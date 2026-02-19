use std::{
    collections::VecDeque,
    fmt::Display,
    fs::File,
    io::{stdout, BufReader, BufWriter, Read, Write},
    path::Path,
    sync::{Arc, Condvar, Mutex},
    thread,
};

use anyhow::{anyhow, Context, Result};
mod chunker;
use crate::chunker::{FastxChunker, SEQ_SENTINEL};

/// Re-export for benchmarks. Not part of the public API.
#[doc(hidden)]
pub mod chunker_for_bench {
    pub use crate::chunker::FastxChunker;
}

#[cfg(test)]
mod worker_test_hooks {
    use std::collections::HashSet;
    use std::sync::{Mutex, OnceLock};

    static WORKER_IDS: OnceLock<Mutex<HashSet<usize>>> = OnceLock::new();

    fn store() -> &'static Mutex<HashSet<usize>> {
        WORKER_IDS.get_or_init(|| Mutex::new(HashSet::new()))
    }

    pub(super) fn reset_worker_ids() {
        store().lock().unwrap().clear();
    }

    pub(super) fn record_worker_id(id: usize) {
        store().lock().unwrap().insert(id);
    }

    pub(super) fn worker_id_count() -> usize {
        store().lock().unwrap().len()
    }
}

type PerThreadVec<T> = Arc<Vec<Mutex<Vec<T>>>>;

/// A constant array mapping ASCII characters to their 2-bit DNA representation index.
///
/// 'A'|'a' -> 0, 'C'|'c' -> 1, 'G'|'g' -> 2, 'T'|'t' -> 3.
/// All other characters map to `usize::MAX`, indicating an invalid base for k-mer counting.
pub const KEY_MAP: [usize; 256] = {
    let mut key_map = [usize::MAX; 256];
    key_map[b'A' as usize] = 0;
    key_map[b'a' as usize] = 0;
    key_map[b'C' as usize] = 1;
    key_map[b'c' as usize] = 1;
    key_map[b'G' as usize] = 2;
    key_map[b'g' as usize] = 2;
    key_map[b'T' as usize] = 3;
    key_map[b't' as usize] = 3;
    key_map
};

/// Returns a reader capable of handling both compressed and uncompressed files (uses
/// [niffler](https://crates.io/crates/niffler) for sniffing and decompression).
///
/// This function uses [niffler](https://crates.io/crates/niffler) to detect compression format
/// (like gzip, bzip2) based on the file's magic numbers. If the file is too short for detection
/// (e.g., empty or very small), it safely falls back to a standard `BufReader`. This allows
/// transparent reading of potentially compressed FASTA/FASTQ files.
pub fn get_compression_agnostic_reader<P: AsRef<Path> + Display>(
    path: P,
) -> Result<Box<dyn Read + Send>> {
    // if the file has fewer than 5 bytes, `niffler` can't sniff the compression format and will
    // return a `FileTooShort` error; this could be due to
    // * an empty file
    // * a file containing only a single FASTA record with the ID consisting only of a single
    //   character and the sequence being empty
    // * a file containing a single sequence with a one-character ID and one-character sequence and
    //   missing newline character at the end
    // we don't want to fail at this stage in these cases and thus handle the `FileTooShort` error
    // separately
    let reader = match niffler::send::from_path(&path) {
        Ok(rdr) => Ok(rdr.0),
        Err(niffler::Error::FileTooShort) => {
            Ok(Box::new(BufReader::new(File::open(&path)?)) as Box<dyn Read + Send>)
        }
        Err(e) => Err(e).context(format!("Failed reading file \"{}\"", path)),
    }?;

    Ok(reader)
}

/// Converts a numerical k-mer index into the corresponding DNA sequence, writing into an existing
/// buffer instead of allocating.
///
/// Given an index (representing a k-mer encoded as a `usize` where each base takes 2 bits) and the
/// k-mer length `k`, this function reconstructs the original k-mer sequence as 'A', 'C', 'G', 'T'
/// in the first `k` bytes of `buf`. The buffer must have length >= k.
pub fn index_to_kmer_buf(index: usize, k: usize, buf: &mut [u8]) {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut idx = index;
    for i in (0..k).rev() {
        buf[i] = bases[idx & 3];
        idx >>= 2;
    }
}

/// Computes the index of the reverse complement directly from a k-mer index using bit manipulation.
///
/// Instead of converting an index to a k-mer string, computing the reverse complement string, and
/// converting back, this function operates directly on the 2-bit encoded index. Each base pair is
/// extracted, complemented (A<->T, C<->G via XOR with 3), and placed in reverse order.
///
/// This is equivalent to computing `dna::revcomp` on the k-mer string and converting back to an
/// index, but avoids all allocations.
pub fn revcomp_index(index: usize, k: usize) -> usize {
    let mut rc_idx = 0;
    let mut idx = index;
    for _ in 0..k {
        // extract last 2 bits, complement (XOR with 3), shift into result
        rc_idx = (rc_idx << 2) | ((idx & 3) ^ 3);
        idx >>= 2;
    }
    rc_idx
}

/// Converts a DNA k-mer sequence (as bytes) into its numerical index representation.
///
/// This function takes a k-mer and translates each base into its 2-bit value, packing them into a
/// `usize`. This index can be used for efficient storage and lookup in count arrays.
pub fn kmer_to_index(kmer: &[u8]) -> usize {
    let mut idx = 0;
    for b in kmer {
        idx = (idx << 2) + KEY_MAP[*b as usize];
    }
    idx
}

/// Enum to hold k-mer counts and configuration (variants for single- and multi-threaded).
///
/// This enum manages k-mer counting, supporting both single-threaded and multi-threaded execution
/// modes. It stores the k-mer length `k`, the counts (either a single vector or one vector per
/// thread), and the bitmask needed for setting the left-shifted bits to 0 during the rolling index
/// calculation.
pub enum KmerCounter {
    SingleThreaded {
        k: usize,
        counts: Vec<usize>,
        mask: usize,
    },
    MultiThreaded {
        k: usize,
        per_thread_counts: PerThreadVec<usize>,
        mask: usize,
        n_threads: usize,
    },
}

impl KmerCounter {
    /// Creates a new `KmerCounter` instance.
    ///
    /// Initializes the counter with the specified k-mer length and the desired number of threads.
    /// If `num_threads` is 1, it creates a `SingleThreaded` counter. Otherwise, it sets up a
    /// `MultiThreaded` counter with thread-local storage for counts.
    pub fn new(k: usize, num_threads: usize) -> Self {
        assert!(k > 0, "k must be >= 1");
        // define mask to set left-shifted bits to 0 (used below)
        let mask = 4usize.pow(k as u32) - 1;

        match num_threads {
            1 => Self::SingleThreaded {
                k,
                counts: vec![0; 4usize.pow(k as u32)],
                mask,
            },
            _ => Self::MultiThreaded {
                k,
                per_thread_counts: Arc::new(
                    (0..num_threads)
                        .map(|_| Mutex::new(vec![0; 4usize.pow(k as u32)]))
                        .collect(),
                ),
                n_threads: num_threads,
                mask,
            },
        }
    }

    /// Counts k-mers from one or more files, using the worker-pool chunking model when threaded.
    pub fn count_paths(
        &mut self,
        paths: &[String],
        chunk_size: usize,
        queue_bytes: usize,
    ) -> Result<()> {
        match self {
            Self::SingleThreaded { k, counts, mask } => {
                for path in paths {
                    let mut file = build_chunker(path, *k, chunk_size)?;
                    while let Some(chunk) = file.chunker.next_chunk()? {
                        count_chunk(*k, counts, *mask, &chunk);
                    }
                }
            }
            Self::MultiThreaded {
                k,
                per_thread_counts,
                mask,
                n_threads,
            } => {
                count_paths_worker_pool(
                    paths,
                    *k,
                    per_thread_counts,
                    *mask,
                    *n_threads,
                    chunk_size,
                    queue_bytes,
                )?;
            }
        }
        Ok(())
    }

    /// Consumes the `KmerCounter` and returns the final counts as a single `Vec<usize>`.
    ///
    /// If the counter was multi-threaded, this method aggregates the counts from all threads into a
    /// single vector. For single-threaded counters, it simply returns the existing vector. Returns
    /// an error if aggregation fails (e.g., due to locking issues or unexpected state).
    pub fn into_vec(self) -> Result<Vec<usize>> {
        match self {
            Self::SingleThreaded {
                k: _,
                counts,
                mask: _,
            } => Ok(counts),
            Self::MultiThreaded {
                k: _,
                per_thread_counts,
                mask: _,
                n_threads: _,
            } => {
                // sum up all the counts (into the first counts vector); we do this in an extra
                // scope to return the borrow before unwrapping the `Arc` below
                {
                    let mut first_thread_counts = per_thread_counts[0].lock().unwrap();
                    for thread_count in per_thread_counts.iter().skip(1) {
                        let thread_count = thread_count.lock().unwrap();
                        for (first, &other) in
                            first_thread_counts.iter_mut().zip(thread_count.iter())
                        {
                            *first += other;
                        }
                    }
                }

                // take ownership of the counts
                let per_thread_counts = Arc::try_unwrap(per_thread_counts).map_err(
                    // `Arc::try_unwrap` returns a `Result` containing the `Arc` itself instead of
                    // an `Error` and we need replace it with an `Error` if we want to use `?`
                    |_| anyhow!("Error combining counts: Arc still has multiple owners"),
                )?;

                let combined_counts = per_thread_counts
                    .into_iter()
                    .next()
                    .context("Error combining counts: no count vector found")?
                    .into_inner()
                    .context("Error combining counts: could not extract vector from `Mutex`")?;

                Ok(combined_counts)
            }
        }
    }

    /// Writes the k-mer counts to the specified output.
    ///
    /// Consumes the `KmerCounter`, aggregates canonical counts if necessary, and writes the results
    /// (tab-delimited) to `outfile` (or standard output if `None`).
    /// - If `canonical` is true, only the lexicographically smaller of a k-mer and its reverse
    ///   complement is reported, with counts summed. Palindromes are reported as is.
    /// - If `write_zeros` is false, k-mers with zero counts are omitted.
    ///
    /// Uses [`revcomp_index`] for allocation-free reverse complement lookup and [`index_to_kmer_buf`]
    /// with [`itoa::Buffer`] to avoid per-k-mer heap allocations in the write loop.
    pub fn write(self, outfile: Option<&String>, canonical: bool, write_zeros: bool) -> Result<()> {
        let k = match self {
            Self::SingleThreaded { k, .. } => k,
            Self::MultiThreaded { k, .. } => k,
        };

        // get the counts vector (this consumes `self`)
        let counts = self.into_vec()?;

        // create output file or write to STDOUT
        let out: Box<dyn Write> = match outfile {
            Some(outfile) => Box::new(File::create(outfile)?),
            None => Box::new(stdout()),
        };

        let mut out_writer = BufWriter::new(out);
        let mut kmer_buf = vec![0u8; k];
        let mut num_buf = itoa::Buffer::new();

        if canonical {
            // look up the reverse complement of the kmer and print the smaller of the two
            // (alongside the sum of the counts)
            for (index, &count) in counts.iter().enumerate() {
                let rc_index = revcomp_index(index, k);
                // only write if this kmer is canonical (i.e. index <= rc_index)
                if index > rc_index {
                    continue;
                }
                // for a palindromic kmer, write the counts as is; otherwise, write the sum of both
                // counts
                let combined_count = if index == rc_index {
                    count
                } else {
                    count + counts[rc_index]
                };
                if combined_count == 0 && !write_zeros {
                    continue;
                }
                index_to_kmer_buf(index, k, &mut kmer_buf);
                out_writer.write_all(&kmer_buf)?;
                out_writer.write_all(b"\t")?;
                out_writer.write_all(num_buf.format(combined_count).as_bytes())?;
                out_writer.write_all(b"\n")?;
            }
        } else {
            for (index, &count) in counts.iter().enumerate() {
                if count == 0 && !write_zeros {
                    continue;
                }
                index_to_kmer_buf(index, k, &mut kmer_buf);
                out_writer.write_all(&kmer_buf)?;
                out_writer.write_all(b"\t")?;
                out_writer.write_all(num_buf.format(count).as_bytes())?;
                out_writer.write_all(b"\n")?;
            }
        }
        Ok(())
    }
}

impl std::ops::Add<KmerCounter> for KmerCounter {
    type Output = Result<KmerCounter>;
    /// Adds the counts from two `KmerCounter` instances.
    ///
    /// This operation is currently only defined for combining two `SingleThreaded` counters with
    /// the same `k` value. It returns a new `SingleThreaded` counter containing the summed counts
    /// or returns an error if the counters are incompatible (different types, different `k`, or
    /// different count vector lengths).
    fn add(self, other: KmerCounter) -> Result<KmerCounter> {
        match (self, other) {
            (
                KmerCounter::SingleThreaded {
                    k,
                    mut counts,
                    mask,
                },
                KmerCounter::SingleThreaded {
                    k: other_k,
                    counts: other_counts,
                    mask: _,
                },
            ) => {
                if k != other_k {
                    return Err(anyhow!("Cannot add KmerCounters with different k values"));
                }
                if counts.len() != other_counts.len() {
                    return Err(anyhow!(
                        "Cannot add KmerCounters with different count lengths"
                    ));
                }
                for (left, &right) in counts.iter_mut().zip(other_counts.iter()) {
                    *left += right;
                }
                Ok(KmerCounter::SingleThreaded { k, counts, mask })
            }
            _ => Err(anyhow!(
                "Cannot add SingleThreaded and MultiThreaded KmerCounters"
            )),
        }
    }
}

struct FileState {
    chunker: FastxChunker<BufReader<Box<dyn Read + Send>>>,
}

struct FilePoolState {
    files: VecDeque<FileState>,
    active_readers: usize,
    done: bool,
}

struct FilePool {
    state: Mutex<FilePoolState>,
    cv: Condvar,
}

impl FilePool {
    fn new(files: Vec<FileState>) -> Self {
        Self {
            state: Mutex::new(FilePoolState {
                files: files.into(),
                active_readers: 0,
                done: false,
            }),
            cv: Condvar::new(),
        }
    }

    fn acquire(&self) -> Option<FileState> {
        let mut state = self.state.lock().unwrap();
        if let Some(file) = state.files.pop_front() {
            state.active_readers += 1;
            return Some(file);
        }
        if state.active_readers == 0 {
            state.done = true;
        }
        None
    }

    fn release(&self, file: Option<FileState>) {
        let mut state = self.state.lock().unwrap();
        if state.active_readers > 0 {
            state.active_readers -= 1;
        }
        if let Some(file) = file {
            state.files.push_back(file);
        }
        if state.files.is_empty() && state.active_readers == 0 {
            state.done = true;
        }
        self.cv.notify_all();
    }

    fn is_done(&self) -> bool {
        self.state.lock().unwrap().done
    }
}

struct ChunkQueueState {
    queue: VecDeque<Vec<u8>>,
    bytes: usize,
    closed: bool,
}

struct ChunkQueue {
    state: Mutex<ChunkQueueState>,
    cv: Condvar,
    capacity_bytes: usize,
}

impl ChunkQueue {
    fn new(capacity_bytes: usize) -> Self {
        Self {
            state: Mutex::new(ChunkQueueState {
                queue: VecDeque::new(),
                bytes: 0,
                closed: false,
            }),
            cv: Condvar::new(),
            capacity_bytes,
        }
    }

    fn try_push(&self, chunk: Vec<u8>) -> Result<(), Vec<u8>> {
        let mut state = self.state.lock().unwrap();
        if state.closed {
            return Ok(());
        }
        let next_bytes = state.bytes + chunk.len();
        if next_bytes > self.capacity_bytes {
            return Err(chunk);
        }
        state.bytes = next_bytes;
        state.queue.push_back(chunk);
        self.cv.notify_one();
        Ok(())
    }

    fn try_pop(&self) -> Option<Vec<u8>> {
        let mut state = self.state.lock().unwrap();
        let chunk = state.queue.pop_front()?;
        state.bytes = state.bytes.saturating_sub(chunk.len());
        Some(chunk)
    }

    fn wait_for_chunk(&self) -> Option<Vec<u8>> {
        let mut state = self.state.lock().unwrap();
        while state.queue.is_empty() && !state.closed {
            state = self.cv.wait(state).unwrap();
        }
        let chunk = state.queue.pop_front()?;
        state.bytes = state.bytes.saturating_sub(chunk.len());
        Some(chunk)
    }

    fn is_empty(&self) -> bool {
        self.state.lock().unwrap().queue.is_empty()
    }

    fn close(&self) {
        let mut state = self.state.lock().unwrap();
        state.closed = true;
        self.cv.notify_all();
    }
}

fn build_chunker(path: &str, k: usize, chunk_size: usize) -> Result<FileState> {
    let reader = get_compression_agnostic_reader(path)?;
    let buf_reader = BufReader::new(reader);
    let chunker = FastxChunker::new(buf_reader, k, chunk_size)?;
    Ok(FileState { chunker })
}

fn count_chunk(k: usize, counts: &mut [usize], mask: usize, chunk: &[u8]) {
    count_kmers_bytes(k, counts, mask, chunk, Some(SEQ_SENTINEL));
}

fn count_kmers_bytes(
    k: usize,
    counts: &mut [usize],
    mask: usize,
    bytes: &[u8],
    reset_on: Option<u8>,
) {
    if k == 0 || bytes.len() < k {
        return;
    }

    let mut skip = k - 1;
    let mut idx = 0usize;

    match reset_on {
        Some(reset) => {
            for &b in bytes {
                if b == reset {
                    skip = k - 1;
                    idx = 0;
                    continue;
                }

                idx = (idx << 2) & mask;
                let add = KEY_MAP[b as usize];
                if add == usize::MAX {
                    skip = k - 1;
                    continue;
                }
                idx += add;

                if skip == 0 {
                    counts[idx] += 1;
                } else {
                    skip -= 1;
                }
            }
        }
        None => {
            for &b in bytes {
                idx = (idx << 2) & mask;
                let add = KEY_MAP[b as usize];
                if add == usize::MAX {
                    skip = k - 1;
                    continue;
                }
                idx += add;

                if skip == 0 {
                    counts[idx] += 1;
                } else {
                    skip -= 1;
                }
            }
        }
    }
}

fn count_paths_worker_pool(
    paths: &[String],
    k: usize,
    per_thread_counts: &mut PerThreadVec<usize>,
    mask: usize,
    n_threads: usize,
    chunk_size: usize,
    queue_bytes: usize,
) -> Result<()> {
    let mut files = Vec::with_capacity(paths.len());
    for path in paths {
        files.push(build_chunker(path, k, chunk_size)?);
    }

    let file_pool = Arc::new(FilePool::new(files));
    let queue_capacity = queue_bytes.max(chunk_size * 2).max(1);
    let chunk_queue = Arc::new(ChunkQueue::new(queue_capacity));

    let mut handles = Vec::with_capacity(n_threads);
    for thread_id in 0..n_threads {
        let file_pool = Arc::clone(&file_pool);
        let chunk_queue = Arc::clone(&chunk_queue);
        let per_thread_counts = Arc::clone(per_thread_counts);

        let handle = thread::spawn(move || -> Result<()> {
            let mut counts = per_thread_counts[thread_id].lock().unwrap();

            loop {
                if let Some(chunk) = chunk_queue.try_pop() {
                    #[cfg(test)]
                    {
                        worker_test_hooks::record_worker_id(thread_id);
                    }
                    count_chunk(k, &mut counts, mask, &chunk);
                    continue;
                }

                if file_pool.is_done() && chunk_queue.is_empty() {
                    chunk_queue.close();
                    break;
                }

                if let Some(mut file) = file_pool.acquire() {
                    loop {
                        match file.chunker.next_chunk()? {
                            Some(chunk) => {
                                if let Err(chunk) = chunk_queue.try_push(chunk) {
                                    // Queue full — count this chunk locally instead of blocking
                                    #[cfg(test)]
                                    {
                                        worker_test_hooks::record_worker_id(thread_id);
                                    }
                                    count_chunk(k, &mut counts, mask, &chunk);
                                    file_pool.release(Some(file));
                                    break;
                                }
                            }
                            None => {
                                file_pool.release(None);
                                break;
                            }
                        }
                    }
                    continue;
                }

                if file_pool.is_done() && chunk_queue.is_empty() {
                    chunk_queue.close();
                    break;
                }

                if let Some(chunk) = chunk_queue.wait_for_chunk() {
                    #[cfg(test)]
                    {
                        worker_test_hooks::record_worker_id(thread_id);
                    }
                    count_chunk(k, &mut counts, mask, &chunk);
                    continue;
                }

                if file_pool.is_done() {
                    chunk_queue.close();
                    break;
                }
            }

            Ok(())
        });

        handles.push(handle);
    }

    for handle in handles {
        handle
            .join()
            .map_err(|_| anyhow!("Worker thread panicked"))??;
    }

    chunk_queue.close();

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::alphabets::dna;
    use std::fs;
    use std::str::from_utf8;

    // helper function for tests.
    fn count_kmers_in_sequence(seq: &[u8], k: usize) -> Vec<usize> {
        let mask = 4usize.pow(k as u32) - 1;
        let mut counts = vec![0; 4usize.pow(k as u32)];
        count_kmers_bytes(k, &mut counts, mask, seq, None);
        counts
    }

    fn kmer_from_index(index: usize, k: usize) -> Vec<u8> {
        let mut buf = vec![0u8; k];
        index_to_kmer_buf(index, k, &mut buf);
        buf
    }

    fn write_temp_fasta(id: &str, seq: &str) -> String {
        let mut path = std::env::temp_dir();
        path.push(format!("smakcr_worker_test_{}.fa", id));
        fs::write(&path, format!(">seq\n{}\n", seq)).unwrap();
        path.to_str().unwrap().to_string()
    }

    #[test]
    fn test_count_kmers_in_sequence() {
        let k = 3;
        let seq: &[u8] = b"AcGTaCgAaaA";

        let counts = count_kmers_in_sequence(seq, k);

        let acg_index = 0b00_01_10; // ACG; comes up twice
        let cgt_index = 0b01_10_11; // CGT
        let gta_index = 0b10_11_00; // GTA
        let tac_index = 0b11_00_01; // TAC
        let cga_index = 0b01_10_00; // CGA
        let gaa_index = 0b10_00_00; // GAA
        let aaa_index = 0b00_00_00; // AAA; comes up twice

        assert_eq!(counts[acg_index], 2);
        assert_eq!(counts[cgt_index], 1);
        assert_eq!(counts[gta_index], 1);
        assert_eq!(counts[tac_index], 1);
        assert_eq!(counts[cga_index], 1);
        assert_eq!(counts[gaa_index], 1);
        assert_eq!(counts[aaa_index], 2);

        // make sure all other counts are 0
        assert_eq!(counts.iter().sum::<usize>(), 9);
    }

    #[test]
    fn test_count_kmers_with_short_sequence() {
        // NOTE: this actually only tests the implementation for `&[u8]` and not for the other types
        let k = 5;
        let counts = count_kmers_in_sequence(b"ACG", k);

        assert!(counts.iter().all(|&count| count == 0));
    }

    #[test]
    fn test_index_to_kmer_buf() {
        let k = 3;

        let acg_index = 0b00_01_10; // ACG
        let cgt_index = 0b01_10_11; // CGT
        let gta_index = 0b10_11_00; // GTA
        let tac_index = 0b11_00_01; // TAC
        let cga_index = 0b01_10_00; // CGA
        let gaa_index = 0b10_00_00; // GAA
        let aaa_index = 0b00_00_00; // AAA

        assert_eq!(kmer_from_index(acg_index, k), b"ACG");
        assert_eq!(kmer_from_index(cgt_index, k), b"CGT");
        assert_eq!(kmer_from_index(gta_index, k), b"GTA");
        assert_eq!(kmer_from_index(tac_index, k), b"TAC");
        assert_eq!(kmer_from_index(cga_index, k), b"CGA");
        assert_eq!(kmer_from_index(gaa_index, k), b"GAA");
        assert_eq!(kmer_from_index(aaa_index, k), b"AAA");
    }

    #[test]
    fn test_unknown_bases() {
        let k = 3;
        let counts = count_kmers_in_sequence(b"AaAxcGtx", k);

        let aaa_index = 0b00_00_00; // AAA
        let cgt_index = 0b01_10_11; // CGT
        assert_eq!(counts[aaa_index], 1);
        assert_eq!(counts[cgt_index], 1);
        // make sure all other counts are 0
        assert_eq!(counts.iter().sum::<usize>(), 2);
    }

    #[test]
    fn test_revcomp_index_basic() {
        // ACG (index 0b00_01_10 = 6) -> revcomp is CGT (index 0b01_10_11 = 27)
        assert_eq!(revcomp_index(0b00_01_10, 3), 0b01_10_11);
        // CGT -> ACG
        assert_eq!(revcomp_index(0b01_10_11, 3), 0b00_01_10);
    }

    #[test]
    fn test_revcomp_index_palindrome() {
        // For k=2, AT (0b00_11 = 3) -> revcomp AT (0b00_11 = 3)
        assert_eq!(revcomp_index(0b00_11, 2), 0b00_11);
        // GC (0b10_01 = 9) -> GC (0b10_01 = 9)
        assert_eq!(revcomp_index(0b10_01, 2), 0b10_01);
    }

    #[test]
    fn test_revcomp_index_single_base() {
        // A (0) -> T (3)
        assert_eq!(revcomp_index(0, 1), 3);
        // C (1) -> G (2)
        assert_eq!(revcomp_index(1, 1), 2);
        // G (2) -> C (1)
        assert_eq!(revcomp_index(2, 1), 1);
        // T (3) -> A (0)
        assert_eq!(revcomp_index(3, 1), 0);
    }

    #[test]
    fn test_revcomp_index_matches_string_method() {
        // exhaustively verify that revcomp_index matches the string-based method for all k-mers
        for k in 1..=6 {
            let n_kmers = 4usize.pow(k as u32);
            for idx in 0..n_kmers {
                let kmer = kmer_from_index(idx, k);
                let rc_string = dna::revcomp(&kmer);
                let expected = kmer_to_index(&rc_string);
                assert_eq!(
                    revcomp_index(idx, k),
                    expected,
                    "mismatch for k={}, idx={}, kmer={}",
                    k,
                    idx,
                    from_utf8(&kmer).unwrap()
                );
            }
        }
    }

    #[test]
    fn test_revcomp_index_involution() {
        // revcomp(revcomp(x)) == x for all k-mers
        for k in 1..=5 {
            let n_kmers = 4usize.pow(k as u32);
            for idx in 0..n_kmers {
                assert_eq!(revcomp_index(revcomp_index(idx, k), k), idx);
            }
        }
    }

    #[test]
    fn test_worker_pool_uses_multiple_workers() {
        let k = 5;
        let n_threads = 4;
        let chunk_size = 64;
        let queue_bytes = 256;

        let seq = "ACGT".repeat(4096);
        let paths = vec![
            write_temp_fasta("a", &seq),
            write_temp_fasta("b", &seq),
            write_temp_fasta("c", &seq),
            write_temp_fasta("d", &seq),
        ];

        worker_test_hooks::reset_worker_ids();
        let mut counter = KmerCounter::new(k, n_threads);
        counter
            .count_paths(&paths, chunk_size, queue_bytes)
            .unwrap();

        let worker_count = worker_test_hooks::worker_id_count();
        for path in paths {
            let _ = fs::remove_file(path);
        }

        assert!(
            worker_count >= 2,
            "expected multiple workers to process chunks, got {}",
            worker_count
        );
    }
}
