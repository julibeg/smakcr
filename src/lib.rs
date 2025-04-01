use std::{
    fmt::Display,
    fs::File,
    io::{stdout, BufRead, BufReader, BufWriter, Read, Write},
    path::Path,
    str::from_utf8,
    sync::{
        atomic::{AtomicUsize, Ordering},
        Arc, Mutex,
    },
};

use anyhow::{anyhow, Context, Error, Result};
use bio::alphabets::dna;
use seq_io::{fasta, fastq, parallel};

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
        Err(e) => Err(e).context(format!("Failed reading FASTA file \"{}\"", path)),
    }?;

    Ok(reader)
}

/// Converts a numerical k-mer index back into the corresponding DNA sequence.
///
/// Given an index (representing a k-mer encoded as a `usize` where each base takes 2 bits) and the
/// k-mer length `k`, this function reconstructs the original k-mer sequence as a `Vec<u8>` using
/// 'A', 'C', 'G', 'T'.
pub fn index_to_kmer(index: usize, k: usize) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut idx = index;
    let mut kmer = vec![b'N'; k];
    for item in kmer.iter_mut() {
        // extract the last two bits
        *item = bases[idx & 3];
        idx >>= 2;
    }
    // we still need to reverse to get the correct order
    kmer.reverse();
    kmer
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

/// An enum abstracting over [seq_io](https://crates.io/crates/seq_io) FASTA and FASTQ readers.
///
/// This allows processing files in either format using a common interface.
pub enum FastxReader {
    Fasta(fasta::Reader<Box<dyn Read + Send>>),
    Fastq(fastq::Reader<Box<dyn Read + Send>>),
}

impl FastxReader {
    /// Creates a `FastxReader` instance by detecting the format of the input file
    /// (compression-agnostic).
    ///
    /// It reads the first byte of the file ('>' indicates FASTA, '@' indicates FASTQ) and returns
    /// an error if the file is empty or the format is unrecognized.
    pub fn from_file<P: AsRef<Path> + Display>(path: P) -> Result<Self> {
        let reader = get_compression_agnostic_reader(&path)?;

        // wrap the reader in a BufReader and fill it to allow peeking at the first byte.
        let mut buf_reader = BufReader::new(reader);
        let buffer = buf_reader
            .fill_buf()
            .context("Failed reading first byte from file")?;
        if buffer.is_empty() {
            return Err(anyhow::anyhow!("File is empty"));
        }
        // Copy the first byte without consuming the buffer.
        let first_byte = buffer[0];

        let reader = match first_byte {
            b'>' => FastxReader::Fasta(fasta::Reader::new(Box::new(buf_reader))),
            b'@' => FastxReader::Fastq(fastq::Reader::new(Box::new(buf_reader))),
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid FASTx file format (first byte is '{}' and neither '>' nor '@'",
                    first_byte as char
                ));
            }
        };

        Ok(reader)
    }
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

    /// Processes the input data to count k-mers.
    ///
    /// Dispatches the counting task to either the single-threaded or parallel implementation based
    /// on the variant of `Self`. The input must implement the `CountKmers` trait.
    pub fn count<T: CountKmers>(&mut self, input: T) -> Result<()> {
        match self {
            Self::SingleThreaded { k, counts, mask } => input.count_kmers(*k, counts, *mask)?,
            Self::MultiThreaded {
                k,
                per_thread_counts,
                mask,
                n_threads,
            } => input.count_kmers_parallel(*k, per_thread_counts, *mask, *n_threads)?,
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
                        for (i, &count) in thread_count.iter().enumerate() {
                            first_thread_counts[i] += count;
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
    pub fn write(self, outfile: Option<&String>, canonical: bool, write_zeros: bool) -> Result<()> {
        // extract `k` from `self`
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

        if canonical {
            // look up the reverse complement of the kmer and print the smaller of the two
            // (alongside the sum of the counts)
            for (index, &count) in counts.iter().enumerate() {
                let kmer = index_to_kmer(index, k);
                // TODO: implement a way to get the index of the reverse complement directly
                let revcomp = dna::revcomp(&kmer);
                // only write if this kmer is canonical (i.e. smaller than the reverse complement)
                if kmer > revcomp {
                    continue;
                }
                // for a palindromic kmer (identical to the reverse complement), write the counts;
                // otherwise, write the sum of both counts
                let combined_count = if kmer == revcomp {
                    // kmer is palindrome
                    count
                } else {
                    // kmer is smaller than its reverse complement
                    count + counts[kmer_to_index(&revcomp)]
                };
                if combined_count == 0 && !write_zeros {
                    continue;
                }
                out_writer
                    .write_all(format!("{}\t{}\n", from_utf8(&kmer)?, combined_count).as_bytes())?;
            }
        } else {
            // also write non-canonical kmers
            for (index, &count) in counts.iter().enumerate() {
                if count == 0 && !write_zeros {
                    continue;
                };
                let kmer = index_to_kmer(index, k);
                out_writer.write_all(format!("{}\t{}\n", from_utf8(&kmer)?, count).as_bytes())?;
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
                if (k != other_k) | (counts.len() != other_counts.len()) {
                    return Err(anyhow!("Cannot add different types of `KmerCounter`"));
                }
                for (i, &count) in other_counts.iter().enumerate() {
                    counts[i] += count;
                }
                Ok(KmerCounter::SingleThreaded { k, counts, mask })
            }
            _ => Err(anyhow!("Cannot add different types of KmerCounters")),
        }
    }
}

/// A helper trait to provide a common interface for `seq_io::fasta::RefRecord` and
/// `seq_io::fastq::RefRecord`.
///
/// The trait enforces a method for checking if a `RefRecord`'s sequence is long enough to contain a
/// k-mer and for iteratively process each base with some function. This avoids code duplication in
/// the `CountKmers` implementation for the two different record types.
trait PerBaseRefRecord {
    /// Checks if the sequence record is potentially too short to contain any k-mers of length `k`.
    ///
    /// Implementations should account for the specific structure of the record format (e.g.,
    /// newlines in FASTA) to provide an accurate check.
    fn too_short(&self, k: usize) -> bool;

    /// Iterates over each base in the sequence record, applying the provided closure `f`.
    ///
    /// This abstracts away the details of how sequence data is stored or accessed
    /// within different record types (e.g., handling line breaks in FASTA).
    fn for_each_base<F>(&self, f: F)
    where
        F: FnMut(u8);
}

impl PerBaseRefRecord for fasta::RefRecord<'_> {
    /// Checks if the FASTA record is potentially too short for k-mers of length `k`.
    ///
    /// Considers potential newlines within the sequence data by checking `full_seq()`
    /// if the raw sequence length `seq().len()` seems borderline.
    fn too_short(&self, k: usize) -> bool {
        // FASTA-specific length check
        if <Self as fasta::Record>::seq(self).len() < k * 2 {
            // the sequence might be too short to contain a kmer. however, when checking we need to
            // make sure to take newlines into account (at worst we could be dealing with a FASTA
            // file with just one base per line and then `full_seq()` would be half the length of
            // `seq()`). note that we don't check `full_seq().len()` right away because
            // `full_seq()` allocates a new vector and we want to avoid that unless we have to
            self.full_seq().len() < k
        } else {
            false
        }
    }

    /// Iterates over bases in a FASTA record, handling potential line breaks.
    fn for_each_base<F>(&self, mut f: F)
    where
        F: FnMut(u8),
    {
        for line in self.seq_lines() {
            for &b in line {
                f(b);
            }
        }
    }
}

impl PerBaseRefRecord for fastq::RefRecord<'_> {
    /// Checks if the FASTQ record is too short for k-mers of length `k`.
    ///
    /// FASTQ sequences are contiguous, so this is a simple length check.
    fn too_short(&self, k: usize) -> bool {
        // FASTQ-specific length check
        <Self as fastq::Record>::seq(self).len() < k
    }

    /// Iterates over bases in a FASTQ record.
    fn for_each_base<F>(&self, mut f: F)
    where
        F: FnMut(u8),
    {
        for &b in <Self as fastq::Record>::seq(self) {
            f(b);
        }
    }
}

/// A trait defining the interface for k-mer counting operations.
///
/// This trait allows different types (like individual sequence records or entire
/// file readers) to implement their own k-mer counting logic, supporting both
/// single-threaded and parallel execution.
pub trait CountKmers {
    // NOTE: a simple way of allowing the trait on both mutable (`FastxReader`) and immutable
    // (`fast{a,q}::RefRecord`) data is to just pass `self` which consumes the object (as we're
    // doing now). this is fine because we're only passing over the data once anyway, but we might
    // still want to do something more sophisticated in the future

    /// Counts k-mers within the implementing type using a single thread.
    ///
    /// Takes the k-mer length `k`, a mutable slice `counts` to store the results,
    /// and the `mask` for rolling hash calculation. Updates the `counts` slice directly.
    fn count_kmers(self, k: usize, counts: &mut [usize], mask: usize) -> Result<()>;

    /// Counts k-mers within the implementing type using multiple threads.
    ///
    /// Takes the k-mer length `k`, a vector of thread-local count vectors `per_thread_counts`,
    /// the `mask` for rolling hash calculation, and the number of threads `n_threads`.
    /// Implementations should distribute the work across threads and update the
    /// corresponding count vectors in `per_thread_counts`.
    ///
    /// Provides a default single-threaded implementation that delegates to `count_kmers`.
    /// Types supporting parallelism should override this method.
    fn count_kmers_parallel(
        self,
        k: usize,
        per_thread_counts: &mut PerThreadVec<usize>,
        mask: usize,
        _n_threads: usize,
    ) -> Result<()>
    where
        Self: Sized,
    {
        // this provides a single-threaded blanket implementation; types that support parallel
        // counting can override it
        self.count_kmers(k, per_thread_counts[0].lock().unwrap().as_mut(), mask)?;

        Ok(())
    }
}

impl<T: PerBaseRefRecord> CountKmers for T {
    /// Implements single-threaded k-mer counting for types that provide base-level access
    /// via the `PerBaseRefRecord` trait (e.g., `fasta::RefRecord`, `fastq::RefRecord`).
    ///
    /// Iterates through the sequence bases, calculates the rolling k-mer index,
    /// and increments the corresponding count in the `counts` slice. Handles invalid
    /// bases ('N' or others not in `KEY_MAP`) by skipping affected k-mers.
    fn count_kmers(self, k: usize, counts: &mut [usize], mask: usize) -> Result<()>
    where
        Self: Sized,
    {
        // we first check if the record is too short to contain a single k-mer
        if self.too_short(k) {
            return Ok(());
        }

        // keep a skip counter to avoid counting kmers with unknown bases (and the first few
        // nucleotides that don't form a complete kmer)
        let mut skip = k - 1;
        let mut idx = 0;

        // make use of `for_each_base()` to iterate over all bases in the sequence
        self.for_each_base(|b| {
            idx = (idx << 2) & mask;
            let add = KEY_MAP[b as usize];
            if add == usize::MAX {
                // unknown base --> skip this kmer and the next `k - 1` kmers
                skip = k - 1;
                return;
            }
            idx += add;
            if skip == 0 {
                // only count the kmer if not skipping
                counts[idx] += 1;
            } else {
                skip -= 1;
            }
        });

        Ok(())
    }
}

impl CountKmers for FastxReader {
    /// Implements single-threaded k-mer counting for a `FastxReader`.
    ///
    /// Iterates through all records (FASTA or FASTQ) provided by the reader
    /// and delegates the counting for each record to its `count_kmers` implementation.
    fn count_kmers(self, k: usize, counts: &mut [usize], mask: usize) -> Result<()> {
        match self {
            Self::Fasta(mut reader) => {
                while let Some(record) = reader.next() {
                    record?.count_kmers(k, counts, mask)?;
                }
            }
            Self::Fastq(mut reader) => {
                while let Some(record) = reader.next() {
                    record?.count_kmers(k, counts, mask)?;
                }
            }
        }
        Ok(())
    }

    /// Implements parallel k-mer counting for a `FastxReader`.
    ///
    /// Uses `seq_io::parallel::read_parallel` to read records concurrently and
    /// distribute the counting task across multiple threads. Each thread updates
    /// its assigned count vector within `per_thread_counts`.
    fn count_kmers_parallel(
        self,
        k: usize,
        per_thread_counts: &mut PerThreadVec<usize>,
        mask: usize,
        n_threads: usize,
    ) -> Result<()> {
        let queue_len = n_threads * 2;

        let record_set_counter = Arc::new(AtomicUsize::new(0));

        // the match arms below would only differ by the type of the `record_set` parameter, so we
        // can use a macro to avoid code duplication
        macro_rules! parallel_read_parallel_impl {
            ($reader:expr, $record_set_type:ty) => {
                parallel::read_parallel(
                    $reader,
                    // use `n_threads - 1` because `parallel::read_parallel` also spawns a reader
                    // thread
                    n_threads as u32 - 1,
                    queue_len,
                    |record_set: &mut $record_set_type| {
                        // get the thread ID and the corresponding counts vector
                        let thread_id =
                            record_set_counter.fetch_add(1, Ordering::SeqCst) % n_threads;

                        let mut counts = per_thread_counts[thread_id].lock().unwrap();

                        for record in record_set.into_iter() {
                            record
                                .count_kmers(k, &mut counts, mask)
                                .expect("Error counting bases in `read_parallel`");
                        }
                    },
                    |record_sets| {
                        // nothing to do, just error propagation
                        while let Some(result) = record_sets.next() {
                            result?;
                        }
                        Ok::<(), Error>(())
                    },
                )?;
            };
        }

        match self {
            Self::Fasta(reader) => {
                parallel_read_parallel_impl!(reader, fasta::RecordSet);
            }
            Self::Fastq(reader) => {
                parallel_read_parallel_impl!(reader, fastq::RecordSet);
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // NOTE: This impl is just for testing convenience, allowing &[u8] to be used
    // where `PerBaseRefRecord` is expected in tests.
    impl PerBaseRefRecord for &[u8] {
        fn too_short(&self, k: usize) -> bool {
            self.len() < k
        }

        fn for_each_base<F>(&self, mut f: F)
        where
            F: FnMut(u8),
        {
            for &b in *self {
                f(b);
            }
        }
    }

    // NOTE: Helper function for tests.
    fn count_kmers_in_sequence(seq: &[u8], k: usize) -> Vec<usize> {
        let mut counter = KmerCounter::new(k, 1);
        // NOTE: Using the test impl of `PerBaseRefRecord` for `&[u8]` here.
        counter.count(seq).unwrap();
        counter.into_vec().unwrap()
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
        // TODO: this actually only tests the implementation for `&[u8]` and for the other types.
        // the test would be more useful if we only used a default implementation of `count_kmers()`
        // that relies on other functions (e.g. `too_short()` and `process_base()`)
        let k = 5;
        let counts = count_kmers_in_sequence(b"ACG", k);

        assert!(counts.iter().all(|&count| count == 0));
    }

    #[test]
    fn test_index_to_kmer() {
        let k = 3;
        let kmer = index_to_kmer(0b00_01_10, k);
        assert_eq!(kmer, b"ACG");

        let acg_index = 0b00_01_10; // ACG
        let cgt_index = 0b01_10_11; // CGT
        let gta_index = 0b10_11_00; // GTA
        let tac_index = 0b11_00_01; // TAC
        let cga_index = 0b01_10_00; // CGA
        let gaa_index = 0b10_00_00; // GAA
        let aaa_index = 0b00_00_00; // AAA

        assert_eq!(index_to_kmer(acg_index, k), b"ACG");
        assert_eq!(index_to_kmer(cgt_index, k), b"CGT");
        assert_eq!(index_to_kmer(gta_index, k), b"GTA");
        assert_eq!(index_to_kmer(tac_index, k), b"TAC");
        assert_eq!(index_to_kmer(cga_index, k), b"CGA");
        assert_eq!(index_to_kmer(gaa_index, k), b"GAA");
        assert_eq!(index_to_kmer(aaa_index, k), b"AAA");
    }

    #[test]
    fn test_unknwn_bases() {
        let k = 3;
        let counts = count_kmers_in_sequence(b"AaAxcGtx", k);

        let aaa_index = 0b00_00_00; // AAA
        let cgt_index = 0b01_10_11; // CGT
        assert_eq!(counts[aaa_index], 1);
        assert_eq!(counts[cgt_index], 1);
        // make sure all other counts are 0
        assert_eq!(counts.iter().sum::<usize>(), 2);
    }

    // TODO: add tests for the whole program (i.e. against an existing file to make sure the output
    //   format etc is as expected)
}
