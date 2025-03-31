use anyhow::{Context, Result};
use bio::io::fasta;
use std::{
    path::Path,
    fmt::Display,
    fs::File,
    io::{BufRead, BufReader, Read},
};

pub fn read_fasta<P: AsRef<Path> + Display>(
    path: P,
) -> Result<impl Iterator<Item = Result<bio::io::fasta::Record>>> {
    let reader = get_compression_agnostic_reader(&path)?;
    let fasta_reader = fasta::Reader::new(reader);

    // add a more specific error context to each record
    let records = fasta_reader
        .records()
        .map(|item| item.context("Could not parse FASTA record"));

    Ok(records)
}

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

fn _print_lines(data: &[u8]) {
    process_lines(data, &mut |line| {
        println!("{}", String::from_utf8_lossy(line));
    });
}

fn _process_lines_2(data: &[u8], f: &mut dyn FnMut(&[u8])) {
    for line in data.split(|&c| c == b'\n') {
        f(line);
    }
}

fn process_lines(data: &[u8], f: &mut dyn FnMut(&[u8])) {
    let mut start = 0;
    let mut line: &[u8];
    while start < data.len() {
        if let Some(pos) = memchr::memchr(b'\n', &data[start..]) {
            line = &data[start..start + pos];
            start += pos + 1;
        } else {
            // no more newlines; we're at the end
            line = &data[start..];
            start = data.len();
        }
        f(line);
    }
}

#[derive(Debug)]
pub struct Chunk {
    data: Vec<u8>,
    seq_start: usize,
    seq_end: usize,
    no_seqs: bool,
    starts_new_seq: bool,
}

impl Chunk {
    pub fn print_seq_data(&self) {
        self.process_seq_data(
            &mut |line| {
                println!("{}", String::from_utf8_lossy(line));
            },
            &mut || {
                println!("==");
            },
        );
    }

    pub fn process_seq_data(
        &self,
        line_func: &mut dyn FnMut(&[u8]),
        seq_start_func: &mut dyn FnMut(),
    ) {
        // prints sequence lines, ignoring headers
        if self.no_seqs {
            return;
        }

        let data = &self.data[self.seq_start..self.seq_end];

        let mut start = 0;

        if self.starts_new_seq {
            seq_start_func();
        }
        // iterate over sequences by looking for '>'
        while start < data.len() {
            if let Some(header_start) = memchr::memchr(b'>', &data[start..]) {
                // we found a header line, let's see if we got the whole header in the chunk and
                // skip it if so
                if let Some(header_end) = memchr::memchr(b'\n', &data[start + header_start..]) {
                    if header_start > 0 {
                        // we got a sequence before the header
                        let seq = &data[start..start + header_start - 1]; // -1 to skip the newline
                        process_lines(seq, line_func);
                    }
                    start = start + header_start + header_end + 1;
                    seq_start_func();
                } else {
                    // the rest of the chunk is in the header line, but we actually should have
                    // skipped any header at the end -> panic
                    panic!(
                        "Found a header line at the end of the chunk; {}",
                        String::from_utf8_lossy(&data[start..])
                    );
                }
            } else {
                // the rest of the chunk is part of a sequence
                process_lines(&data[start..], line_func);
                start = data.len();
            }
        }
    }
}

pub struct ChunkReader<R: Read> {
    reader: R,
    chunk_size: usize,
    prev_chunk_ended_in_header: bool,
    eof: bool,
}

impl<R: Read> ChunkReader<R> {
    pub fn new(reader: R, chunk_size: usize) -> Self {
        ChunkReader {
            reader,
            chunk_size,
            prev_chunk_ended_in_header: true, // Assume we start with a header
            eof: false,
        }
    }
}

impl ChunkReader<Box<dyn Read>> {
    pub fn from_file<P: AsRef<Path> + Display>(path: P, chunk_size: usize) -> Result<Self> {
        let reader = get_compression_agnostic_reader(&path)?;
        Ok(ChunkReader::<Box<dyn Read>>::new(reader, chunk_size))
    }
}

impl<R: Read> Iterator for ChunkReader<R> {
    type Item = Result<Chunk>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.eof {
            return None;
        }

        // TODO: we might want to avoid initialising this if profiling shows that it's indeed a
        // bottleneck (but let's keep it for now)
        let mut chunk_buffer = vec![0u8; self.chunk_size];

        // Read up to chunk_size bytes
        let bytes_read = match self.reader.read(&mut chunk_buffer) {
            Ok(0) => {
                self.eof = true;
                return None;
            }
            Ok(n) => n,
            Err(e) => {
                return Some(Err(
                    anyhow::Error::new(e).context("Failed to read chunk from file")
                ))
            }
        };

        // Resize the buffer to the actual number of bytes read
        chunk_buffer.truncate(bytes_read);

        let mut seq_start = 0;
        let mut seq_end = chunk_buffer.len();
        let mut no_seqs = false;
        let mut starts_new_seq = false;

        // check if the chunk starts with or within a header
        if self.prev_chunk_ended_in_header || chunk_buffer[0] == b'>' {
            if let Some(pos) = memchr::memchr(b'\n', &chunk_buffer) {
                // we're starting in a header and it ends before the end of the chunk
                seq_start = pos + 1;
                starts_new_seq = true;
            } else {
                // no newline found; the whole chunk is within the header line
                seq_end = 0;
                no_seqs = true;
                self.prev_chunk_ended_in_header = true;
            }
        }

        // to see if we end in a header line, search from the end of the chunk to find the last
        // newline and check if there's a '>' right after
        if !no_seqs {
            if let Some(last_nl_pos) = memchr::memrchr(b'\n', &chunk_buffer) {
                if last_nl_pos + 1 < chunk_buffer.len() && chunk_buffer[last_nl_pos + 1] == b'>' {
                    seq_end = last_nl_pos;
                    self.prev_chunk_ended_in_header = true;
                } else {
                    self.prev_chunk_ended_in_header = false;
                }
            }
        }

        let chunk = Chunk {
            data: chunk_buffer,
            seq_start,
            seq_end,
            starts_new_seq,
            no_seqs,
        };

        Some(Ok(chunk))
    }
}

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

pub fn kmer_to_index(kmer: &Vec<u8>, key_map: [usize; 256]) -> usize {
    let mut idx = 0;
    for b in kmer {
        idx = (idx << 2) + key_map[*b as usize];
    }
    idx
}

pub enum FastxReader {
    Fasta(seq_io::fasta::Reader<Box<dyn Read + Send>>),
    Fastq(seq_io::fastq::Reader<Box<dyn Read + Send>>),
}

impl FastxReader {
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
            b'>' => FastxReader::Fasta(seq_io::fasta::Reader::new(Box::new(buf_reader))),
            b'@' => FastxReader::Fastq(seq_io::fastq::Reader::new(Box::new(buf_reader))),
            _ => {
                return Err(anyhow::anyhow!(
                    "Invalid FASTx file format (first byte is {} and neither '>' nor '@'",
                    first_byte as char
                ));
            }
        };

        Ok(reader)
    }
}

