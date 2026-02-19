use std::io::BufRead;

use anyhow::{anyhow, Context, Result};
use memchr::memchr;

pub(crate) const SEQ_SENTINEL: u8 = 0xFF;

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum FastxFormat {
    Fasta,
    Fastq,
}

fn detect_fastx_format<R: BufRead>(reader: &mut R) -> Result<FastxFormat> {
    let buffer = reader
        .fill_buf()
        .context("Failed reading first byte from file")?;
    if buffer.is_empty() {
        return Err(anyhow!("File is empty"));
    }
    match buffer[0] {
        b'>' => Ok(FastxFormat::Fasta),
        b'@' => Ok(FastxFormat::Fastq),
        _ => Err(anyhow!(
            "Invalid FASTx file format (first byte is '{}' and neither '>' nor '@'",
            buffer[0] as char
        )),
    }
}

/// Append `span` to `out`, stripping newlines (\n and \r) using memchr for bulk copies.
#[inline]
fn append_stripping_newlines(out: &mut Vec<u8>, span: &[u8]) {
    let mut pos = 0;
    while pos < span.len() {
        // Find next \n or \r
        match memchr::memchr2(b'\n', b'\r', &span[pos..]) {
            Some(offset) => {
                if offset > 0 {
                    out.extend_from_slice(&span[pos..pos + offset]);
                }
                pos += offset + 1;
            }
            None => {
                out.extend_from_slice(&span[pos..]);
                break;
            }
        }
    }
}

/// FASTA parser state
struct FastaParser {
    /// True once we've seen the first '>' header
    seen_header: bool,
    /// True if current sequence has had any bases emitted
    seq_has_bases: bool,
}

impl FastaParser {
    fn new() -> Self {
        Self {
            seen_header: false,
            seq_has_bases: false,
        }
    }

    /// Fill `chunk` from the reader. Returns true if there is more data, false at EOF.
    fn fill_chunk<R: BufRead>(
        &mut self,
        reader: &mut R,
        chunk: &mut Vec<u8>,
        limit: usize,
    ) -> Result<bool> {
        loop {
            if chunk.len() >= limit {
                return Ok(true);
            }

            let buf = reader.fill_buf().context("FASTA read error")?;
            if buf.is_empty() {
                return Ok(false);
            }

            let buf_len = buf.len();
            let mut consumed = 0;

            while consumed < buf_len && chunk.len() < limit {
                if !self.seen_header {
                    // Skip to end of first header line
                    match memchr(b'\n', &buf[consumed..]) {
                        Some(offset) => {
                            consumed += offset + 1;
                            self.seen_header = true;
                        }
                        None => {
                            consumed = buf_len;
                        }
                    }
                    continue;
                }

                // Look for next '>' (new record) in remaining buffer
                match memchr(b'>', &buf[consumed..]) {
                    Some(offset) => {
                        // Everything before '>' is sequence data from current record
                        if offset > 0 {
                            append_stripping_newlines(chunk, &buf[consumed..consumed + offset]);
                            if !self.seq_has_bases {
                                self.seq_has_bases = true;
                            }
                        }
                        consumed += offset; // now pointing at '>'

                        // Emit sentinel between sequences
                        if self.seq_has_bases {
                            chunk.push(SEQ_SENTINEL);
                        }
                        self.seq_has_bases = false;

                        // Skip header line (from '>' to '\n')
                        match memchr(b'\n', &buf[consumed..]) {
                            Some(nl_offset) => {
                                consumed += nl_offset + 1;
                                self.seen_header = true;
                            }
                            None => {
                                // Header spans beyond this buffer; consume what we have,
                                // mark that we need to skip rest of header
                                consumed = buf_len;
                                self.seen_header = false;
                            }
                        }
                    }
                    None => {
                        // No '>' found — rest of buffer is sequence data
                        append_stripping_newlines(chunk, &buf[consumed..]);
                        if !self.seq_has_bases {
                            self.seq_has_bases = true;
                        }
                        consumed = buf_len;
                    }
                }
            }

            reader.consume(consumed);
        }
    }
}

/// FASTQ parser state — 4-line state machine
struct FastqParser {
    /// 0=header, 1=sequence, 2=plus, 3=quality
    line_state: u8,
}

impl FastqParser {
    fn new() -> Self {
        Self { line_state: 0 }
    }

    fn fill_chunk<R: BufRead>(
        &mut self,
        reader: &mut R,
        chunk: &mut Vec<u8>,
        limit: usize,
    ) -> Result<bool> {
        loop {
            if chunk.len() >= limit {
                return Ok(true);
            }

            let buf = reader.fill_buf().context("FASTQ read error")?;
            if buf.is_empty() {
                return Ok(false);
            }

            let buf_len = buf.len();
            let mut consumed = 0;

            while consumed < buf_len && chunk.len() < limit {
                match self.line_state {
                    0 => {
                        // Header line — skip to newline
                        match memchr(b'\n', &buf[consumed..]) {
                            Some(offset) => {
                                consumed += offset + 1;
                                self.line_state = 1;
                            }
                            None => {
                                consumed = buf_len;
                            }
                        }
                    }
                    1 => {
                        // Sequence line — copy bases, append sentinel at end
                        match memchr(b'\n', &buf[consumed..]) {
                            Some(offset) => {
                                let line = &buf[consumed..consumed + offset];
                                // Strip \r if present
                                let line = if line.last() == Some(&b'\r') {
                                    &line[..line.len() - 1]
                                } else {
                                    line
                                };
                                chunk.extend_from_slice(line);
                                chunk.push(SEQ_SENTINEL);
                                consumed += offset + 1;
                                self.line_state = 2;
                            }
                            None => {
                                // Sequence spans beyond buffer — copy what we have
                                let span = &buf[consumed..];
                                let span = if span.last() == Some(&b'\r') {
                                    &span[..span.len() - 1]
                                } else {
                                    span
                                };
                                chunk.extend_from_slice(span);
                                consumed = buf_len;
                            }
                        }
                    }
                    2 => {
                        // Plus line — skip
                        match memchr(b'\n', &buf[consumed..]) {
                            Some(offset) => {
                                consumed += offset + 1;
                                self.line_state = 3;
                            }
                            None => {
                                consumed = buf_len;
                            }
                        }
                    }
                    _ => {
                        // Quality line — skip
                        match memchr(b'\n', &buf[consumed..]) {
                            Some(offset) => {
                                consumed += offset + 1;
                                self.line_state = 0;
                            }
                            None => {
                                consumed = buf_len;
                            }
                        }
                    }
                }
            }

            reader.consume(consumed);
        }
    }
}

enum FastxStream<R: BufRead> {
    Fasta { reader: R, parser: FastaParser },
    Fastq { reader: R, parser: FastqParser },
}

impl<R: BufRead> FastxStream<R> {
    fn fill_chunk(&mut self, chunk: &mut Vec<u8>, limit: usize) -> Result<bool> {
        match self {
            Self::Fasta { reader, parser } => parser.fill_chunk(reader, chunk, limit),
            Self::Fastq { reader, parser } => parser.fill_chunk(reader, chunk, limit),
        }
    }
}

pub struct FastxChunker<R: BufRead> {
    stream: FastxStream<R>,
    k: usize,
    chunk_size: usize,
    tail: Vec<u8>,
    done: bool,
}

impl<R: BufRead> FastxChunker<R> {
    pub fn new(mut reader: R, k: usize, chunk_size: usize) -> Result<Self> {
        if chunk_size < k {
            return Err(anyhow!("Chunk size ({}) must be >= k ({})", chunk_size, k));
        }

        let format = detect_fastx_format(&mut reader)?;
        let stream = match format {
            FastxFormat::Fasta => FastxStream::Fasta {
                reader,
                parser: FastaParser::new(),
            },
            FastxFormat::Fastq => FastxStream::Fastq {
                reader,
                parser: FastqParser::new(),
            },
        };

        Ok(Self {
            stream,
            k,
            chunk_size,
            tail: Vec::new(),
            done: false,
        })
    }

    pub fn next_chunk(&mut self) -> Result<Option<Vec<u8>>> {
        if self.done {
            return Ok(None);
        }

        let mut chunk = Vec::with_capacity(self.chunk_size + self.k);

        // Prepend tail from previous chunk (overlap for k-mer continuity)
        if !self.tail.is_empty() {
            chunk.extend_from_slice(&self.tail);
        }

        let had_more = self.stream.fill_chunk(&mut chunk, self.chunk_size)?;

        if !had_more {
            self.done = true;
        }

        if chunk.is_empty() {
            return Ok(None);
        }

        // Extract tail for next chunk: last (k-1) bytes, but respect sentinels
        if self.k > 1 {
            // Find position of last sentinel in the chunk
            let tail_start = chunk.len().saturating_sub(self.k - 1);
            self.tail.clear();
            self.tail.extend_from_slice(&chunk[tail_start..]);
        }

        Ok(Some(chunk))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::count_kmers_bytes;
    use std::io::{BufReader, Cursor};

    fn count_via_chunks(input: &[u8], k: usize, chunk_size: usize) -> Vec<usize> {
        let reader = BufReader::new(Cursor::new(input));
        let mut chunker = FastxChunker::new(reader, k, chunk_size).unwrap();
        let mask = 4usize.pow(k as u32) - 1;
        let mut counts = vec![0; 4usize.pow(k as u32)];
        while let Some(chunk) = chunker.next_chunk().unwrap() {
            count_kmers_bytes(k, &mut counts, mask, &chunk, Some(SEQ_SENTINEL));
        }
        counts
    }

    fn count_all_at_once(input: &[u8], k: usize) -> Vec<usize> {
        // Use a very large chunk size so everything fits in one chunk
        count_via_chunks(input, k, 1_000_000)
    }

    #[test]
    fn test_fasta_basic_counts() {
        let input = b">seq1\nACGTAC\n>seq2\nTTAA\n";
        let k = 3;
        let counts = count_all_at_once(input, k);

        // Verify specific k-mers from seq1: ACG, CGT, GTA, TAC
        let acg = 0b00_01_10;
        let cgt = 0b01_10_11;
        let gta = 0b10_11_00;
        let tac = 0b11_00_01;
        assert_eq!(counts[acg], 1);
        assert_eq!(counts[cgt], 1);
        assert_eq!(counts[gta], 1);
        assert_eq!(counts[tac], 1);

        // seq2: TTA, TAA
        let tta = 0b11_11_00;
        let taa = 0b11_00_00;
        assert_eq!(counts[tta], 1);
        assert_eq!(counts[taa], 1);

        assert_eq!(counts.iter().sum::<usize>(), 6);
    }

    #[test]
    fn test_fasta_chunked_matches_single() {
        let input = b">seq1\nACGTACGTACGT\n>seq2\nTTAAGGCC\n";
        for k in 1..=5 {
            let expected = count_all_at_once(input, k);
            // Try small chunk sizes to exercise overlap logic
            for chunk_size in [k, k + 1, 5, 8, 13, 100] {
                let got = count_via_chunks(input, k, chunk_size);
                assert_eq!(
                    expected, got,
                    "mismatch for k={}, chunk_size={}",
                    k, chunk_size
                );
            }
        }
    }

    #[test]
    fn test_fasta_multiline_sequence() {
        let input = b">seq\nACGT\nACGT\n";
        let k = 3;
        let counts = count_all_at_once(input, k);

        // Should be equivalent to "ACGTACGT" — 6 3-mers
        assert_eq!(counts.iter().sum::<usize>(), 6);

        // The junction k-mer "TAC" should be present
        let tac = 0b11_00_01;
        assert_eq!(counts[tac], 1);
    }

    #[test]
    fn test_fasta_multiline_chunked() {
        let input = b">seq\nACGT\nACGT\nACGT\n";
        for k in 1..=4 {
            let expected = count_all_at_once(input, k);
            for chunk_size in [k, k + 2, 7, 20] {
                let got = count_via_chunks(input, k, chunk_size);
                assert_eq!(
                    expected, got,
                    "multiline mismatch k={}, chunk_size={}",
                    k, chunk_size
                );
            }
        }
    }

    #[test]
    fn test_fastq_basic_counts() {
        let input = b"@r1\nACGTA\n+\n!!!!!\n@r2\nTT\n+\n!!\n";
        let k = 3;
        let counts = count_all_at_once(input, k);

        // r1: ACG, CGT, GTA = 3 k-mers
        // r2: too short (only TT, length 2 < k=3)
        assert_eq!(counts.iter().sum::<usize>(), 3);
    }

    #[test]
    fn test_fastq_chunked_matches_single() {
        let input = b"@r1\nACGTACGT\n+\n!!!!!!!!\n@r2\nTTAAGGCC\n+\n!!!!!!!!\n";
        for k in 1..=5 {
            let expected = count_all_at_once(input, k);
            for chunk_size in [k, k + 1, 5, 8, 20] {
                let got = count_via_chunks(input, k, chunk_size);
                assert_eq!(
                    expected, got,
                    "FASTQ mismatch k={}, chunk_size={}",
                    k, chunk_size
                );
            }
        }
    }

    #[test]
    fn test_format_detection_fasta() {
        let input = b">seq\nACGT\n";
        let mut reader = BufReader::new(Cursor::new(input));
        assert_eq!(
            detect_fastx_format(&mut reader).unwrap(),
            FastxFormat::Fasta
        );
    }

    #[test]
    fn test_format_detection_fastq() {
        let input = b"@read\nACGT\n+\n!!!!\n";
        let mut reader = BufReader::new(Cursor::new(input));
        assert_eq!(
            detect_fastx_format(&mut reader).unwrap(),
            FastxFormat::Fastq
        );
    }

    #[test]
    fn test_format_detection_empty() {
        let input = b"";
        let mut reader = BufReader::new(Cursor::new(input));
        assert!(detect_fastx_format(&mut reader).is_err());
    }

    #[test]
    fn test_empty_sequences() {
        let input = b">seq1\n>seq2\nACGT\n";
        let k = 3;
        let counts = count_all_at_once(input, k);
        // Only seq2 has bases: ACG, CGT
        assert_eq!(counts.iter().sum::<usize>(), 2);
    }

    #[test]
    fn test_single_base_kmer() {
        let input = b">seq\nACGT\n";
        let counts = count_all_at_once(input, 1);
        // A=1, C=1, G=1, T=1
        assert_eq!(counts[0], 1); // A
        assert_eq!(counts[1], 1); // C
        assert_eq!(counts[2], 1); // G
        assert_eq!(counts[3], 1); // T
    }

    #[test]
    fn test_sequence_shorter_than_k() {
        let input = b">seq\nAC\n";
        let k = 5;
        let counts = count_all_at_once(input, k);
        assert_eq!(counts.iter().sum::<usize>(), 0);
    }
}
