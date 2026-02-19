//! Parser benchmark: compare base-counting throughput across FASTX parsers.
//!
//! Parsers compared:
//!   - **smakcr (FastxChunker)** — custom memchr-accelerated chunking parser
//!   - **seq_io** — popular zero-copy FASTA/FASTQ parser
//!   - **needletail** — another popular Rust FASTX parser
//!
//! Workload: count total number of valid DNA bases (A/C/G/T, case-insensitive)
//! per file. This isolates parser overhead from any k-mer counting logic.
//!
//! Test scenarios (set via BENCH_FILE env var):
//!   - Single long sequence  (e.g. chrM, chr1)
//!   - Many sequences in one file (e.g. chrUn_combined)
//!   - Default: chrM
//!
//! Run with:
//!   cargo bench --bench parser                               # default: chrM
//!   BENCH_FILE=tests/fixtures/chrUn_combined.fa.gz cargo bench --bench parser
//!   BENCH_FILE=benchmark/human/per-chromosome/chr1.fa.gz cargo bench --bench parser
//!
//! Filter specific benchmarks:
//!   cargo bench --bench parser -- smakcr
//!   cargo bench --bench parser -- seq_io
//!   cargo bench --bench parser -- needletail
//!   cargo bench --bench parser -- chunk_size

use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use flate2::read::GzDecoder;

const DEFAULT_FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/chrM.fa.gz");

/// Load a FASTA/FASTQ file into memory, decompressing gzip if needed.
fn load_file(path: &str) -> Vec<u8> {
    let file = File::open(path).unwrap_or_else(|e| panic!("cannot open {}: {}", path, e));
    let mut buf_reader = BufReader::new(file);

    // Peek at first two bytes — gzip magic number is 0x1f 0x8b
    let header = buf_reader.fill_buf().unwrap();
    let is_gzip = header.len() >= 2 && header[0] == 0x1f && header[1] == 0x8b;

    let mut data = Vec::new();
    if is_gzip {
        GzDecoder::new(buf_reader)
            .read_to_end(&mut data)
            .unwrap();
    } else {
        buf_reader.read_to_end(&mut data).unwrap();
    }
    data
}

fn file_stem(path: &str) -> String {
    let name = Path::new(path)
        .file_name()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    name.strip_suffix(".fa")
        .or_else(|| name.strip_suffix(".fasta"))
        .or_else(|| name.strip_suffix(".fna"))
        .or_else(|| name.strip_suffix(".fq"))
        .or_else(|| name.strip_suffix(".fastq"))
        .unwrap_or(name)
        .to_string()
}

fn human_bytes(n: u64) -> String {
    if n >= 1 << 30 {
        format!("{:.1} GiB", n as f64 / (1u64 << 30) as f64)
    } else if n >= 1 << 20 {
        format!("{:.1} MiB", n as f64 / (1u64 << 20) as f64)
    } else if n >= 1 << 10 {
        format!("{:.1} KiB", n as f64 / (1u64 << 10) as f64)
    } else {
        format!("{} B", n)
    }
}

// ---------------------------------------------------------------------------
// Parser implementations — each returns total valid DNA bases
// ---------------------------------------------------------------------------

fn count_bases_smakcr(data: &[u8], chunk_size: usize) -> u64 {
    use smakcr::KEY_MAP;

    let reader = BufReader::new(std::io::Cursor::new(data));
    let mut chunker =
        smakcr::chunker_for_bench::FastxChunker::new(reader, 1, chunk_size).unwrap();

    let mut total: u64 = 0;
    while let Some(chunk) = chunker.next_chunk().unwrap() {
        for &b in &chunk {
            if KEY_MAP[b as usize] != usize::MAX {
                total += 1;
            }
        }
    }
    total
}

fn count_bases_seq_io(data: &[u8]) -> u64 {
    use smakcr::KEY_MAP;

    let mut reader = seq_io::fasta::Reader::new(std::io::Cursor::new(data));
    let mut total: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record.unwrap();
        for seq_line in record.seq_lines() {
            for &b in seq_line {
                if KEY_MAP[b as usize] != usize::MAX {
                    total += 1;
                }
            }
        }
    }
    total
}

fn count_bases_needletail(data: &[u8]) -> u64 {
    use smakcr::KEY_MAP;

    let mut reader = needletail::parse_fastx_reader(std::io::Cursor::new(data)).unwrap();
    let mut total: u64 = 0;

    while let Some(record) = reader.next() {
        let record = record.unwrap();
        let seq = record.seq();
        for &b in seq.iter() {
            if KEY_MAP[b as usize] != usize::MAX {
                total += 1;
            }
        }
    }
    total
}

// ---------------------------------------------------------------------------
// Tuning helpers
// ---------------------------------------------------------------------------

/// Configure group for large files: fewer samples, longer measurement.
fn tune_for_size(group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>, bytes: u64) {
    if bytes > 500 << 20 {
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(30));
        group.warm_up_time(Duration::from_secs(5));
    } else if bytes > 50 << 20 {
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(15));
    } else if bytes > 5 << 20 {
        group.sample_size(20);
        group.measurement_time(Duration::from_secs(10));
    }
}

// ---------------------------------------------------------------------------
// Benchmark groups
// ---------------------------------------------------------------------------

fn bench_parsers(c: &mut Criterion) {
    let path = std::env::var("BENCH_FILE").unwrap_or_else(|_| DEFAULT_FILE.to_string());
    let label = file_stem(&path);

    eprintln!("Loading {}...", path);
    let data = load_file(&path);
    let data_len = data.len() as u64;
    eprintln!("Loaded {} ({} decompressed)", label, human_bytes(data_len));

    // --- Accuracy check ---
    eprintln!("Verifying parser accuracy...");
    let count_smakcr = count_bases_smakcr(&data, 1 << 20);
    let count_seqio = count_bases_seq_io(&data);
    let count_needletail = count_bases_needletail(&data);

    eprintln!("  smakcr:     {} bases", count_smakcr);
    eprintln!("  seq_io:     {} bases", count_seqio);
    eprintln!("  needletail: {} bases", count_needletail);

    assert_eq!(
        count_smakcr, count_seqio,
        "MISMATCH: smakcr ({}) != seq_io ({})",
        count_smakcr, count_seqio,
    );
    assert_eq!(
        count_smakcr, count_needletail,
        "MISMATCH: smakcr ({}) != needletail ({})",
        count_smakcr, count_needletail,
    );
    eprintln!("  All parsers agree on {} bases", count_smakcr);

    // --- Throughput benchmarks ---
    let mut group = c.benchmark_group(format!("parser/{}", label));
    group.throughput(Throughput::Bytes(data_len));
    tune_for_size(&mut group, data_len);

    group.bench_function("smakcr", |b| {
        b.iter(|| count_bases_smakcr(&data, 1 << 20))
    });

    group.bench_function("seq_io", |b| b.iter(|| count_bases_seq_io(&data)));

    group.bench_function("needletail", |b| {
        b.iter(|| count_bases_needletail(&data))
    });

    group.finish();
}

fn bench_chunk_sizes(c: &mut Criterion) {
    let path = std::env::var("BENCH_FILE").unwrap_or_else(|_| DEFAULT_FILE.to_string());
    let label = file_stem(&path);
    let data = load_file(&path);
    let data_len = data.len() as u64;

    let mut group = c.benchmark_group(format!("parser/{}/smakcr_chunk_size", label));
    group.throughput(Throughput::Bytes(data_len));
    tune_for_size(&mut group, data_len);

    for chunk_exp in [12, 14, 16, 18, 20] {
        let chunk_size = 1usize << chunk_exp;
        let label = if chunk_size >= 1024 {
            format!("{}K", chunk_size / 1024)
        } else {
            format!("{}B", chunk_size)
        };
        group.bench_with_input(
            BenchmarkId::from_parameter(label),
            &chunk_size,
            |b, &cs| b.iter(|| count_bases_smakcr(&data, cs)),
        );
    }

    group.finish();
}

criterion_group!(benches, bench_parsers, bench_chunk_sizes);
criterion_main!(benches);
