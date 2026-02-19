//! End-to-end k-mer counting benchmark across file layouts and thread counts.
//!
//! Scenarios:
//!   1. **single_seq**  — one file, one long sequence (e.g. chrM or chr1)
//!   2. **multi_seq**   — one file, many sequences (chrUn_combined)
//!   3. **multi_file**  — many files, one sequence each (chrUn per-scaffold)
//!
//! Each scenario is tested at 1, 4, and 8 threads (configurable via BENCH_THREADS).
//! Workload: full 5-mer counting via `KmerCounter::count_paths`.
//!
//! Run with:
//!   cargo bench --bench e2e
//!   cargo bench --bench e2e -- single_seq       # only the single-sequence scenario
//!   cargo bench --bench e2e -- multi_file/t4    # only multi-file at 4 threads
//!
//! Override the single-sequence file (default chrM):
//!   BENCH_SINGLE_SEQ=benchmark/human/per-chromosome/chr1.fa.gz cargo bench --bench e2e

use std::path::Path;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

// ---------------------------------------------------------------------------
// Paths
// ---------------------------------------------------------------------------

const CHRM_PATH: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/chrM.fa.gz");

const CHRUN_COMBINED: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/chrUn_combined.fa.gz"
);

const CHRUN_PER_FILE_DIR: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/tests/fixtures/chrUn-one-seq-per-file"
);

fn thread_counts() -> Vec<usize> {
    match std::env::var("BENCH_THREADS") {
        Ok(val) => val.split(',').map(|s| s.trim().parse().unwrap()).collect(),
        Err(_) => vec![1, 4, 8],
    }
}

/// Collect all .fa.gz files in a directory, sorted for determinism.
fn collect_files(dir: &str) -> Vec<String> {
    let mut files: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .filter_map(|e| {
            let e = e.unwrap();
            let name = e.file_name().to_str().unwrap().to_string();
            if name.ends_with(".fa.gz") || name.ends_with(".fa") {
                Some(e.path().to_str().unwrap().to_string())
            } else {
                None
            }
        })
        .collect();
    files.sort();
    files
}

/// Total bytes of all files (compressed, for throughput baseline).
fn total_file_bytes(paths: &[String]) -> u64 {
    paths
        .iter()
        .map(|p| std::fs::metadata(p).unwrap().len())
        .sum()
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

/// Run smakcr k-mer counting end-to-end. Returns the total count (for verification).
fn run_smakcr(paths: &[String], k: usize, threads: usize, chunk_size: usize) -> u64 {
    let mut counter = smakcr::KmerCounter::new(k, threads);
    counter
        .count_paths(paths, chunk_size, 64 * 1024 * 1024)
        .unwrap();
    let counts = counter.into_vec().unwrap();
    counts.iter().sum::<usize>() as u64
}

/// Configure group for large files.
fn tune_for_size(
    group: &mut criterion::BenchmarkGroup<criterion::measurement::WallTime>,
    bytes: u64,
) {
    if bytes > 500 << 20 {
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(60));
        group.warm_up_time(Duration::from_secs(10));
    } else if bytes > 50 << 20 {
        group.sample_size(10);
        group.measurement_time(Duration::from_secs(20));
        group.warm_up_time(Duration::from_secs(5));
    } else if bytes > 5 << 20 {
        group.sample_size(20);
        group.measurement_time(Duration::from_secs(10));
    }
}

// ---------------------------------------------------------------------------
// Scenario 1: single file, one long sequence
// ---------------------------------------------------------------------------

fn bench_single_seq(c: &mut Criterion) {
    let path =
        std::env::var("BENCH_SINGLE_SEQ").unwrap_or_else(|_| CHRM_PATH.to_string());
    let label = file_stem(&path);
    let paths = vec![path.clone()];
    let file_bytes = total_file_bytes(&paths);
    let k = 5;
    let chunk_size = 1 << 20;

    eprintln!(
        "[single_seq] {} ({})",
        label,
        human_bytes(file_bytes)
    );

    // Verify: single-threaded == multi-threaded
    let count_t1 = run_smakcr(&paths, k, 1, chunk_size);
    let count_t4 = run_smakcr(&paths, k, 4, chunk_size);
    eprintln!(
        "  Accuracy check: t1={} t4={} {}",
        count_t1,
        count_t4,
        if count_t1 == count_t4 { "OK" } else { "MISMATCH!" }
    );
    assert_eq!(count_t1, count_t4, "single_seq count mismatch across threads");

    let mut group = c.benchmark_group(format!("e2e/single_seq/{}", label));
    group.throughput(Throughput::Bytes(file_bytes));
    tune_for_size(&mut group, file_bytes);

    for &t in &thread_counts() {
        group.bench_with_input(BenchmarkId::from_parameter(format!("t{}", t)), &t, |b, &t| {
            b.iter(|| run_smakcr(&paths, k, t, chunk_size))
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Scenario 2: single file, many sequences
// ---------------------------------------------------------------------------

fn bench_multi_seq(c: &mut Criterion) {
    let paths = vec![CHRUN_COMBINED.to_string()];
    let file_bytes = total_file_bytes(&paths);
    let k = 5;
    let chunk_size = 1 << 20;

    eprintln!(
        "[multi_seq] chrUn_combined ({})",
        human_bytes(file_bytes)
    );

    let count_t1 = run_smakcr(&paths, k, 1, chunk_size);
    let count_t4 = run_smakcr(&paths, k, 4, chunk_size);
    eprintln!(
        "  Accuracy check: t1={} t4={} {}",
        count_t1,
        count_t4,
        if count_t1 == count_t4 { "OK" } else { "MISMATCH!" }
    );
    assert_eq!(count_t1, count_t4, "multi_seq count mismatch across threads");

    let mut group = c.benchmark_group("e2e/multi_seq/chrUn_combined");
    group.throughput(Throughput::Bytes(file_bytes));
    tune_for_size(&mut group, file_bytes);

    for &t in &thread_counts() {
        group.bench_with_input(BenchmarkId::from_parameter(format!("t{}", t)), &t, |b, &t| {
            b.iter(|| run_smakcr(&paths, k, t, chunk_size))
        });
    }

    group.finish();
}

// ---------------------------------------------------------------------------
// Scenario 3: many files, one sequence each
// ---------------------------------------------------------------------------

fn bench_multi_file(c: &mut Criterion) {
    let paths = collect_files(CHRUN_PER_FILE_DIR);
    let file_bytes = total_file_bytes(&paths);
    let k = 5;
    let chunk_size = 1 << 20;

    eprintln!(
        "[multi_file] chrUn per-scaffold ({} files, {})",
        paths.len(),
        human_bytes(file_bytes)
    );

    let count_t1 = run_smakcr(&paths, k, 1, chunk_size);
    let count_t4 = run_smakcr(&paths, k, 4, chunk_size);
    eprintln!(
        "  Accuracy check: t1={} t4={} {}",
        count_t1,
        count_t4,
        if count_t1 == count_t4 { "OK" } else { "MISMATCH!" }
    );
    assert_eq!(count_t1, count_t4, "multi_file count mismatch across threads");

    let mut group = c.benchmark_group("e2e/multi_file/chrUn_per_scaffold");
    group.throughput(Throughput::Bytes(file_bytes));
    tune_for_size(&mut group, file_bytes);

    for &t in &thread_counts() {
        group.bench_with_input(BenchmarkId::from_parameter(format!("t{}", t)), &t, |b, &t| {
            b.iter(|| run_smakcr(&paths, k, t, chunk_size))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_single_seq, bench_multi_seq, bench_multi_file);
criterion_main!(benches);
