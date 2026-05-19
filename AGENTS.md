# AGENTS.md

This file provides guidance to coding agents when working with code in this repository.

## What is smakcr?

**Sma**ll **K**-mer **C**ounting with **R**ust. A bioinformatics tool for counting k-mers in FASTA/FASTQ files, optimized for small k-values (k < 13). Faster than Jellyfish and KMC for small k.

## Build & Test Commands

```bash
cargo build --release          # release build (LTO + codegen-units=1)
cargo test                     # run all tests (unit + integration)
cargo test --test integration  # run only integration tests
cargo test -- test_name        # run a single test by name
cargo clippy                   # lint
cargo fmt                      # format
```

## Running

```bash
cargo run --release -- input.fa -k 5              # count 5-mers, single-threaded
cargo run --release -- input.fa -k 5 -t 4          # 4 threads
cargo run --release -- input.fa -k 5 -c            # canonical k-mers
cargo run --release -- input.fa -k 5 -z            # include zero-count k-mers
cargo run --release -- input.fa -k 5 -o output.tsv # write to file
```

## Architecture

The codebase is split into `src/lib.rs` (core library) and `src/bin/smakcr.rs` (CLI).

### Core types in `lib.rs`

- **`KmerCounter`** — enum with `SingleThreaded` / `MultiThreaded` variants. Holds k, a bitmask for rolling index calculation, and count vectors. Call `.count(input)` then `.write()` or `.into_vec()`.
- **`FastxReader`** — enum wrapping `seq_io` FASTA/FASTQ readers. Auto-detects format from first byte (`>` vs `@`). Handles compressed files via `niffler`.
- **`CountKmers`** trait — implemented for `FastxReader` (file-level counting) and for individual `RefRecord` types (record-level). The parallel path uses `seq_io::parallel::read_parallel`.
- **`PerBaseRefRecord`** trait — abstracts base iteration over FASTA records (handles multi-line sequences without allocating via `seq_lines()`) and FASTQ records.

### K-mer encoding

Bases are encoded as 2 bits (A=0, C=1, G=2, T=3) via the compile-time `KEY_MAP` lookup table. K-mers are packed into `usize` indices, giving O(1) count array access. Invalid bases (N, etc.) set a skip counter to avoid counting affected k-mers.

### Multi-threading strategy

- **Single file, multiple threads**: uses `seq_io::parallel::read_parallel` with per-thread count vectors (mutex-protected), merged at the end.
- **Multiple files, multiple threads**: uses `rayon` to process files in parallel with per-file `SingleThreaded` counters, reduced via the `Add` impl on `KmerCounter`.

## Test data

All test data lives in `tests/fixtures/` to be self-contained:

- `test.fna` — small FASTA with edge cases (soft/hard masking, non-canonical bases, multi-line sequences)
- `chrM.fa.gz` — human mitochondrial genome (single sequence)
- `chrUn_combined.fa.gz` — 39 unplaced scaffolds in one file (exercises single-file multi-threaded path)
- `chrUn-one-seq-per-file/` — same 39 scaffolds as individual files (exercises multi-file rayon path)
- `*_kmc.tsv` — KMC reference outputs used as ground truth in integration tests

KMC binaries for regenerating references are at `benchmark/kmc/bin/`. To regenerate a reference:
```bash
# non-canonical (-b disables canonical mode in KMC)
benchmark/kmc/bin/kmc -b -t1 -sf1 -sp1 -sr1 -fm -k5 -ci1 -cs9999999999999 INPUT TMPPREFIX TMPDIR
benchmark/kmc/bin/kmc_tools transform TMPPREFIX dump OUTPUT.tsv

# canonical (omit -b)
benchmark/kmc/bin/kmc -t1 -sf1 -sp1 -sr1 -fm -k5 -ci1 -cs9999999999999 INPUT TMPPREFIX TMPDIR
benchmark/kmc/bin/kmc_tools transform TMPPREFIX dump OUTPUT.tsv
```
