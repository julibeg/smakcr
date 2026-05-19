# Parallel Reader Branch Handoff

## Status

Branch: `parallel-reader`

Goal: finish and validate the custom FASTX chunk reader plus worker-pool counting path, then decide whether/how to merge it into `main`.

Current implementation is committed on this branch. The remaining context files in this working tree are mostly untracked Markdown notes/reviews. This handoff consolidates the useful parts so a different machine/session does not need all scratch artifacts.

Important branch facts:

- `parallel-reader` contains the implementation commit `71c7bc7 implement custom Fastx-reader and concurrency model` plus `bf3196f add AGENTS.md`.
- It is local-only unless pushed explicitly: `git push -u origin parallel-reader`.
- It is based before latest `main` commit `5cc6f92 upgrade deps`; rebase/merge may be needed before PR.
- `main`/`origin/main` still use the older `seq_io`/Rayon path and do not have `src/chunker.rs`.

## Files carrying the implementation

Committed on `parallel-reader`:

- `src/chunker.rs`
  - Custom `FastxChunker`.
  - Detects FASTA vs FASTQ from first byte.
  - Uses `BufRead::fill_buf()` and `memchr`/`memchr2` for bulk parsing.
  - Emits `SEQ_SENTINEL = 0xFF` between sequence records.
  - Produces fixed-size-ish chunks with `k - 1` overlap via `tail`.
- `src/lib.rs`
  - Wires `FastxChunker` into `KmerCounter::count_paths`.
  - Adds worker-pool path for `n_threads > 1`.
  - Adds `FilePool`, `ChunkQueue`, `count_chunk`, and shared `count_kmers_bytes`.
  - Removes old `CountKmers`/`seq_io::parallel::read_parallel` counting path.
- `src/bin/smakcr.rs`
  - CLI always constructs `KmerCounter` and calls `count_paths`.
  - Adds `--chunk-size` and `--queue-mb` knobs.
- `Cargo.toml`
  - Adds runtime `memchr`.
  - Keeps `seq_io` only as a dev-dependency for benchmarks.
- `benches/parser.rs`
  - Parser-only benchmark comparing `FastxChunker`, `seq_io`, and `needletail`.
- `benches/e2e.rs`
  - End-to-end benchmark.
- `tests/integration.rs`
  - Mostly formatting plus branch test adjustments.
- `AGENTS.md`
  - Currently stale: still describes the old `seq_io`/Rayon architecture. Update before merging.

## Relationship between old plans/reviews

### `PLAN-CONCURRENCY.md`

Stale/superseded.

It describes the original desired architecture: worker pool, chunk queue, sequence-boundary sentinel, `k - 1` overlap, and unified handling of single large files and many small files. The `parallel-reader` branch is the implementation of that intent, but not exactly the same design: the branch moved parsing/chunking into `src/chunker.rs` instead of keeping a `FileState` reader state machine entirely inside `src/lib.rs`.

Keep only as historical context unless you want provenance.

### `REVIEW-worker-pool-chunking.md`

Partly still relevant.

Useful findings to keep in mind:

- FASTQ parser is a 4-line state machine and does not support multi-line FASTQ records. This is a behavior change vs `seq_io`. Either document/validate single-line-only FASTQ or implement proper multiline FASTQ parsing.
- `k == 0` has mostly been addressed by `KmerCounter::new` asserting `k > 0`, but CLI still accepts `-k 0` and then panics. Prefer user-facing validation/error instead of panic.
- `ChunkQueue::try_push` returns `Ok(())` if closed, silently dropping a chunk. This should probably be an error or at least a debug assertion.

Some findings are already fixed in current code:

- Duplicate counting logic was consolidated into `count_kmers_bytes`.
- `KmerCounter::new` now asserts `k > 0`.

### `REVIEW-fastx-chunker-concurrency.md`

Partly stale, partly still useful.

The review was written against a slower custom chunker that used byte-at-a-time parsing and intermediate buffers. The current `src/chunker.rs` already implements the biggest recommended fix: `fill_buf()` plus `memchr` bulk parsing. So do not blindly follow every finding.

Still useful findings:

- Worker-pool termination logic is hard to reason about. Current code still checks `file_pool.is_done()` and `chunk_queue.is_empty()` separately in multiple places.
- `FilePool::acquire` can set `done` when no files are queued and no readers are active. Re-check this state machine carefully around races and shutdown.
- Workers lock their per-thread count vector for their whole lifetime. This avoids repeated locking but is odd because each vector is already thread-specific. Consider changing the storage shape or taking ownership per thread.
- Hand-rolled `ChunkQueue` has mutex/condvar complexity and possible contention. A bounded crossbeam channel or simpler dedicated reader/worker split may be easier to reason about.
- Per-chunk `Vec` allocation may remain a perf issue. Consider reusable buffers only after correctness and simpler concurrency are solid.

Already implemented/currently obsolete findings:

- Replace byte-at-a-time parser with bulk parsing: done.
- Remove `FastxSymbolStream`, `FastaStream`/`FastqStream`, pending `VecDeque`: done.
- Remove old `CountKmers` trait/read_parallel path: done.
- Move `bio` to dev-only: effectively done.
- Remove `rayon`: done.

### `PLAN-OPTIMIZE-FASTX-CHUNKER.md`

Partly already implemented.

Current status by step:

- Step 1, Cargo cleanup/build issue: mostly done. Runtime deps are lean; `seq_io` is dev-only. CLI still should validate `k >= 1`, `threads >= 1`, `chunk_size >= k`, `queue_mb > 0` with clean errors.
- Step 2, memchr-accelerated chunker rewrite: mostly done in `src/chunker.rs`.
- Step 2e, tests: partly done. `src/chunker.rs` unit tests already compare final k-mer counts across small chunk sizes for FASTA/FASTQ. Remaining gap is integration/CLI coverage across tiny chunk sizes and real fixture paths.
- Step 3, worker pool queue cleanup: the specific plan items are mostly done (`try_push` uses `notify_one`, `try_pop` no longer notifies, and queue-full chunks are counted locally). Broader shutdown/backpressure simplification is still not done: current code still has custom `ChunkQueue`/`FilePool` with mutex/condvar.
- Step 4, remove vestigial `CountKmers`: done.
- Step 5, final cleanup/verification: not done.

### `docs/ARCHITECTURE.md`

Useful as a short architecture writeup, but currently untracked and partly stale. If keeping docs, update it together with this handoff and `AGENTS.md` so all architecture docs agree. Known stale points: it mentions `pending` chunks even though queue-full chunks are counted locally, claims FASTQ whitespace filtering though current FASTQ parsing only strips trailing `\r` from sequence lines, and says overlap is only within the same sequence while current `tail` is simply the last `k - 1` bytes and may include sentinels.

### `SIMD-KMER-FINDINGS.md`

Not directly required for finishing `parallel-reader`. Useful future performance context: full SIMD exact k-mer counting is hard; prioritize parser/chunking, output path, and targeted byte classification experiments before any full SIMD rewrite.

## Current architecture summary

### Chunker

`FastxChunker` is format-aware and returns sequence-only chunks:

- FASTA: skips header lines, copies sequence spans, strips `\n`/`\r`, inserts sentinel before each new record.
- FASTQ: currently assumes records are exactly four lines: header, sequence, plus, quality. It copies the sequence line and inserts sentinel after it.
- Chunk overlap: `next_chunk()` prepends the previous `k - 1` bytes to the next chunk. Sentinels are preserved in chunks and `count_kmers_bytes` resets on sentinel, so k-mers should not cross sequence boundaries.

Main correctness edge to review: tail handling around sentinels and chunk boundaries, especially tiny `--chunk-size`, empty records, records shorter than `k`, and chunks ending immediately before/after sentinels.

### Worker pool

`count_paths_worker_pool` currently:

1. Builds one `FastxChunker` per input file.
2. Stores them in `FilePool`.
3. Starts `n_threads` workers.
4. Each worker:
   - counts queued chunks first,
   - otherwise acquires a file and produces chunks,
   - pushes chunks to `ChunkQueue` until queue capacity is hit,
   - if queue is full, counts the current chunk locally and releases the file back to the pool.
5. Each worker counts into its own per-thread vector; final merge happens in `into_vec()`.

This model can cover both one big file and many small files, but the shutdown/backpressure state machine needs careful review.

## Highest-priority remaining work

1. Validate correctness against old implementation/KMC fixtures.
   - Force very small chunks (`--chunk-size 1`, `2`, `k`, `k + 1`, etc. where allowed). Note: the parser may copy whole FASTA spans/FASTQ sequence lines past the nominal limit, so design tests with controlled input/buffer boundaries rather than assuming every chunk is exactly byte-sized.
   - Test single-sequence FASTA, multi-record FASTA, multi-line FASTA, short records, records with N/lowercase, sequence boundary at chunk edge.
   - Test single-line FASTQ explicitly.
   - Decide/document/implement multiline FASTQ behavior.
2. Simplify or harden worker-pool shutdown.
   - Avoid separate non-atomic `is_done()` + `is_empty()` decisions where possible.
   - Do not silently drop chunks when queue is closed.
   - Consider replacing `ChunkQueue` with `crossbeam::channel::bounded` or a clearer single reader / worker split.
3. Clean CLI validation.
   - `-k 0` should be a clean clap/user error, not an assertion panic.
   - Add a sane upper bound for `k`; `4usize.pow(k as u32)` can overflow or attempt huge allocations.
   - Validate `--chunk-size >= k`, `--queue-mb > 0`, and `--threads >= 1`.
4. Update docs/comments.
   - `AGENTS.md` is stale and still says `seq_io::parallel::read_parallel`/Rayon.
   - `docs/ARCHITECTURE.md` is useful but has the stale points listed above.
   - `tests/integration.rs` still has comments naming the old `seq_io::parallel` and Rayon paths; update wording to the worker-pool path.
   - Either commit updated `docs/ARCHITECTURE.md` or fold the architecture into README/AGENTS.
5. Run verification once the environment can download/build crates.
   - `cargo fmt`
   - `cargo test`
   - `cargo test --test integration`
   - `cargo clippy`
   - Optional: `cargo bench --bench parser`, `cargo bench --bench e2e` on representative FASTA/FASTQ.

## Validation note from this session

I attempted `cargo test --quiet`, but the environment failed before compiling because Cargo could not write into the global registry cache:

```text
failed to create directory /home/user/.cargo/registry/src/.../anes-0.1.6
Read-only file system (os error 30)
```

Per project instructions, I did not try to fix the environment. Re-run tests on a writable Cargo setup.

## What to commit for continuation

Minimum useful continuation commit:

- `PARALLEL-READER-HANDOFF.md` (this file)
- Updated `AGENTS.md` if you fix the stale architecture section
- Optional: updated/committed `docs/ARCHITECTURE.md` if you want a shorter architecture doc in-repo

Old untracked Markdown can then be skipped if this handoff is enough:

- `PLAN-CONCURRENCY.md` — stale/superseded
- `PLAN-OPTIMIZE-FASTX-CHUNKER.md` — partially implemented, now summarized here
- `REVIEW-fastx-chunker-concurrency.md` — partially stale, summarized here
- `REVIEW-worker-pool-chunking.md` — summarized here
- `SIMD-KMER-FINDINGS.md` — optional future perf context
- `CODE_REVIEW.md` — older/general review, not specific to this branch
- `tmp.md` — original prompt/provenance, not needed

Do not commit as-is unless intentionally curating them:

- `benchmark/` — huge, includes raw genomes, KMC binaries, logs, and generated results
- `out.csv`, `test.out`, `tmp.out` — generated output
- `.vscode/` — local editor settings
- `src/lib-old.rs`, `src/bin/legacy-versions/` — old snapshots/reference code
