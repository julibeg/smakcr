# Smakcr Chunking And Concurrency Architecture

## Overview

On the `parallel-reader` branch, smakcr uses a custom streaming FASTA/FASTQ
chunker for both single-threaded and multi-threaded counting. The chunker emits
sequence-only byte chunks and uses a sentinel byte to mark sequence boundaries so
k-mers never cross records.

For multi-threaded counting, worker threads share a pool of open chunkers and a
bounded chunk queue. Workers alternate between counting queued chunks and reading
more chunks from input files.

## Chunker

The chunker lives in `src/chunker.rs` and provides `FastxChunker`.

Important details:

- Format is detected from the first byte: `>` for FASTA, `@` for FASTQ.
- Parsing works directly on `BufRead::fill_buf()` buffers.
- `memchr`/`memchr2` are used to find headers, line endings, and record
  boundaries efficiently.
- `SEQ_SENTINEL` is `0xFF` and is inserted at sequence boundaries.
- `next_chunk()` prepends the previous chunk's last `k - 1` bytes to preserve
  k-mer continuity across chunk splits.
- The tail is currently just the last `k - 1` bytes of the previous chunk. It may
  include sentinels; counting remains correct because sentinel-aware counting
  resets state when it sees `SEQ_SENTINEL`.
- `chunk_size` is a target/limit checked before reads, not a strict final chunk
  length. A copied FASTA span or FASTQ sequence line can make a chunk exceed the
  nominal size.

### FASTA behavior

The FASTA parser:

- skips header lines,
- copies sequence spans,
- strips `\n` and `\r`,
- inserts a sentinel before starting a new record if the previous record emitted
  bases.

Other whitespace characters such as spaces and tabs are not stripped specially;
they flow into the counting path and are treated as invalid bases by `KEY_MAP`.

### FASTQ behavior

The FASTQ parser is currently a simple four-line state machine:

1. header line,
2. sequence line,
3. plus line,
4. quality line.

It copies the sequence line, strips a trailing `\r` if present, appends a
sentinel, then skips the plus and quality lines.

Current limitation: multi-line FASTQ records are not supported. If that behavior
is desired, the parser needs to accumulate sequence lines until `+`, then skip
quality bytes until the matching sequence length has been consumed.

## Worker Pool

The worker pool is implemented in `src/lib.rs` around `count_paths_worker_pool`.
It consists of:

- `FilePool`: a mutex/condvar-protected pool of `FileState` objects. Each
  `FileState` owns an open `FastxChunker`, so file/chunker state is preserved
  when a file is released and reacquired by another worker.
- `ChunkQueue`: a mutex/condvar-protected FIFO of `Vec<u8>` chunks with a byte
  capacity limit.
- `n_threads` worker threads, each counting into its own per-thread count vector.

## Reading And Backpressure

Each worker loop roughly does this:

1. Try to pop and count a queued chunk.
2. If no chunk is queued, check whether all file work is done and the queue is
   empty.
3. If possible, acquire a file from `FilePool` and read chunks from it.
4. Try to push each chunk to `ChunkQueue`.
5. If the queue is full, count that chunk locally in the same worker, release the
   file back to `FilePool`, and return to the top of the loop.
6. If the file reaches EOF, release it as finished.
7. If no file is available, wait on the queue condvar for more chunks or
   shutdown.

There is no separate `pending` chunk stored in `FileState`. Queue-full chunks are
processed immediately by the worker that produced them.

No file offsets are stored. The reader state lives inside `FastxChunker` and is
kept by moving the `FileState` in and out of `FilePool`.

## Counting

Chunk counting uses `count_chunk`, which delegates to `count_kmers_bytes` with
`Some(SEQ_SENTINEL)` as a reset byte.

Counting behavior:

- A sentinel resets the rolling k-mer index and sets `skip = k - 1`.
- Invalid bases, including any non-ACGT byte, also reset/skip so invalid data
  cannot contribute to k-mers.
- Lowercase and uppercase A/C/G/T are accepted through `KEY_MAP`.
- Per-worker count vectors are merged in `KmerCounter::into_vec()`.

## Single-Threaded Path

The single-threaded path also uses `FastxChunker`. It reads each input file
sequentially and counts each chunk with the same sentinel-aware logic used by the
worker pool.

## Known Limitations And Follow-Up Work

- CLI validation should reject invalid values cleanly instead of relying on
  panics/assertions, especially `k == 0`, `threads == 0`, too-large `k`, and
  `chunk_size < k`.
- Multi-line FASTQ support is not implemented.
- Worker-pool shutdown uses several separate `is_done()` / `is_empty()` checks;
  this should be simplified or carefully reviewed for races.
- `ChunkQueue::try_push` silently returns `Ok(())` if the queue is closed, which
  could hide logic bugs.
- The hand-rolled queue/file-pool could be replaced with a simpler channel-based
  design if correctness or contention becomes problematic.
- Tests should include forced tiny chunk sizes through the CLI/integration path,
  plus FASTA/FASTQ boundary cases.
