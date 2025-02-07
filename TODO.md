## Short term

- why is writing for k>9 currently so slow?
  - we need better handling of canonical kmers
- implement with hashmap for k > 9
  - benchmark for ks up to 13 or so and see where the break-even point is
- benchmark against
  - jellyfish
  - krust
  - kmc3
- more testing
  - test larger files (e.g. 10kb) and compare with expected hash of counts file or similar
  - automate testing (with pre-commit or github actions)

## Longer term

- add thread for decompression
- add multi-threading
- RNA and protein support
- restructure so that library can be used in other projects
