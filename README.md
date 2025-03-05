# SmaKCR
**Sma**ll **K**-mer **C**ounting with **R**ust.

Most k-mer counters are optimized for large k.
However, some (admittedly rather niche) use cases require only short k-mers.
SmaKCR is a simple k-mer counter that is faster than the most popular alternatives (Jellyfish and kmc) for k < 13.

The difference to Jellyfish is large, the difference to kmc is more moderate, but will improve with further development.
