# Changelog

All notable changes to **PhyloTrim** will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-11-16

### Added
- Minimalist 5-flag CLI (`-i`, `-o`, `-g`, `-c/--codon`, `-t`) focused exclusively on column-gap and 3-bp codon trimming.
- Automatic DNA / RNA / Protein alphabet detection with strict canonical symbol validation.
- L1-cache tiled (`8,192` columns) + 255-row overflow-free `u8` SIMD gap counting engine.
- Zero-copy fast path (`Cow::Borrowed`) when all columns pass and run-length coalesced parallel `memcpy` extraction.
- LPT (Largest-Processing-Time-First) directory scheduling across Rayon worker threads.
- Native FASTA, PHYLIP (sequential/interleaved), NEXUS, and `.gz` reading/writing with atomic file commits.
- Comprehensive unit and integration test suite plus GitHub Actions CI (`rustfmt`, `clippy`, `cargo test --release`).
