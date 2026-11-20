<div align="center">

# PhyloTrim

**Ultra-fast, minimalist column-gap and 3-bp codon MSA trimmer written in Rust**

[![CI](https://github.com/Jhe1004/PhyloTrim/actions/workflows/ci.yml/badge.svg)](https://github.com/Jhe1004/PhyloTrim/actions/workflows/ci.yml)
[![Version](https://img.shields.io/badge/version-0.2.0-blue.svg)](CHANGELOG.md)
[![Rust](https://img.shields.io/badge/rust-2021%20edition-orange.svg)](https://www.rust-lang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#installation)

**English** | [简体中文](README.zh-CN.md)

[Features](#why-phylotrim) • [Installation](#installation) • [Quick Start](#quick-start) • [CLI Reference](#cli-reference) • [Documentation](#documentation) • [Performance](#performance-snapshot)

</div>

---

## Why PhyloTrim?

Phylogenomic pipelines routinely trim thousands of single-copy ortholog alignments or mega-base supermatrices before tree inference. **PhyloTrim** (`phylotrim`) strips away heuristic bloat and focuses on doing **two jobs with mathematical determinism and hardware-level speed**:

| Capability | Design Guarantee |
| :--- | :--- |
| **Column-Gap Trimming (Default)** | Removes any column whose fraction of non-standard/gap characters strictly exceeds `-g` (`gap_count / n_taxa > max_gap`). |
| **Frame-Safe Codon Trimming (`-c`)** | Evaluates aligned coding sequences in non-overlapping 3-bp triplets; if **any** of the 3 positions exceeds `-g`, the entire codon is dropped to preserve reading frame. |
| **Never Drops Taxa** | Strictly column-only trimming — every input taxon row and its original ordering are preserved 100%. |
| **Auto Nucleotide / Protein Detection** | Automatically distinguishes DNA/RNA vs. Protein alignments without manual flags, applying exact canonical alphabet rules. |
| **Zero-Dependency Native Speed** | Single statically linked Rust binary with LPT directory scheduling, L1-cache tiled SIMD counting, and run-length coalesced zero-copy output. |

---

## Installation

### From Source (Recommended)

Requires a recent [Rust toolchain](https://rustup.rs/) (`cargo >= 1.70`):

```bash
git clone https://github.com/Jhe1004/PhyloTrim.git
cd PhyloTrim
cargo build --release
# Binary is built at target/release/phylotrim
```

Or install directly into `$HOME/.cargo/bin`:

```bash
cargo install --git https://github.com/Jhe1004/PhyloTrim.git
```

Verify installation:

```bash
phylotrim --version
```

---

## Quick Start

```bash
# 1. Trim a single alignment: drop columns with >20% non-standard/gap symbols
phylotrim -i alignment.fasta -o alignment.trimmed.fasta -g 0.2

# 2. Codon mode: drop an entire 3-bp codon if any of its 3 columns exceeds 30% gaps
phylotrim -i cds.fasta -o cds.trimmed.fasta -g 0.3 --codon

# 3. Batch-trim an entire directory of thousands of loci across all CPU cores (supports .gz)
phylotrim -i raw_loci/ -o trimmed_loci/ -g 0.2

# 4. Trim a chromosome-scale gzipped alignment with 16 worker threads
phylotrim -i chr1_supermatrix.fas.gz -o chr1_trimmed.fas.gz -g 0.1 -t 16
```

---

## CLI Reference

`phylotrim` intentionally exposes only **5 options**:

```text
Usage: phylotrim [OPTIONS] --input <PATH> --output <PATH>

Options:
  -i, --input <PATH>     Input alignment file or directory (FASTA / PHYLIP / NEXUS, optional .gz)
  -o, --output <PATH>    Output alignment file or directory (optional .gz)
  -g, --max-gap <FLOAT>  Maximum allowed gap / non-standard symbol fraction per column [0.0..=1.0] [default: 0.2]
  -c, --codon            Trim by 3-bp codons (drop entire codon if any of its 3 columns exceeds --max-gap)
  -t, --threads <INT>    Number of worker threads (0 = use all available CPU cores) [default: 0]
  -h, --help             Print help
  -V, --version          Print version
```

---

## Character Classification Rules

Sequence type is detected automatically from the alignment content:

| Sequence Type | Valid Canonical Characters (Case-Insensitive) | Counted as Gap / Non-Standard |
| :--- | :--- | :--- |
| **DNA / RNA** | `A, T, C, G, U` | `-`, `.`, `?`, `N`, and all IUPAC degenerate codes (`R, Y, S, W, K, M, B, D, H, V`, etc.) |
| **Protein** | 20 canonical amino acids (`A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y`) | `-`, `.`, `?`, `X`, `*`, `B`, `Z`, `J`, `U`, `O`, etc. |

> **Safety guard**: Passing `--codon` on a Protein alignment exits with an explicit error rather than silently corrupting amino acid columns.

---

## Performance Snapshot

End-to-end wall-clock throughput (including disk I/O, parsing, trimming, and atomic writing) measured on a local Apple Silicon workstation (`phylotrim v0.2.0`, default `-t 0`):

| Workload | Matrix Dimensions / File Count | Mode | Wall-Clock Time |
| :--- | :--- | :--- | :---: |
| Single huge matrix (38.1 MB FASTA) | 200 taxa × 200,000 sites | Column (`-g 0.2`) | **0.030 s** |
| Single huge matrix (38.1 MB FASTA) | 200 taxa × 200,000 sites | Codon (`-g 0.2 --codon`) | **0.025 s** |
| Multi-locus directory (~75 MB total) | 1,000 files × (50 taxa × 1,500 bp) | Column (`-g 0.2`) | **0.110 s** |

See [Algorithm & Performance Architecture](docs/ALGORITHM_AND_PERFORMANCE.md) for how L1-cache tiling, branchless SIMD accumulation, and LPT scheduling achieve GB/s-class throughput.

---

## Documentation

- **[User Guide (English)](docs/USER_GUIDE.md)** — Supported formats, directory vs. single-file behavior, codon edge cases, and phylogenomic pipeline integration (MAFFT → PhyloTrim → IQ-TREE / ASTRAL).
- **[Algorithm & Performance Architecture](docs/ALGORITHM_AND_PERFORMANCE.md)** — Internal memory layout, L1-cache tiling, 255-row `u8` SIMD sub-batching, run-length coalescing, and atomic I/O design.
- **[用户手册（简体中文）](docs/zh-CN/USER_GUIDE.md)** — 中文完整使用说明与常见问题解答。
- **[Changelog](CHANGELOG.md)** — Release history and version notes.

---

## License

PhyloTrim is released under the [MIT License](LICENSE).
