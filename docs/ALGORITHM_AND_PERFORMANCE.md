# Algorithm & Performance Architecture

[Back to README](../README.md)

This document explains how **PhyloTrim** achieves GB/s-class end-to-end alignment trimming throughput while maintaining a tiny, auditable codebase.

---

## 1. Execution Pipeline Overview

```mermaid
flowchart LR
    A["Input Path (File or Dir)"] --> B{"Directory?"}
    B -- Yes --> C["Collect Supported Files & LPT Sort (Largest First)"]
    C --> D["Rayon Work-Stealing Pool (File-Level Parallelism)"]
    B -- No --> E["Single-File Dispatch"]
    D --> F["Flat Contiguous Row-Major AlignmentMatrix"]
    E --> F
    F --> G["Auto SeqKind Detection (64 KB Sample)"]
    G --> H["L1-Cache Tiled (8,192 cols) + 255-Row u8 SIMD Gap Counting"]
    H --> I{"Codon Mode?"}
    I -- Yes --> J["3-bp Triplet Keep Mask"]
    I -- No --> K["Per-Column Keep Mask"]
    J --> L["Run-Length Coalescing (Contiguous Kept Spans)"]
    K --> L
    L --> M{"kept == orig?"}
    M -- Yes --> N["Zero-Copy Fast Path (Cow::Borrowed)"]
    M -- No --> O["Parallel Contiguous memcpy (extract_runs)"]
    N --> P["Atomic Tempfile + Rename Write"]
    O --> P
```

---

## 2. Memory Layout: Contiguous Row-Major Matrix

Instead of storing `Vec<Vec<u8>>` or per-character heap objects, PhyloTrim parses alignments directly into a single flat byte buffer:

```rust
pub struct AlignmentMatrix {
    pub names: Vec<String>,
    pub data: Vec<u8>,   // row r at [r * n_sites .. (r + 1) * n_sites]
    pub n_taxa: usize,
    pub n_sites: usize,
    pub seq_kind: SeqKind,
}
```

- **Single contiguous allocation** eliminates pointer chasing and allocator contention across threads.
- **Flat FASTA parser**: uses `memchr` to scan newlines and appends sequence bytes into one contiguous buffer before constructing the padded matrix.

---

## 3. L1-Cache Tiling & 255-Row `u8` SIMD Sub-Batching

Naive column-wise counting (`for col { for row { ... } }`) incurs cache misses on every step when `n_sites` is large. Conversely, naive row-wise accumulation into `Vec<u32>` spills out of L1 cache once `n_sites > 8,192`.

PhyloTrim resolves both bottlenecks simultaneously in `src/trim.rs`:

1. **Column Tiling (`TILE_COLS = 8,192`)**:
   - The `u8` local sub-accumulator occupies exactly **8 KB** (comfortably inside every modern CPU's 32–128 KB L1 data cache).
   - For large matrices (`n_taxa * usable_sites >= 512 KB` and `usable_sites > 8,192`), tiles are distributed across CPU cores via Rayon `par_chunks_mut`.
2. **255-Row Overflow-Free `u8` Micro-Batches**:
   - Within each tile, rows are processed in batches of at most `255` rows.
   - Because each row contributes `0` or `1`, the `u8` accumulator can never overflow within a 255-row batch, allowing 4× denser SIMD vector width (`16 × u8` in NEON / `32 × u8` in AVX2) before flushing into `u32` output counts.
3. **Branchless Nucleotide Classification**:
   - For DNA/RNA, `is_nucleotide_gap_byte(b)` performs a register-only bitmask (`b & !0x20`) and 5-way comparison with zero memory table lookups, enabling LLVM to auto-vectorize the inner loop cleanly.
   - For Protein, a 64-byte cache-line-aligned 256-entry `GapLut` (`#[repr(align(64))]`) resolves each byte in one L1 hit.

---

## 4. Zero-Copy Fast Path & Run-Length Coalescing

After computing the boolean `keep` mask:

- **Zero-Copy Fast Path**: If `kept_sites == orig_sites`, `trim_alignment` returns `Cow::Borrowed(matrix)` with zero allocation and zero byte copying.
- **Run-Length Coalescing**: Adjacent kept columns are merged into `(start_col, length)` spans. Instead of copying 1 byte at a time, each taxon row copies contiguous slices via `copy_from_slice` (`memcpy`), parallelized across rows when the output exceeds 512 KB.

---

## 5. LPT Scheduling & Atomic I/O

- **LPT (Largest-Processing-Time-First) Scheduling**: In directory mode, files are sorted by byte size descending before entering the Rayon work-stealing pool. Huge alignments start immediately on turn 0, preventing tail-latency stragglers.
- **Atomic File Writes**: Output is written to a hidden `.phylotrim.tmp.<pid>.<tid>.<seq>` file in the destination directory and committed via `std::fs::rename`, guaranteeing crash-safe outputs.
