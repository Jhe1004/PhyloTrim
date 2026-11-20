# PhyloTrim User Guide

[简体中文版](zh-CN/USER_GUIDE.md) | [Back to README](../README.md)

---

## 1. Supported Formats & Compression

PhyloTrim automatically sniffs alignment format from file content (with UTF-8 BOM stripping) and natively supports gzip compression (`.gz`):

| Format | Recognized Extensions (Directory Mode) | Notes |
| :--- | :--- | :--- |
| **FASTA** | `.fa`, `.fas`, `.fasta`, `.fna`, `.ffn`, `.faa`, `.frn`, `.aln` | Multi-line or single-line sequences; unequal row lengths are right-padded with `-` during parsing. |
| **PHYLIP** | `.phy`, `.phylip` | Supports both Sequential and Interleaved variants. |
| **NEXUS** | `.nex`, `.nexus` | Supports standard `DATA` / `CHARACTERS` blocks, quoted taxon names, and interleaved matrices. |
| **Gzip** | Any of the above + `.gz` | Detected via gzip magic bytes (`1F 8B`) on input; written compressed when input was `.gz` (directory mode) or output ends with `.gz` (single-file mode). |

### Format Conversion in Single-File Mode

In single-file mode, PhyloTrim infers the output format from the output file extension (falling back to the input format if unrecognized):

```bash
# Read gzipped FASTA, trim at 20% max gap, and write PHYLIP directly
phylotrim -i gene1.fas.gz -o gene1.trimmed.phy -g 0.2
```

In directory mode, each file retains its original filename, format, and `.gz` status in the output directory, and **no extra summary/log files are written into the output directory** so downstream batch tools can glob cleanly.

---

## 2. Exact Threshold & Boundary Semantics

Given `n_taxa` sequences and threshold `max_gap` $\in [0.0, 1.0]$:

- A column with `gap_count` non-standard characters is **kept** if and only if:
  $$\frac{\text{gap\_count}}{n_{\text{taxa}}} \le \text{max\_gap}$$
- **Zero-tolerance (`-g 0.0`)**: Only columns with `gap_count == 0` are retained.
- **Full-retention (`-g 1.0`)**: Every column is retained (triggers the zero-copy fast path).
- **All-gap rows are never removed**: Even if a taxon consists entirely of `-` after trimming (or all columns are dropped to length `0`), every taxon header is preserved in its exact original order.

---

## 3. Codon Mode (`-c / --codon`)

When `--codon` is enabled:

1. Input must be a nucleotide alignment (DNA or RNA); protein alignments are rejected with a clear error.
2. Trailing bases beyond the last complete triplet (`n_sites % 3`) are excluded so the output length is guaranteed to be a multiple of 3.
3. For each triplet `(3k, 3k+1, 3k+2)`, the triplet is **kept** if and only if **all three** columns individually satisfy `gap_count <= max_allowed_gaps`. Otherwise all 3 columns are dropped together.

---

## 4. Pipeline Integration Examples

### Multi-Locus Gene Tree Workflow (MAFFT → PhyloTrim → IQ-TREE → ASTRAL)

```bash
# 1. Trim all MAFFT alignments in parallel
phylotrim -i 01_mafft_loci/ -o 02_trimmed_loci/ -g 0.2

# 2. Infer gene trees with IQ-TREE
mkdir -p 03_genetrees
for f in 02_trimmed_loci/*.fas; do
  iqtree2 -s "$f" -m MFP -B 1000 -T AUTO --prefix "03_genetrees/$(basename "$f" .fas)"
done
```

### Coding Sequence (CDS) Workflow

```bash
# Preserve reading frame for downstream dN/dS (PAML/HyPhy) or codon-partitioned ML inference
phylotrim -i aligned_cds/ -o trimmed_cds/ -g 0.3 --codon
```

---

## 5. Frequently Asked Questions

**Q: Why does `phylotrim` treat `N` and IUPAC ambiguity codes (`R`, `Y`, `S`, `W`, `K`, `M`, ...) as gaps in DNA alignments?**  
A: In phylogenetic gap-proportion filtering, degenerate or undetermined bases carry missing/ambiguous character states rather than resolved nucleotides. Counting only canonical `A/T/C/G/U` as valid ensures strict, reproducible column occupancy.

**Q: What happens if `-i` and `-o` point to the same path?**  
A: `phylotrim` validates canonical paths upfront and exits with an error to prevent accidental in-place overwrite.

**Q: Is output writing crash-safe?**  
A: Yes. Every output file is written to a hidden temporary file in the target directory and atomically moved into place via `rename(2)`, preventing truncated or half-written files if a job is interrupted.
