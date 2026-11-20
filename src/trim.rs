use crate::alphabet::{is_nucleotide_gap_byte, GapLut};
use crate::matrix::AlignmentMatrix;
use rayon::prelude::*;
use std::borrow::Cow;

/// Width of each column tile (`8,192` columns = 8 KB in `u8` sub-accumulator, 32 KB in `u32` counts),
/// fitting comfortably inside a single CPU core's L1 data cache.
const TILE_COLS: usize = 8_192;

/// Outcome of trimming a single alignment.
pub struct TrimOutcome<'a> {
    pub matrix: Cow<'a, AlignmentMatrix>,
    pub orig_sites: usize,
    pub kept_sites: usize,
}

/// Computes the exact maximum allowed gap count per column for `n_taxa` rows and `max_gap` in `[0.0, 1.0]`.
/// A column is kept iff `gap_count <= max_allowed_gaps` (i.e. `gap_count / n_taxa <= max_gap`).
#[inline]
pub fn max_allowed_gaps(n_taxa: usize, max_gap: f64) -> u32 {
    if max_gap <= 0.0 {
        0
    } else if max_gap >= 1.0 {
        n_taxa as u32
    } else {
        ((max_gap + 1e-12) * (n_taxa as f64))
            .floor()
            .min(n_taxa as f64) as u32
    }
}

/// Counts non-standard/gap characters per column across `[0..usable_sites]` using L1-cache-resident
/// row-major tiles and 255-row `u8` SIMD sub-batches.
///
/// Automatically parallelizes across tiles when the matrix is large (`>= 512 KB`).
pub fn count_column_gaps(matrix: &AlignmentMatrix, usable_sites: usize) -> Vec<u32> {
    let mut counts = vec![0u32; usable_sites];
    if usable_sites == 0 || matrix.n_taxa == 0 {
        return counts;
    }

    let n_taxa = matrix.n_taxa;
    let stride = matrix.n_sites;
    let data = &matrix.data;
    let is_nuc = matrix.seq_kind.is_nucleotide();
    let prot_lut = if is_nuc {
        None
    } else {
        Some(GapLut::new(matrix.seq_kind))
    };

    let process_tile = |(tile_idx, out_tile): (usize, &mut [u32])| {
        let col_start = tile_idx * TILE_COLS;
        let tile_len = out_tile.len();
        let mut local_u8 = [0u8; TILE_COLS];

        let mut row_base = 0usize;
        while row_base < n_taxa {
            let batch_end = (row_base + 255).min(n_taxa);
            local_u8[..tile_len].fill(0);

            if is_nuc {
                for r in row_base..batch_end {
                    let off = r * stride + col_start;
                    let row_slice = &data[off..off + tile_len];
                    for (&b, acc) in row_slice.iter().zip(local_u8[..tile_len].iter_mut()) {
                        *acc += is_nucleotide_gap_byte(b);
                    }
                }
            } else if let Some(ref lut) = prot_lut {
                for r in row_base..batch_end {
                    let off = r * stride + col_start;
                    let row_slice = &data[off..off + tile_len];
                    for (&b, acc) in row_slice.iter().zip(local_u8[..tile_len].iter_mut()) {
                        *acc += lut.is_gap_byte(b);
                    }
                }
            }

            for (dst, &src) in out_tile.iter_mut().zip(local_u8[..tile_len].iter()) {
                *dst += src as u32;
            }
            row_base = batch_end;
        }
    };

    if n_taxa * usable_sites >= 512 * 1024 && usable_sites > TILE_COLS {
        counts
            .par_chunks_mut(TILE_COLS)
            .enumerate()
            .for_each(process_tile);
    } else {
        counts
            .chunks_mut(TILE_COLS)
            .enumerate()
            .for_each(process_tile);
    }

    counts
}

/// Coalesces a boolean keep mask into contiguous `(start_col, length)` runs and total kept sites.
fn coalesce_runs(keep: &[bool]) -> (Vec<(usize, usize)>, usize) {
    let mut runs = Vec::new();
    let mut total_kept = 0usize;
    let mut idx = 0usize;
    let n = keep.len();

    while idx < n {
        if keep[idx] {
            let start = idx;
            idx += 1;
            while idx < n && keep[idx] {
                idx += 1;
            }
            let len = idx - start;
            runs.push((start, len));
            total_kept += len;
        } else {
            idx += 1;
        }
    }

    (runs, total_kept)
}

/// Executes fast column-gap or 3-bp codon trimming on `matrix`.
///
/// - Never drops or alters sequence rows (`n_taxa` and taxon order are 100% preserved).
/// - Zero-copy fast path (`Cow::Borrowed`) when all columns pass the threshold.
pub fn trim_alignment<'a>(
    matrix: &'a AlignmentMatrix,
    max_gap: f64,
    codon_mode: bool,
) -> TrimOutcome<'a> {
    let orig_sites = matrix.n_sites;
    let usable_sites = if codon_mode {
        (orig_sites / 3) * 3
    } else {
        orig_sites
    };

    if usable_sites == 0 {
        if orig_sites == 0 {
            return TrimOutcome {
                matrix: Cow::Borrowed(matrix),
                orig_sites: 0,
                kept_sites: 0,
            };
        }
        return TrimOutcome {
            matrix: Cow::Owned(matrix.extract_runs(&[], 0)),
            orig_sites,
            kept_sites: 0,
        };
    }

    let max_gaps = max_allowed_gaps(matrix.n_taxa, max_gap);
    let gap_counts = count_column_gaps(matrix, usable_sites);

    let mut keep = vec![false; usable_sites];
    if codon_mode {
        for codon_idx in 0..(usable_sites / 3) {
            let base = codon_idx * 3;
            if gap_counts[base] <= max_gaps
                && gap_counts[base + 1] <= max_gaps
                && gap_counts[base + 2] <= max_gaps
            {
                keep[base] = true;
                keep[base + 1] = true;
                keep[base + 2] = true;
            }
        }
    } else {
        for (flag, &cnt) in keep.iter_mut().zip(gap_counts.iter()) {
            *flag = cnt <= max_gaps;
        }
    }

    let (runs, kept_sites) = coalesce_runs(&keep);

    // Zero-copy fast path when every column in the original alignment is kept.
    if kept_sites == orig_sites {
        return TrimOutcome {
            matrix: Cow::Borrowed(matrix),
            orig_sites,
            kept_sites,
        };
    }

    let trimmed = matrix.extract_runs(&runs, kept_sites);
    TrimOutcome {
        matrix: Cow::Owned(trimmed),
        orig_sites,
        kept_sites,
    }
}
