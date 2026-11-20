use crate::alphabet::{detect_seq_kind, SeqKind};
use anyhow::{bail, Result};
use rayon::prelude::*;

/// Contiguous row-major byte representation of a multiple sequence alignment.
#[derive(Debug, Clone)]
pub struct AlignmentMatrix {
    pub names: Vec<String>,
    pub data: Vec<u8>,
    pub n_taxa: usize,
    pub n_sites: usize,
    pub seq_kind: SeqKind,
}

impl AlignmentMatrix {
    /// Builds an `AlignmentMatrix` directly from an already contiguous row-major `data` buffer.
    pub fn from_flat(names: Vec<String>, data: Vec<u8>, n_sites: usize) -> Result<Self> {
        if names.is_empty() {
            bail!("Alignment contains no sequences");
        }
        let n_taxa = names.len();
        if data.len() != n_taxa * n_sites {
            bail!(
                "Contiguous buffer size mismatch: expected {}, got {}",
                n_taxa * n_sites,
                data.len()
            );
        }
        let seq_kind = detect_seq_kind(&data);
        Ok(Self {
            names,
            data,
            n_taxa,
            n_sites,
            seq_kind,
        })
    }

    /// Constructs a validated matrix from per-taxon byte vectors, padding shorter rows with `'-'`.
    pub fn from_records(names: Vec<String>, seqs: Vec<Vec<u8>>) -> Result<Self> {
        if names.is_empty() {
            bail!("Alignment contains no sequences");
        }
        if names.len() != seqs.len() {
            bail!("Mismatched taxon names and sequence records");
        }

        let n_taxa = names.len();
        let max_len = seqs.iter().map(|s| s.len()).max().unwrap_or(0);
        let mut data = Vec::with_capacity(n_taxa * max_len);

        for seq in &seqs {
            data.extend_from_slice(seq);
            if seq.len() < max_len {
                data.resize(data.len() + (max_len - seq.len()), b'-');
            }
        }

        let seq_kind = detect_seq_kind(&data);
        Ok(Self {
            names,
            data,
            n_taxa,
            n_sites: max_len,
            seq_kind,
        })
    }

    #[inline(always)]
    pub fn row(&self, r: usize) -> &[u8] {
        let start = r * self.n_sites;
        &self.data[start..start + self.n_sites]
    }

    /// Extracts a column-trimmed matrix using run-length coalesced `copy_from_slice` (`memcpy`).
    /// Automatically parallelizes across rows when the output matrix is large (`>= 512 KB`).
    pub fn extract_runs(&self, runs: &[(usize, usize)], new_sites: usize) -> Self {
        let n_taxa = self.n_taxa;
        let orig_sites = self.n_sites;
        let total_out = n_taxa * new_sites;
        let mut out_data = vec![0u8; total_out];

        if new_sites > 0 {
            let copy_row = |(r, dst_row): (usize, &mut [u8])| {
                let src_row = &self.data[r * orig_sites..(r + 1) * orig_sites];
                let mut offset = 0usize;
                for &(start, len) in runs {
                    dst_row[offset..offset + len].copy_from_slice(&src_row[start..start + len]);
                    offset += len;
                }
            };

            if total_out >= 512 * 1024 && n_taxa >= 4 {
                out_data
                    .par_chunks_mut(new_sites)
                    .enumerate()
                    .for_each(copy_row);
            } else {
                out_data
                    .chunks_mut(new_sites)
                    .enumerate()
                    .for_each(copy_row);
            }
        }

        Self {
            names: self.names.clone(),
            data: out_data,
            n_taxa,
            n_sites: new_sites,
            seq_kind: self.seq_kind,
        }
    }
}
