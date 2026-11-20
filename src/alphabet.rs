/// Automatically detected biological sequence kind.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqKind {
    Dna,
    Rna,
    Protein,
}

impl SeqKind {
    #[inline]
    pub fn is_nucleotide(self) -> bool {
        matches!(self, SeqKind::Dna | SeqKind::Rna)
    }

    pub fn nexus_datatype(self) -> &'static str {
        match self {
            SeqKind::Dna => "DNA",
            SeqKind::Rna => "RNA",
            SeqKind::Protein => "PROTEIN",
        }
    }
}

/// Cache-line-aligned 256-entry lookup table mapping each ASCII byte to `0` (valid standard symbol)
/// or `1` (gap / missing / ambiguous / any non-standard symbol).
#[repr(align(64))]
#[derive(Debug, Clone)]
pub struct GapLut {
    table: [u8; 256],
    pub kind: SeqKind,
}

impl GapLut {
    pub fn new(kind: SeqKind) -> Self {
        let mut table = [1u8; 256];
        match kind {
            SeqKind::Dna | SeqKind::Rna => {
                // Only standard nucleotides A, T, C, G, U (case-insensitive) are valid (0).
                // Everything else (-, ., ?, N, IUPAC degenerate codes, etc.) counts as gap (1).
                for &b in b"ATCGUatcgu" {
                    table[b as usize] = 0;
                }
            }
            SeqKind::Protein => {
                // Only the 20 canonical amino acids (case-insensitive) are valid (0).
                // Everything else (-, ., ?, X, *, B, Z, J, U, O, etc.) counts as gap (1).
                for &b in b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy" {
                    table[b as usize] = 0;
                }
            }
        }
        Self { table, kind }
    }

    #[inline(always)]
    pub fn is_gap_byte(&self, b: u8) -> u8 {
        self.table[b as usize]
    }
}

/// Branchless register-only check for standard nucleotide bytes (`A, T, C, G, U`, case-insensitive).
/// Returns `0` for valid standard nucleotides and `1` for any gap/non-standard byte.
/// Enables LLVM auto-vectorization (NEON / AVX2) over contiguous row slices without memory lookups.
#[inline(always)]
pub fn is_nucleotide_gap_byte(b: u8) -> u8 {
    let u = b & !0x20;
    u8::from(!(u == b'A' || u == b'T' || u == b'C' || u == b'G' || u == b'U'))
}

/// Automatically detects whether alignment bytes represent DNA, RNA, or Protein.
///
/// Uses the 6 common amino acids absent from all IUPAC nucleotide codes (`E, F, I, L, P, Q`,
/// which typically comprise ~35% of real proteins) so nucleotide alignments rich in `N` or
/// degenerate IUPAC bases are never misclassified as Protein.
pub fn detect_seq_kind(data: &[u8]) -> SeqKind {
    const SAMPLE_LIMIT: usize = 65_536;
    let mut non_gap = 0usize;
    let mut acgtun = 0usize;
    let mut prot_exclusive = 0usize;
    let mut t_count = 0usize;
    let mut u_count = 0usize;

    for &b in data {
        match b {
            b'-' | b'?' | b'.' | b'*' | b' ' | b'\t' | b'\r' | b'\n' => continue,
            _ => {}
        }
        non_gap += 1;
        let upper = b & !0x20;
        match upper {
            b'A' | b'C' | b'G' | b'N' => acgtun += 1,
            b'T' => {
                acgtun += 1;
                t_count += 1;
            }
            b'U' => {
                acgtun += 1;
                u_count += 1;
            }
            b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' => prot_exclusive += 1,
            _ => {}
        }
        if non_gap >= SAMPLE_LIMIT {
            break;
        }
    }

    if non_gap == 0 {
        return SeqKind::Dna;
    }

    let prot_ratio = prot_exclusive as f64 / non_gap as f64;
    let acgtun_ratio = acgtun as f64 / non_gap as f64;

    if prot_ratio >= 0.02 || (prot_exclusive > 0 && acgtun_ratio < 0.75) {
        SeqKind::Protein
    } else if u_count > t_count {
        SeqKind::Rna
    } else {
        SeqKind::Dna
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_seq_kind_detection() {
        assert_eq!(detect_seq_kind(b"ATCGATCG--NN??"), SeqKind::Dna);
        assert_eq!(detect_seq_kind(b"AUCGAUCG--NN??"), SeqKind::Rna);
        // DNA rich in IUPAC degenerate symbols should still be detected as Dna, not Protein
        assert_eq!(detect_seq_kind(b"ATCGRYSWKMBDHVNNNN--"), SeqKind::Dna);
        assert_eq!(detect_seq_kind(b"MKTIIALSYIFCLVFAQDLP"), SeqKind::Protein);
        assert_eq!(detect_seq_kind(b"MNNELPQIFNNSKLV"), SeqKind::Protein);
    }

    #[test]
    fn test_gap_rules_nucleotide_and_protein() {
        let nuc_lut = GapLut::new(SeqKind::Dna);
        for &b in b"ATCGUatcgu" {
            assert_eq!(nuc_lut.is_gap_byte(b), 0);
            assert_eq!(is_nucleotide_gap_byte(b), 0);
        }
        for &b in b"-?.NnRYSWKMBDHVryswkmbdhvX*" {
            assert_eq!(nuc_lut.is_gap_byte(b), 1);
            assert_eq!(is_nucleotide_gap_byte(b), 1);
        }

        let prot_lut = GapLut::new(SeqKind::Protein);
        for &b in b"ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy" {
            assert_eq!(prot_lut.is_gap_byte(b), 0);
        }
        for &b in b"-?.Xx*BZJUObzjuo" {
            assert_eq!(prot_lut.is_gap_byte(b), 1);
        }
    }
}
