use flate2::write::GzEncoder;
use flate2::Compression;
use phylotrim::cli::CliArgs;
use phylotrim::engine::run;
use phylotrim::io::read_alignment;
use std::fs;
use std::io::Write;
use tempfile::tempdir;

#[test]
fn test_single_column_gap_threshold_and_non_atcg_symbols() {
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("dna.fas");
    let out_f = dir.path().join("out.fas");

    // 4 taxa, 6 columns:
    // col 0: A A A A -> 0/4 = 0.00 (keep at g=0.25)
    // col 1: T - T T -> 1/4 = 0.25 (keep at g=0.25)
    // col 2: C N C C -> 1/4 = 0.25 ('N' counts as gap -> keep at g=0.25, drop at g=0.0)
    // col 3: G R ? G -> 2/4 = 0.50 ('R' and '?' both count as gap -> drop at g=0.25)
    // col 4: u u u R -> 1/4 = 0.25 ('u' is valid nucleotide, 'R' is gap -> keep at g=0.25)
    // col 5: - - - - -> 4/4 = 1.00 (drop at g=0.25)
    // Note: Taxon4 is sparse/mostly gaps, must NEVER be dropped!
    let content = "\
>Taxon1
ATCGU-
>Taxon2
A-NRu-
>Taxon3
ATC?u-
>Taxon4
ATGGR-
";
    fs::write(&in_f, content).unwrap();

    let args = CliArgs {
        input: in_f.clone(),
        output: out_f.clone(),
        max_gap: 0.25,
        codon: false,
        threads: 2,
    };
    let stats = run(&args).unwrap();
    assert_eq!(stats.files_processed, 1);
    assert_eq!(stats.total_orig_sites, 6);
    assert_eq!(stats.total_kept_sites, 4); // cols 0, 1, 2, 4

    let (mat, _, _) = read_alignment(&out_f).unwrap();
    assert_eq!(mat.n_taxa, 4);
    assert_eq!(mat.n_sites, 4);
    assert_eq!(mat.row(0), b"ATCU");
    assert_eq!(mat.row(1), b"A-Nu");
    assert_eq!(mat.row(2), b"ATCu");
    assert_eq!(mat.row(3), b"ATGR");

    // Zero-tolerance max_gap = 0.0 keeps only col 0
    let out_zero = dir.path().join("out_zero.fas");
    run(&CliArgs {
        input: in_f,
        output: out_zero.clone(),
        max_gap: 0.0,
        codon: false,
        threads: 1,
    })
    .unwrap();
    let (mat0, _, _) = read_alignment(&out_zero).unwrap();
    assert_eq!(mat0.n_sites, 1);
    assert_eq!(mat0.row(0), b"A");
}

#[test]
fn test_all_gap_row_is_never_dropped() {
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("keep_rows.fas");
    let out_f = dir.path().join("out.fas");

    let content = "\
>Sp1
ATCG
>Sp2
ATCG
>Sp3
ATCG
>AllGapTaxon
----
";
    fs::write(&in_f, content).unwrap();

    run(&CliArgs {
        input: in_f,
        output: out_f.clone(),
        max_gap: 0.3, // 1/4 = 0.25 <= 0.3 -> all 4 columns kept, AllGapTaxon kept intact!
        codon: false,
        threads: 1,
    })
    .unwrap();

    let (mat, _, _) = read_alignment(&out_f).unwrap();
    assert_eq!(mat.n_taxa, 4);
    assert_eq!(mat.names[3], "AllGapTaxon");
    assert_eq!(mat.row(3), b"----");
}

#[test]
fn test_codon_mode_drops_entire_triplet_if_any_column_exceeds_threshold() {
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("cds.fas");
    let out_f = dir.path().join("cds_out.fas");

    // 4 taxa, 11 bp (3 complete codons + 2 trailing bases ignored):
    // Codon 0 (cols 0..3): gaps = [0, 1, 0] -> max 0.25 <= 0.25 -> KEEP
    // Codon 1 (cols 3..6): gaps = [0, 0, 2] -> col 5 has 2/4 = 0.50 > 0.25 -> DROP WHOLE CODON
    // Codon 2 (cols 6..9): gaps = [1, 0, 1] -> max 0.25 <= 0.25 -> KEEP
    // Trailing 2 bp (cols 9..11): ignored in codon mode
    let content = "\
>T1
ATGAAATTTAA
>T2
A-GAA-TT-AA
>T3
ATGAANTTAAA
>T4
ATGAAAGGGAA
";
    fs::write(&in_f, content).unwrap();

    let stats = run(&CliArgs {
        input: in_f,
        output: out_f.clone(),
        max_gap: 0.25,
        codon: true,
        threads: 2,
    })
    .unwrap();

    assert_eq!(stats.total_orig_sites, 11);
    assert_eq!(stats.total_kept_sites, 6); // Codon 0 + Codon 2

    let (mat, _, _) = read_alignment(&out_f).unwrap();
    assert_eq!(mat.n_sites % 3, 0);
    assert_eq!(mat.row(0), b"ATGTTT");
    assert_eq!(mat.row(1), b"A-GTT-");
    assert_eq!(mat.row(2), b"ATGTTA");
    assert_eq!(mat.row(3), b"ATGGGG");
}

#[test]
fn test_protein_auto_detection_and_codon_guard() {
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("prot.faa");
    let out_f = dir.path().join("prot_out.faa");

    // 4 taxa, 6 aa:
    // col 0: M M M M (0 gap) -> keep
    // col 1: N N N N (Asparagine is a standard AA -> 0 gap!) -> keep
    // col 2: E X E E ('X' is non-standard -> 1/4 = 0.25) -> keep at 0.25
    // col 3: L - X L (2/4 = 0.50) -> drop at 0.25
    // col 4: P Q F I (0 gap) -> keep
    // col 5: B Z * ? (all 4 are non-standard/gap -> 4/4 = 1.00) -> drop
    let content = "\
>P1
MNELPB
>P2
MNX-QZ
>P3
MNEXF*
>P4
MNELL?
";
    fs::write(&in_f, content).unwrap();

    // --codon on protein must return an error
    let err = run(&CliArgs {
        input: in_f.clone(),
        output: out_f.clone(),
        max_gap: 0.25,
        codon: true,
        threads: 1,
    });
    assert!(err.is_err());

    // Normal column gap mode on protein
    run(&CliArgs {
        input: in_f,
        output: out_f.clone(),
        max_gap: 0.25,
        codon: false,
        threads: 1,
    })
    .unwrap();

    let (mat, _, _) = read_alignment(&out_f).unwrap();
    assert_eq!(mat.n_sites, 4); // cols 0, 1, 2, 4
    assert_eq!(mat.row(0), b"MNEP");
    assert_eq!(mat.row(1), b"MNXQ");
    assert_eq!(mat.row(2), b"MNEF");
    assert_eq!(mat.row(3), b"MNEL");
}

#[test]
fn test_large_alignment_multi_tile_and_over_255_taxa_simd_batches() {
    // 300 taxa (> 255 u8 batch boundary!) x 24,000 sites (~7.2 MB, 3 L1 tiles, parallel execution)
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("huge.fas");
    let out_col = dir.path().join("huge_col.fas");
    let out_cod = dir.path().join("huge_cod.fas");

    let n_taxa = 300usize;
    let n_sites = 24_000usize;
    let mut fasta_bytes = Vec::with_capacity(n_taxa * (n_sites + 16));

    // Construct deterministic pattern:
    // Every column `c`: first `(c % 101)` taxa have '-', rest have 'A'/'C'/'G'/'T'.
    // With n_taxa = 300 and max_gap = 0.20 -> max_allowed_gaps = 60.
    // So column `c` passes iff `(c % 101) <= 60`.
    let bases = *b"ACGT";
    for r in 0..n_taxa {
        writeln!(fasta_bytes, ">Sp_{r:04}").unwrap();
        let mut row = Vec::with_capacity(n_sites);
        for c in 0..n_sites {
            let n_gaps = c % 101;
            if r < n_gaps {
                row.push(if (r + c) % 2 == 0 { b'-' } else { b'N' });
            } else {
                row.push(bases[(r + c) & 3]);
            }
        }
        fasta_bytes.extend_from_slice(&row);
        fasta_bytes.push(b'\n');
    }
    fs::write(&in_f, &fasta_bytes).unwrap();

    let expected_col_kept: Vec<usize> = (0..n_sites).filter(|&c| (c % 101) <= 60).collect();
    let mut expected_cod_kept = Vec::new();
    for k in 0..(n_sites / 3) {
        let c0 = 3 * k;
        if (c0 % 101) <= 60 && ((c0 + 1) % 101) <= 60 && ((c0 + 2) % 101) <= 60 {
            expected_cod_kept.extend_from_slice(&[c0, c0 + 1, c0 + 2]);
        }
    }

    let s1 = run(&CliArgs {
        input: in_f.clone(),
        output: out_col.clone(),
        max_gap: 0.20,
        codon: false,
        threads: 4,
    })
    .unwrap();
    assert_eq!(s1.total_kept_sites as usize, expected_col_kept.len());

    let s2 = run(&CliArgs {
        input: in_f,
        output: out_cod.clone(),
        max_gap: 0.20,
        codon: true,
        threads: 4,
    })
    .unwrap();
    assert_eq!(s2.total_kept_sites as usize, expected_cod_kept.len());

    let (m_col, _, _) = read_alignment(&out_col).unwrap();
    assert_eq!(m_col.n_taxa, n_taxa);
    assert_eq!(m_col.n_sites, expected_col_kept.len());
    // Spot-check taxon 0 and taxon 299
    for &r in &[0usize, 59, 60, 255, 299] {
        let expected_row: Vec<u8> = expected_col_kept
            .iter()
            .map(|&c| {
                if r < (c % 101) {
                    if (r + c) % 2 == 0 {
                        b'-'
                    } else {
                        b'N'
                    }
                } else {
                    bases[(r + c) & 3]
                }
            })
            .collect();
        assert_eq!(m_col.row(r), expected_row.as_slice());
    }
}

#[test]
fn test_batch_directory_multiformat_and_gzip_clean_output() {
    let dir = tempdir().unwrap();
    let in_dir = dir.path().join("in_loci");
    let out_dir = dir.path().join("out_loci");
    fs::create_dir_all(&in_dir).unwrap();

    // 1) Gzipped FASTA
    {
        let f = fs::File::create(in_dir.join("gene1.fas.gz")).unwrap();
        let mut gz = GzEncoder::new(f, Compression::fast());
        gz.write_all(b">A\nATG---\n>B\nATG-AA\n>C\nATGAAA\n")
            .unwrap();
        gz.finish().unwrap();
    }
    // 2) PHYLIP
    fs::write(in_dir.join("gene2.phy"), "3 4\nA  ATCG\nB  AT-G\nC  AT-G\n").unwrap();
    // 3) NEXUS
    fs::write(
        in_dir.join("gene3.nex"),
        "#NEXUS\nBEGIN DATA;\nDIMENSIONS NTAX=3 NCHAR=4;\nFORMAT DATATYPE=DNA GAP=- MISSING=?;\nMATRIX\nA ATCG\nB ATC-\nC ATC-\n;\nEND;\n",
    )
    .unwrap();

    let stats = run(&CliArgs {
        input: in_dir,
        output: out_dir.clone(),
        max_gap: 0.4, // allows at most 1 gap out of 3 taxa (1/3 = 0.333 <= 0.4; 2/3 = 0.667 > 0.4)
        codon: false,
        threads: 3,
    })
    .unwrap();
    assert_eq!(stats.files_processed, 3);

    // Ensure NO extra report/summary files pollute out_dir
    let mut names: Vec<String> = fs::read_dir(&out_dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().to_string())
        .collect();
    names.sort();
    assert_eq!(names, vec!["gene1.fas.gz", "gene2.phy", "gene3.nex"]);

    let (m1, _, gz1) = read_alignment(&out_dir.join("gene1.fas.gz")).unwrap();
    assert!(gz1);
    assert_eq!(m1.n_sites, 5); // col 3 has 2 gaps -> dropped; col 4 & 5 have 1 gap -> kept

    let (m2, _, _) = read_alignment(&out_dir.join("gene2.phy")).unwrap();
    assert_eq!(m2.n_sites, 3); // col 2 has 2 gaps -> dropped

    let (m3, _, _) = read_alignment(&out_dir.join("gene3.nex")).unwrap();
    assert_eq!(m3.n_sites, 3); // col 3 has 2 gaps -> dropped
}

#[test]
fn test_zero_copy_all_kept_and_all_dropped_edge_cases() {
    let dir = tempdir().unwrap();

    // 1) All columns pass -> zero-copy fast path preserves full alignment verbatim
    let in_keep = dir.path().join("all_keep.fas");
    let out_keep = dir.path().join("all_keep_out.fas");
    fs::write(&in_keep, ">S1\nATCGAT\n>S2\nATCGAT\n>S3\nATCGAT\n").unwrap();

    let stats_keep = run(&CliArgs {
        input: in_keep,
        output: out_keep.clone(),
        max_gap: 0.0,
        codon: true,
        threads: 1,
    })
    .unwrap();
    assert_eq!(stats_keep.total_orig_sites, 6);
    assert_eq!(stats_keep.total_kept_sites, 6);
    let (m_keep, _, _) = read_alignment(&out_keep).unwrap();
    assert_eq!(m_keep.n_taxa, 3);
    assert_eq!(m_keep.n_sites, 6);
    assert_eq!(m_keep.row(0), b"ATCGAT");

    // 2) All columns exceed threshold -> 0 kept sites, but all taxa are still preserved
    let in_drop = dir.path().join("all_drop.fas");
    let out_drop = dir.path().join("all_drop_out.fas");
    fs::write(&in_drop, ">S1\nA-N?\n>S2\n-T-?\n>S3\n--G-\n").unwrap();

    let stats_drop = run(&CliArgs {
        input: in_drop,
        output: out_drop.clone(),
        max_gap: 0.0,
        codon: false,
        threads: 1,
    })
    .unwrap();
    assert_eq!(stats_drop.total_orig_sites, 4);
    assert_eq!(stats_drop.total_kept_sites, 0);
    let raw_out = fs::read_to_string(&out_drop).unwrap();
    assert!(raw_out.contains(">S1\n"));
    assert!(raw_out.contains(">S2\n"));
    assert!(raw_out.contains(">S3\n"));
}

#[test]
fn test_cli_validation_guards() {
    let dir = tempdir().unwrap();
    let in_f = dir.path().join("sample.fas");
    let out_f = dir.path().join("out.fas");
    fs::write(&in_f, ">A\nATCG\n>B\nATCG\n").unwrap();

    // Same input and output path must fail
    let same_err = run(&CliArgs {
        input: in_f.clone(),
        output: in_f.clone(),
        max_gap: 0.2,
        codon: false,
        threads: 1,
    });
    assert!(same_err.is_err());

    // Out-of-range and NaN max_gap must fail
    for bad_gap in [-0.01, 1.01, f64::NAN, f64::INFINITY] {
        let err = run(&CliArgs {
            input: in_f.clone(),
            output: out_f.clone(),
            max_gap: bad_gap,
            codon: false,
            threads: 1,
        });
        assert!(err.is_err(), "expected error for max_gap={bad_gap}");
    }
}
