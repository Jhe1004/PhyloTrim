use crate::io::strip_utf8_bom;
use crate::matrix::AlignmentMatrix;
use anyhow::{bail, Result};
use std::collections::HashSet;
use std::io::Write;

/// Parses a relaxed PHYLIP alignment (supports both Sequential and Interleaved layouts).
pub fn parse_phylip(bytes: &[u8]) -> Result<AlignmentMatrix> {
    let text = std::str::from_utf8(strip_utf8_bom(bytes))
        .map_err(|_| anyhow::anyhow!("PHYLIP file contains invalid UTF-8"))?;

    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
    let header_line = lines
        .next()
        .ok_or_else(|| anyhow::anyhow!("PHYLIP file is empty"))?;

    let header_tokens: Vec<&str> = header_line.split_whitespace().collect();
    if header_tokens.len() < 2 {
        bail!("Invalid PHYLIP header line: '{header_line}'");
    }
    let n_taxa: usize = header_tokens[0]
        .parse()
        .map_err(|_| anyhow::anyhow!("Invalid PHYLIP taxon count '{}'", header_tokens[0]))?;
    let n_sites: usize = header_tokens[1]
        .parse()
        .map_err(|_| anyhow::anyhow!("Invalid PHYLIP site count '{}'", header_tokens[1]))?;

    if n_taxa == 0 || n_sites == 0 {
        bail!("PHYLIP dimensions must be positive, got {n_taxa} x {n_sites}");
    }

    let body_lines: Vec<&str> = lines.collect();
    if body_lines.len() < n_taxa {
        bail!(
            "PHYLIP body has {} non-empty lines, fewer than declared {n_taxa} taxa",
            body_lines.len()
        );
    }

    if body_lines.len() == n_taxa {
        return parse_single_block(&body_lines, n_taxa, n_sites);
    }

    if let Ok(mat) = try_parse_interleaved(&body_lines, n_taxa, n_sites) {
        return Ok(mat);
    }
    if let Ok(mat) = try_parse_sequential(&body_lines, n_taxa, n_sites) {
        return Ok(mat);
    }

    bail!("Failed to parse PHYLIP alignment ({n_taxa} taxa x {n_sites} sites) as either interleaved or sequential")
}

fn parse_single_block(lines: &[&str], n_taxa: usize, n_sites: usize) -> Result<AlignmentMatrix> {
    let mut names = Vec::with_capacity(n_taxa);
    let mut flat = Vec::with_capacity(n_taxa * n_sites);
    let mut seen = HashSet::with_capacity(n_taxa);

    for &line in lines {
        let (name, frag) = split_phylip_name_and_seq(line)?;
        if !seen.insert(name.clone()) {
            bail!("Duplicate taxon name '{name}' in PHYLIP alignment");
        }
        let start = flat.len();
        for b in frag.bytes() {
            if !b.is_ascii_whitespace() {
                flat.push(b);
            }
        }
        let added = flat.len() - start;
        if added != n_sites {
            bail!("PHYLIP taxon '{name}' has length {added}, expected {n_sites}");
        }
        names.push(name);
    }

    AlignmentMatrix::from_flat(names, flat, n_sites)
}

fn try_parse_interleaved(lines: &[&str], n_taxa: usize, n_sites: usize) -> Result<AlignmentMatrix> {
    if !lines.len().is_multiple_of(n_taxa) {
        bail!("Line count is not a multiple of n_taxa for interleaved PHYLIP");
    }

    let mut names = Vec::with_capacity(n_taxa);
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n_taxa);
    let mut seen = HashSet::with_capacity(n_taxa);

    for &line in &lines[..n_taxa] {
        let (name, frag) = split_phylip_name_and_seq(line)?;
        if !seen.insert(name.clone()) {
            bail!("Duplicate taxon name '{name}'");
        }
        let mut s = Vec::with_capacity(n_sites);
        for b in frag.bytes() {
            if !b.is_ascii_whitespace() {
                s.push(b);
            }
        }
        names.push(name);
        seqs.push(s);
    }

    for chunk in lines[n_taxa..].chunks_exact(n_taxa) {
        for (i, &line) in chunk.iter().enumerate() {
            for b in line.bytes() {
                if !b.is_ascii_whitespace() {
                    seqs[i].push(b);
                }
            }
        }
    }

    for (name, seq) in names.iter().zip(seqs.iter()) {
        if seq.len() != n_sites {
            bail!(
                "Interleaved PHYLIP taxon '{name}' length {} != {n_sites}",
                seq.len()
            );
        }
    }

    AlignmentMatrix::from_records(names, seqs)
}

fn try_parse_sequential(lines: &[&str], n_taxa: usize, n_sites: usize) -> Result<AlignmentMatrix> {
    let mut names = Vec::with_capacity(n_taxa);
    let mut seqs: Vec<Vec<u8>> = Vec::with_capacity(n_taxa);
    let mut seen = HashSet::with_capacity(n_taxa);
    let mut idx = 0usize;

    for _ in 0..n_taxa {
        if idx >= lines.len() {
            bail!("Unexpected end of sequential PHYLIP");
        }
        let (name, frag) = split_phylip_name_and_seq(lines[idx])?;
        idx += 1;
        if !seen.insert(name.clone()) {
            bail!("Duplicate taxon name '{name}'");
        }
        let mut seq = Vec::with_capacity(n_sites);
        for b in frag.bytes() {
            if !b.is_ascii_whitespace() {
                seq.push(b);
            }
        }
        while seq.len() < n_sites {
            if idx >= lines.len() {
                bail!("Unexpected end of sequence for '{name}'");
            }
            for b in lines[idx].bytes() {
                if !b.is_ascii_whitespace() {
                    seq.push(b);
                }
            }
            idx += 1;
        }
        if seq.len() != n_sites {
            bail!(
                "Sequential PHYLIP taxon '{name}' length {} != {n_sites}",
                seq.len()
            );
        }
        names.push(name);
        seqs.push(seq);
    }

    if idx != lines.len() {
        bail!("Trailing lines left over after sequential PHYLIP parsing");
    }

    AlignmentMatrix::from_records(names, seqs)
}

fn split_phylip_name_and_seq(line: &str) -> Result<(String, &str)> {
    let mut parts = line.splitn(2, char::is_whitespace);
    let name = parts.next().unwrap_or("").trim();
    let seq_part = parts.next().unwrap_or("").trim();
    if name.is_empty() || seq_part.is_empty() {
        bail!("Invalid PHYLIP taxon line: '{line}'");
    }
    Ok((name.to_string(), seq_part))
}

/// Writes an `AlignmentMatrix` in relaxed sequential PHYLIP format.
pub fn write_phylip<W: Write>(writer: &mut W, matrix: &AlignmentMatrix) -> Result<()> {
    writeln!(writer, "{} {}", matrix.n_taxa, matrix.n_sites)?;
    let max_name_len = matrix
        .names
        .iter()
        .map(|n| n.len())
        .max()
        .unwrap_or(10)
        .max(10);
    let spaces = vec![b' '; max_name_len + 2];
    for i in 0..matrix.n_taxa {
        let name = &matrix.names[i];
        writer.write_all(name.as_bytes())?;
        let pad = (max_name_len + 2).saturating_sub(name.len());
        writer.write_all(&spaces[..pad])?;
        writer.write_all(matrix.row(i))?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}
