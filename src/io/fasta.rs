use crate::io::strip_utf8_bom;
use crate::matrix::AlignmentMatrix;
use anyhow::{bail, Result};
use memchr::{memchr, memchr3};
use std::collections::HashSet;
use std::io::Write;

/// Parses a FASTA alignment using `memchr` line scanning directly into a single flat buffer.
/// When all sequences already have equal length (standard MSA), avoids any secondary reallocation
/// or `memcpy` completely.
pub fn parse_fasta(bytes: &[u8]) -> Result<AlignmentMatrix> {
    let data = strip_utf8_bom(bytes);
    let mut names: Vec<String> = Vec::new();
    let mut seen_names: HashSet<String> = HashSet::new();
    let mut flat: Vec<u8> = Vec::with_capacity(data.len());
    let mut spans: Vec<(usize, usize)> = Vec::new();

    let mut current_start = 0usize;
    let mut in_record = false;
    let mut pos = 0usize;
    let len = data.len();

    while pos < len {
        let line_end = memchr(b'\n', &data[pos..])
            .map(|rel| pos + rel)
            .unwrap_or(len);
        let mut line = &data[pos..line_end];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        pos = line_end + 1;

        let trimmed = trim_ascii_whitespace(line);
        if trimmed.is_empty() {
            continue;
        }

        if trimmed[0] == b'>' {
            if in_record {
                spans.push((current_start, flat.len() - current_start));
            }
            let header = std::str::from_utf8(&trimmed[1..]).unwrap_or("").trim();
            if header.is_empty() {
                bail!("FASTA record has empty header '>'");
            }
            let id = header
                .split_whitespace()
                .next()
                .unwrap_or(header)
                .to_string();
            if !seen_names.insert(id.clone()) {
                bail!("Duplicate taxon name '{id}' found in FASTA alignment");
            }
            names.push(id);
            current_start = flat.len();
            in_record = true;
        } else {
            if !in_record {
                bail!("FASTA sequence data appeared before first '>' header");
            }
            if memchr3(b' ', b'\t', b'\r', trimmed).is_none() {
                flat.extend_from_slice(trimmed);
            } else {
                for &b in trimmed {
                    if !b.is_ascii_whitespace() {
                        flat.push(b);
                    }
                }
            }
        }
    }

    if in_record {
        spans.push((current_start, flat.len() - current_start));
    }

    if names.is_empty() {
        bail!("No sequences found in FASTA input");
    }

    let max_len = spans.iter().map(|&(_, l)| l).max().unwrap_or(0);
    let min_len = spans.iter().map(|&(_, l)| l).min().unwrap_or(0);

    // Zero-copy fast path for standard equal-length alignments.
    if min_len == max_len {
        return AlignmentMatrix::from_flat(names, flat, max_len);
    }

    // Pad ragged sequences with '-' to `max_len`.
    let mut padded = Vec::with_capacity(names.len() * max_len);
    for &(st, slen) in &spans {
        padded.extend_from_slice(&flat[st..st + slen]);
        if slen < max_len {
            padded.resize(padded.len() + (max_len - slen), b'-');
        }
    }
    AlignmentMatrix::from_flat(names, padded, max_len)
}

/// Writes an `AlignmentMatrix` in unwrapped single-line-per-taxon FASTA format.
pub fn write_fasta<W: Write>(writer: &mut W, matrix: &AlignmentMatrix) -> Result<()> {
    for i in 0..matrix.n_taxa {
        writer.write_all(b">")?;
        writer.write_all(matrix.names[i].as_bytes())?;
        writer.write_all(b"\n")?;
        writer.write_all(matrix.row(i))?;
        writer.write_all(b"\n")?;
    }
    Ok(())
}

fn trim_ascii_whitespace(mut s: &[u8]) -> &[u8] {
    while let Some((&b, rest)) = s.split_first() {
        if b.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    while let Some((&b, rest)) = s.split_last() {
        if b.is_ascii_whitespace() {
            s = rest;
        } else {
            break;
        }
    }
    s
}
