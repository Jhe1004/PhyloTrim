use crate::io::strip_utf8_bom;
use crate::matrix::AlignmentMatrix;
use anyhow::{bail, Result};
use std::collections::HashMap;
use std::io::Write;

/// Parses a NEXUS alignment (`BEGIN DATA` or `BEGIN CHARACTERS` block with `MATRIX ... ;`).
/// Supports single-quoted taxon names, `[...]` comments, and interleaved blocks.
pub fn parse_nexus(bytes: &[u8]) -> Result<AlignmentMatrix> {
    let text = std::str::from_utf8(strip_utf8_bom(bytes))
        .map_err(|_| anyhow::anyhow!("NEXUS file contains invalid UTF-8"))?;

    let cleaned = strip_nexus_comments(text);
    let upper = cleaned.to_ascii_uppercase();

    let matrix_kw_pos = find_keyword_boundary(&upper, "MATRIX")
        .ok_or_else(|| anyhow::anyhow!("No MATRIX section found in NEXUS file"))?;
    let after_matrix = &cleaned[matrix_kw_pos + "MATRIX".len()..];

    let semicolon_pos = find_unquoted_semicolon(after_matrix)
        .ok_or_else(|| anyhow::anyhow!("NEXUS MATRIX block is missing terminating ';'"))?;
    let matrix_body = &after_matrix[..semicolon_pos];

    let mut names: Vec<String> = Vec::new();
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut name_to_idx: HashMap<String, usize> = HashMap::new();

    for raw_line in matrix_body.lines() {
        let line = raw_line.trim();
        if line.is_empty() {
            continue;
        }
        let (taxon, seq_str) = parse_nexus_matrix_line(line)?;
        let seq_bytes: Vec<u8> = seq_str
            .bytes()
            .filter(|b| !b.is_ascii_whitespace())
            .collect();
        if seq_bytes.is_empty() {
            continue;
        }

        if let Some(&idx) = name_to_idx.get(&taxon) {
            seqs[idx].extend_from_slice(&seq_bytes);
        } else {
            let idx = names.len();
            name_to_idx.insert(taxon.clone(), idx);
            names.push(taxon);
            seqs.push(seq_bytes);
        }
    }

    if names.is_empty() {
        bail!("NEXUS MATRIX block contains no sequences");
    }

    AlignmentMatrix::from_records(names, seqs)
}

/// Writes an `AlignmentMatrix` in standard NEXUS format.
pub fn write_nexus<W: Write>(writer: &mut W, matrix: &AlignmentMatrix) -> Result<()> {
    writeln!(writer, "#NEXUS")?;
    writeln!(writer, "BEGIN DATA;")?;
    writeln!(
        writer,
        "\tDIMENSIONS NTAX={} NCHAR={};",
        matrix.n_taxa, matrix.n_sites
    )?;
    writeln!(
        writer,
        "\tFORMAT DATATYPE={} GAP=- MISSING=?;",
        matrix.seq_kind.nexus_datatype()
    )?;
    writeln!(writer, "MATRIX")?;
    for i in 0..matrix.n_taxa {
        let formatted_name = format_nexus_name(&matrix.names[i]);
        write!(writer, "\t{formatted_name} ")?;
        writer.write_all(matrix.row(i))?;
        writer.write_all(b"\n")?;
    }
    writeln!(writer, ";")?;
    writeln!(writer, "END;")?;
    Ok(())
}

fn format_nexus_name(name: &str) -> String {
    if name
        .chars()
        .any(|c| c.is_whitespace() || "()[];,=".contains(c))
    {
        let escaped = name.replace('\'', "''");
        format!("'{escaped}'")
    } else {
        name.to_string()
    }
}

fn strip_nexus_comments(input: &str) -> String {
    let mut out = String::with_capacity(input.len());
    let mut in_quotes = false;
    let mut comment_depth = 0usize;

    for ch in input.chars() {
        if comment_depth == 0 {
            if ch == '\'' {
                in_quotes = !in_quotes;
                out.push(ch);
            } else if !in_quotes && ch == '[' {
                comment_depth = 1;
            } else {
                out.push(ch);
            }
        } else if ch == '[' {
            comment_depth += 1;
        } else if ch == ']' {
            comment_depth -= 1;
        } else if ch == '\n' {
            out.push('\n');
        }
    }
    out
}

fn find_keyword_boundary(upper: &str, keyword: &str) -> Option<usize> {
    let mut offset = 0usize;
    while let Some(pos) = upper[offset..].find(keyword) {
        let abs = offset + pos;
        let before_ok = abs == 0
            || !upper.as_bytes()[abs - 1].is_ascii_alphanumeric()
                && upper.as_bytes()[abs - 1] != b'_';
        let end = abs + keyword.len();
        let after_ok = end == upper.len()
            || !upper.as_bytes()[end].is_ascii_alphanumeric() && upper.as_bytes()[end] != b'_';
        if before_ok && after_ok {
            return Some(abs);
        }
        offset = abs + keyword.len();
    }
    None
}

fn find_unquoted_semicolon(s: &str) -> Option<usize> {
    let mut in_quotes = false;
    for (idx, ch) in s.char_indices() {
        if ch == '\'' {
            in_quotes = !in_quotes;
        } else if !in_quotes && ch == ';' {
            return Some(idx);
        }
    }
    None
}

fn parse_nexus_matrix_line(line: &str) -> Result<(String, &str)> {
    if let Some(after_open) = line.strip_prefix('\'') {
        let mut name = String::new();
        let mut chars = after_open.char_indices();
        while let Some((idx, ch)) = chars.next() {
            if ch == '\'' {
                if let Some((_, next_ch)) = chars.clone().next() {
                    if next_ch == '\'' {
                        name.push('\'');
                        chars.next();
                        continue;
                    }
                }
                let rest = after_open[idx + 1..].trim();
                return Ok((name, rest));
            }
            name.push(ch);
        }
        bail!("Unterminated single-quoted taxon name in NEXUS MATRIX: {line}");
    } else {
        let mut parts = line.splitn(2, char::is_whitespace);
        let name = parts.next().unwrap_or("").trim();
        let seq = parts.next().unwrap_or("").trim();
        if name.is_empty() {
            bail!("Empty taxon name in NEXUS MATRIX");
        }
        Ok((name.to_string(), seq))
    }
}
