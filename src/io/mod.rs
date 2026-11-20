pub mod fasta;
pub mod nexus;
pub mod phylip;

use crate::matrix::AlignmentMatrix;
use anyhow::{bail, Context, Result};
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::fs::{self, File};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

static TMP_COUNTER: AtomicU64 = AtomicU64::new(0);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentFormat {
    Fasta,
    Phylip,
    Nexus,
}

impl AlignmentFormat {
    pub fn default_extension(self) -> &'static str {
        match self {
            AlignmentFormat::Fasta => "fas",
            AlignmentFormat::Phylip => "phy",
            AlignmentFormat::Nexus => "nex",
        }
    }
}

/// Detects format from raw alignment text bytes (FASTA, PHYLIP, or NEXUS).
pub fn detect_format(bytes: &[u8]) -> Result<AlignmentFormat> {
    let text = strip_utf8_bom(bytes);
    let trimmed = trim_leading_ascii_whitespace(text);
    if trimmed.is_empty() {
        bail!("Alignment input is empty");
    }

    if trimmed[0] == b'>' {
        return Ok(AlignmentFormat::Fasta);
    }

    if trimmed.len() >= 6 && trimmed[..6].eq_ignore_ascii_case(b"#NEXUS") {
        return Ok(AlignmentFormat::Nexus);
    }

    if let Some(first_line) = trimmed.split(|&b| b == b'\n').next() {
        let line_str = String::from_utf8_lossy(first_line);
        let tokens: Vec<&str> = line_str.split_whitespace().collect();
        if tokens.len() >= 2
            && tokens[0].parse::<usize>().is_ok_and(|n| n > 0)
            && tokens[1].parse::<usize>().is_ok_and(|n| n > 0)
        {
            return Ok(AlignmentFormat::Phylip);
        }
    }

    bail!("Unable to detect alignment format (expected FASTA, PHYLIP, or NEXUS)")
}

/// Reads an alignment file (transparently handling `.gz` compression) and returns
/// `(matrix, format, was_gzipped)`.
pub fn read_alignment(path: &Path) -> Result<(AlignmentMatrix, AlignmentFormat, bool)> {
    let (raw_bytes, was_gz) = read_file_maybe_gz(path)
        .with_context(|| format!("Failed to read alignment file '{}'", path.display()))?;
    let format = detect_format(&raw_bytes)
        .with_context(|| format!("Failed to detect format of '{}'", path.display()))?;

    let matrix = match format {
        AlignmentFormat::Fasta => fasta::parse_fasta(&raw_bytes),
        AlignmentFormat::Phylip => phylip::parse_phylip(&raw_bytes),
        AlignmentFormat::Nexus => nexus::parse_nexus(&raw_bytes),
    }
    .with_context(|| format!("Failed to parse alignment '{}'", path.display()))?;

    Ok((matrix, format, was_gz))
}

/// Writes an alignment matrix atomically, respecting output format and optional gzip compression.
pub fn write_alignment(
    path: &Path,
    matrix: &AlignmentMatrix,
    format: AlignmentFormat,
    gzip: bool,
) -> Result<()> {
    let approx_bytes = matrix
        .n_taxa
        .saturating_mul(matrix.n_sites.saturating_add(32));
    let buf_cap = approx_bytes.clamp(8 * 1024, 256 * 1024);

    write_atomic(path, buf_cap, |writer| {
        if gzip {
            let mut gz = GzEncoder::new(writer, Compression::fast());
            write_format_to(&mut gz, matrix, format)?;
            gz.finish()?;
        } else {
            write_format_to(writer, matrix, format)?;
        }
        Ok(())
    })
}

fn write_format_to<W: Write>(
    writer: &mut W,
    matrix: &AlignmentMatrix,
    format: AlignmentFormat,
) -> Result<()> {
    match format {
        AlignmentFormat::Fasta => fasta::write_fasta(writer, matrix),
        AlignmentFormat::Phylip => phylip::write_phylip(writer, matrix),
        AlignmentFormat::Nexus => nexus::write_nexus(writer, matrix),
    }
}

/// Writes bytes atomically to `path` via a unique hidden temp file in the same directory.
pub fn write_atomic<F>(path: &Path, buf_capacity: usize, write_fn: F) -> Result<()>
where
    F: FnOnce(&mut BufWriter<File>) -> Result<()>,
{
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    if !parent.exists() {
        fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create parent directory '{}'", parent.display()))?;
    }

    let file_name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("output");
    let seq = TMP_COUNTER.fetch_add(1, Ordering::Relaxed);
    let tmp_path = parent.join(format!(
        ".phylotrim.tmp.{}.{}.{}",
        std::process::id(),
        seq,
        file_name
    ));

    let write_result = (|| -> Result<()> {
        let file = File::create(&tmp_path)
            .with_context(|| format!("Failed to create temporary file '{}'", tmp_path.display()))?;
        let mut buf = BufWriter::with_capacity(buf_capacity, file);
        write_fn(&mut buf)?;
        buf.flush()?;
        Ok(())
    })();

    if let Err(err) = write_result {
        let _ = fs::remove_file(&tmp_path);
        return Err(err);
    }

    if let Err(err) = fs::rename(&tmp_path, path) {
        let _ = fs::remove_file(&tmp_path);
        return Err(err).with_context(|| {
            format!(
                "Failed to atomically move '{}' to '{}'",
                tmp_path.display(),
                path.display()
            )
        });
    }

    Ok(())
}

/// Resolves output path and gzip flag for single-file or batch-directory outputs.
pub fn resolve_output_target(
    input_path: &Path,
    input_was_gz: bool,
    output_arg: &Path,
    output_is_dir: bool,
) -> (PathBuf, bool) {
    if output_is_dir {
        let fname = input_path
            .file_name()
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from("output.fas"));
        let gz = input_was_gz
            || fname
                .extension()
                .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
        (output_arg.join(fname), gz)
    } else {
        let gz = output_arg
            .extension()
            .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
        (output_arg.to_path_buf(), gz)
    }
}

/// Checks whether a directory entry has a supported alignment extension (with optional `.gz`).
pub fn is_supported_alignment_file(path: &Path) -> bool {
    let Some(file_name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    if file_name.starts_with('.') {
        return false;
    }
    let without_gz = file_name
        .strip_suffix(".gz")
        .or_else(|| file_name.strip_suffix(".GZ"))
        .unwrap_or(file_name);

    let Some((_, ext)) = without_gz.rsplit_once('.') else {
        return false;
    };

    matches!(
        ext.to_ascii_lowercase().as_str(),
        "fa" | "fas"
            | "fasta"
            | "fna"
            | "ffn"
            | "faa"
            | "frn"
            | "aln"
            | "phy"
            | "phylip"
            | "nex"
            | "nexus"
    )
}

fn read_file_maybe_gz(path: &Path) -> Result<(Vec<u8>, bool)> {
    let mut file = File::open(path)?;
    let hint = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
    let mut raw = Vec::with_capacity(hint);
    file.read_to_end(&mut raw)?;

    let is_gz_magic = raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b;
    let is_gz_ext = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));

    if is_gz_magic || is_gz_ext {
        let mut decoder = MultiGzDecoder::new(&raw[..]);
        let mut decompressed = Vec::with_capacity(raw.len().saturating_mul(3));
        decoder
            .read_to_end(&mut decompressed)
            .context("Failed to decompress gzip stream")?;
        Ok((decompressed, true))
    } else {
        Ok((raw, false))
    }
}

pub(crate) fn strip_utf8_bom(bytes: &[u8]) -> &[u8] {
    bytes.strip_prefix(&[0xEF, 0xBB, 0xBF]).unwrap_or(bytes)
}

fn trim_leading_ascii_whitespace(mut bytes: &[u8]) -> &[u8] {
    while let Some((&b, rest)) = bytes.split_first() {
        if b.is_ascii_whitespace() {
            bytes = rest;
        } else {
            break;
        }
    }
    bytes
}
