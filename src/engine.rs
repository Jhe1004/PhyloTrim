use crate::cli::CliArgs;
use crate::io::{
    is_supported_alignment_file, read_alignment, resolve_output_target, write_alignment,
    AlignmentFormat,
};
use crate::trim::trim_alignment;
use anyhow::{bail, Context, Result};
use rayon::prelude::*;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::time::Instant;

/// Summary statistics returned by a trimming run.
#[derive(Debug, Clone, Copy)]
pub struct RunStats {
    pub files_processed: usize,
    pub total_orig_sites: u64,
    pub total_kept_sites: u64,
    pub elapsed_secs: f64,
}

struct FileTask {
    input_path: PathBuf,
    size_bytes: u64,
}

/// Runs single-file or batch-directory trimming using the configured Rayon thread pool.
pub fn run(args: &CliArgs) -> Result<RunStats> {
    let start = Instant::now();
    args.validate()?;

    if !args.input.exists() {
        bail!("Input path does not exist: '{}'", args.input.display());
    }

    if args.threads > 0 {
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build()
            .context("Failed to initialize Rayon worker thread pool")?;
        pool.install(|| dispatch(args, start))
    } else {
        dispatch(args, start)
    }
}

fn dispatch(args: &CliArgs, start: Instant) -> Result<RunStats> {
    if args.input.is_file() {
        run_single_file(args, start)
    } else if args.input.is_dir() {
        run_directory(args, start)
    } else {
        bail!(
            "Input path '{}' is neither a regular file nor a directory",
            args.input.display()
        )
    }
}

fn path_looks_like_dir(path: &Path) -> bool {
    path.is_dir() || path.to_string_lossy().ends_with('/') || path.to_string_lossy().ends_with('\\')
}

fn infer_explicit_output_format(
    output_path: &Path,
    default_fmt: AlignmentFormat,
) -> AlignmentFormat {
    let Some(fname) = output_path.file_name().and_then(|n| n.to_str()) else {
        return default_fmt;
    };
    let without_gz = fname
        .strip_suffix(".gz")
        .or_else(|| fname.strip_suffix(".GZ"))
        .unwrap_or(fname);
    let Some((_, ext)) = without_gz.rsplit_once('.') else {
        return default_fmt;
    };
    match ext.to_ascii_lowercase().as_str() {
        "fa" | "fas" | "fasta" | "fna" | "ffn" | "faa" | "frn" | "aln" => AlignmentFormat::Fasta,
        "phy" | "phylip" => AlignmentFormat::Phylip,
        "nex" | "nexus" => AlignmentFormat::Nexus,
        _ => default_fmt,
    }
}

fn process_one_file(
    in_path: &Path,
    out_arg: &Path,
    out_is_dir: bool,
    max_gap: f64,
    codon: bool,
) -> Result<(u64, u64)> {
    let (matrix, in_fmt, was_gz) = read_alignment(in_path)?;
    if codon && !matrix.seq_kind.is_nucleotide() {
        bail!(
            "Codon mode (--codon) requires nucleotide alignments, but '{}' was detected as Protein",
            in_path.display()
        );
    }

    let outcome = trim_alignment(&matrix, max_gap, codon);
    let (target_path, gz_out) = resolve_output_target(in_path, was_gz, out_arg, out_is_dir);
    let out_fmt = if out_is_dir {
        in_fmt
    } else {
        infer_explicit_output_format(&target_path, in_fmt)
    };

    write_alignment(&target_path, &outcome.matrix, out_fmt, gz_out)?;
    Ok((outcome.orig_sites as u64, outcome.kept_sites as u64))
}

fn run_single_file(args: &CliArgs, start: Instant) -> Result<RunStats> {
    let out_is_dir = path_looks_like_dir(&args.output);
    if out_is_dir {
        fs::create_dir_all(&args.output).with_context(|| {
            format!(
                "Failed to create output directory '{}'",
                args.output.display()
            )
        })?;
    }

    let (orig_sites, kept_sites) = process_one_file(
        &args.input,
        &args.output,
        out_is_dir,
        args.max_gap,
        args.codon,
    )?;

    Ok(RunStats {
        files_processed: 1,
        total_orig_sites: orig_sites,
        total_kept_sites: kept_sites,
        elapsed_secs: start.elapsed().as_secs_f64(),
    })
}

fn run_directory(args: &CliArgs, start: Instant) -> Result<RunStats> {
    if args.output.is_file() {
        bail!(
            "Output '{}' is an existing file, but input '{}' is a directory",
            args.output.display(),
            args.input.display()
        );
    }
    fs::create_dir_all(&args.output).with_context(|| {
        format!(
            "Failed to create output directory '{}'",
            args.output.display()
        )
    })?;

    let mut tasks: Vec<FileTask> = Vec::new();
    for entry_res in fs::read_dir(&args.input)
        .with_context(|| format!("Failed to read input directory '{}'", args.input.display()))?
    {
        let entry = entry_res?;
        let path = entry.path();
        if path.is_file() && is_supported_alignment_file(&path) {
            let size_bytes = entry.metadata().map(|m| m.len()).unwrap_or(0);
            tasks.push(FileTask {
                input_path: path,
                size_bytes,
            });
        }
    }

    if tasks.is_empty() {
        bail!(
            "No supported alignment files (.fa/.fas/.fasta/.phy/.nex[.gz]) found in '{}'",
            args.input.display()
        );
    }

    // Largest-Processing-Time-First (LPT) ordering eliminates straggler tail latency when
    // huge alignments and small alignments coexist in the same directory.
    tasks.sort_unstable_by(|a, b| {
        b.size_bytes
            .cmp(&a.size_bytes)
            .then_with(|| a.input_path.cmp(&b.input_path))
    });

    let total_orig = AtomicU64::new(0);
    let total_kept = AtomicU64::new(0);

    tasks.par_iter().try_for_each(|task| -> Result<()> {
        let (orig, kept) = process_one_file(
            &task.input_path,
            &args.output,
            true,
            args.max_gap,
            args.codon,
        )?;
        total_orig.fetch_add(orig, Ordering::Relaxed);
        total_kept.fetch_add(kept, Ordering::Relaxed);
        Ok(())
    })?;

    Ok(RunStats {
        files_processed: tasks.len(),
        total_orig_sites: total_orig.load(Ordering::Relaxed),
        total_kept_sites: total_kept.load(Ordering::Relaxed),
        elapsed_secs: start.elapsed().as_secs_f64(),
    })
}
