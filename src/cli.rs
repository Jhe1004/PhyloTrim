use anyhow::{bail, Result};
use clap::Parser;
use std::path::PathBuf;

/// Minimalist, ultra-fast column-gap and codon MSA trimmer.
#[derive(Debug, Clone, Parser)]
#[command(
    name = "phylotrim",
    version,
    about = "Minimalist, ultra-fast column-gap and 3-bp codon MSA trimmer in Rust"
)]
pub struct CliArgs {
    /// Input alignment file or directory (FASTA / PHYLIP / NEXUS, optional .gz)
    #[arg(short = 'i', long = "input", value_name = "PATH")]
    pub input: PathBuf,

    /// Output alignment file or directory (optional .gz)
    #[arg(short = 'o', long = "output", value_name = "PATH")]
    pub output: PathBuf,

    /// Maximum allowed gap / non-standard symbol fraction per column [0.0..=1.0]
    #[arg(
        short = 'g',
        long = "max-gap",
        value_name = "FLOAT",
        default_value_t = 0.2
    )]
    pub max_gap: f64,

    /// Trim by 3-bp codons (drop entire codon if any of its 3 columns exceeds --max-gap)
    #[arg(short = 'c', long = "codon")]
    pub codon: bool,

    /// Number of worker threads (0 = use all available CPU cores)
    #[arg(short = 't', long = "threads", value_name = "INT", default_value_t = 0)]
    pub threads: usize,
}

impl CliArgs {
    pub fn validate(&self) -> Result<()> {
        if !self.max_gap.is_finite() || !(0.0..=1.0).contains(&self.max_gap) {
            bail!(
                "Invalid --max-gap {}: must be a number between 0.0 and 1.0",
                self.max_gap
            );
        }
        if let (Ok(in_canon), Ok(out_canon)) =
            (self.input.canonicalize(), self.output.canonicalize())
        {
            if in_canon == out_canon {
                bail!(
                    "Input and output paths point to the same location ('{}')",
                    in_canon.display()
                );
            }
        }
        Ok(())
    }
}
