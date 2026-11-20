use anyhow::Result;
use clap::Parser;
use phylotrim::cli::CliArgs;
use phylotrim::engine::run;

fn main() -> Result<()> {
    let args = CliArgs::parse();
    let stats = run(&args)?;

    let pct = if stats.total_orig_sites == 0 {
        0.0
    } else {
        (stats.total_kept_sites as f64 / stats.total_orig_sites as f64) * 100.0
    };

    eprintln!(
        "[phylotrim] Done: {} file(s), {} -> {} sites ({:.1}% kept) in {:.3}s",
        stats.files_processed,
        stats.total_orig_sites,
        stats.total_kept_sites,
        pct,
        stats.elapsed_secs
    );
    Ok(())
}
