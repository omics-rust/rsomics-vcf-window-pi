use std::io::{Write, stdout};
use std::path::PathBuf;
use std::process::ExitCode;

use clap::Parser;
use rsomics_common::{CommonFlags, ToolMeta};
use serde_json::json;

use rsomics_vcf_window_pi::{compute_window_pi, header};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

/// Nucleotide diversity (π) in fixed-width windows from VCF.
///
/// Reads a VCF (plain or gzip), partitions each chromosome into non-overlapping
/// windows of --window-size bp, and reports per-window π — byte-identical to
/// vcftools --window-pi.
#[derive(Parser, Debug)]
#[command(name = "rsomics-vcf-window-pi", version, about, long_about = None)]
pub struct Cli {
    /// Input VCF (plain or .gz); required.
    #[arg(value_name = "VCF")]
    pub vcf: PathBuf,

    /// Window size in bp.
    #[arg(long, default_value = "10000", value_name = "BP")]
    pub window_size: u64,

    #[command(flatten)]
    pub common: CommonFlags,
}

impl Cli {
    pub fn run(self) -> ExitCode {
        match self.execute() {
            Ok(()) => ExitCode::SUCCESS,
            Err(e) => {
                eprintln!("error: {e}");
                ExitCode::FAILURE
            }
        }
    }

    fn execute(self) -> rsomics_common::Result<()> {
        let rows = compute_window_pi(&self.vcf, self.window_size)?;
        if self.common.json {
            let env = json!({
                "schema_version": rsomics_common::SCHEMA_VERSION,
                "tool": META.name,
                "tool_version": META.version,
                "status": "ok",
                "result": { "rows": rows },
            });
            println!(
                "{}",
                serde_json::to_string(&env)
                    .map_err(|e| { rsomics_common::RsomicsError::InvalidInput(e.to_string()) })?
            );
        } else {
            let mut out = stdout().lock();
            out.write_all(header().as_bytes())
                .map_err(rsomics_common::RsomicsError::Io)?;
            for row in &rows {
                out.write_all(row.to_text().as_bytes())
                    .map_err(rsomics_common::RsomicsError::Io)?;
            }
        }
        Ok(())
    }
}

#[test]
fn cli_debug_assert() {
    use clap::CommandFactory;
    Cli::command().debug_assert();
}
