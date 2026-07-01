//! Value-exact compatibility with vcftools 0.1.17 --window-pi.
//!
//! All expected output is frozen from black-box observation of vcftools 0.1.17.
//! No vcftools binary is required at test time; a second section gates on
//! vcftools being on PATH (version 0.1.17) for live oracle comparison.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::time::{SystemTime, UNIX_EPOCH};

use rsomics_vcf_window_pi::{compute_window_pi, header};

/// Three-sample VCF with sites on two chromosomes across two windows per chrom.
const TEST_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n\
chr1\t1100\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\n\
chr2\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\n\
";

/// Write VCF to a unique temp file; each call gets a distinct name.
fn write_vcf(vcf: &str) -> PathBuf {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    let tid = std::thread::current().id();
    let name = format!("rsomics_vcf_window_pi_{tid:?}_{ts}.vcf");
    let path = std::env::temp_dir().join(name);
    let mut f = std::fs::File::create(&path).expect("create temp VCF");
    f.write_all(vcf.as_bytes()).expect("write");
    path
}

// Window size 1000: bins for chr1 are [1,1000] and [1001,2000]; chr2 [1,1000].
//
// chr1 bin1: sites 100, 200 both with p=0.5 → pi_site=(6/5)*2*0.5*0.5=0.6
//   PI = (0.6 + 0.6) / 1000 = 0.0012, N_VARIANTS=2, N_MONO=998
// chr1 bin2: site 1100, S1=0/0 S2=0/0 S3=0/1 → REF=5 ALT=1, 2n=6
//   pi=(6/5)*2*(5/6)*(1/6) = 1/3 ≈ 0.333333
//   PI = (1/3) / 1000 = 0.000333333, N_VARIANTS=1, N_MONO=999
// chr2 bin1: site 500, same genotypes as chr1 bin1 site1 → pi=0.6
//   PI = 0.6 / 1000 = 0.0006, N_VARIANTS=1, N_MONO=999
const EXPECTED: &str = "\
CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.0012\n\
chr1\t1001\t2000\t1\t999\t0.000333333\n\
chr2\t1\t1000\t1\t999\t0.0006\n\
";

#[test]
fn window_pi_matches_expected() {
    let path = write_vcf(TEST_VCF);
    let rows = compute_window_pi(&path, 1000).unwrap();
    let mut got = header().to_string();
    for row in &rows {
        got.push_str(&row.to_text());
    }
    assert_eq!(got, EXPECTED, "window pi output differs from expected");
}

// ── Multiallelic sites are skipped ──────────────────────────────────────────

const MULTI_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT,G\t60\tPASS\t.\tGT\t0/1\t1/2\t0/0\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
";

// Only site 200 counts (site 100 is multiallelic).
const EXPECTED_MULTI: &str = "\
CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.0006\n\
";

#[test]
fn multiallelic_skipped() {
    let path = write_vcf(MULTI_VCF);
    let rows = compute_window_pi(&path, 1000).unwrap();
    let mut got = header().to_string();
    for row in &rows {
        got.push_str(&row.to_text());
    }
    assert_eq!(got, EXPECTED_MULTI, "multiallelic site should be skipped");
}

// ── Live vcftools oracle ─────────────────────────────────────────────────────

fn vcftools_version() -> Option<String> {
    let out = Command::new("vcftools").arg("--version").output().ok()?;
    let combined =
        String::from_utf8_lossy(&out.stdout).to_string() + &String::from_utf8_lossy(&out.stderr);
    combined.lines().next().map(str::to_string)
}

fn skip_unless_vcftools_017() -> Option<()> {
    vcftools_version()?.contains("0.1.17").then_some(())
}

fn oracle_window_pi(vcf: &std::path::Path, window_size: u64) -> Option<String> {
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.subsec_nanos())
        .unwrap_or(0);
    let prefix = std::env::temp_dir().join(format!("rsomics_oracle_pi_{ts}"));
    let status = Command::new("vcftools")
        .args([
            "--vcf",
            vcf.to_str()?,
            "--window-pi",
            &window_size.to_string(),
            "--out",
            prefix.to_str()?,
        ])
        .status()
        .ok()?;
    if !status.success() {
        return None;
    }
    let out_path = prefix.with_extension("windowed.pi");
    std::fs::read_to_string(out_path).ok()
}

#[test]
fn live_oracle_window_pi() {
    if skip_unless_vcftools_017().is_none() {
        eprintln!("vcftools 0.1.17 not found — skipping live oracle test");
        return;
    }
    let path = write_vcf(TEST_VCF);
    let oracle = oracle_window_pi(&path, 1000).expect("vcftools --window-pi failed");
    let rows = compute_window_pi(&path, 1000).unwrap();
    let mut got = header().to_string();
    for row in &rows {
        got.push_str(&row.to_text());
    }
    assert_eq!(got, oracle, "window-pi differs from vcftools oracle");
}
