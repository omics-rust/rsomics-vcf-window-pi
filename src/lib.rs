//! Nucleotide diversity (π) in fixed-width windows from VCF.
//!
//! Implements vcftools `--window-pi` exactly: non-overlapping bins of `window_size` bp
//! per chromosome, Tajima's estimator `π_site = (2n/(2n-1)) * 2p(1-p)` per biallelic site,
//! window π = Σ(π_site) / window_size.
//!
//! Only biallelic sites with at least one called genotype contribute to π.
//! Monomorphic sites (all REF or all ALT) contribute 0 to the sum (2p(1-p) = 0).

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

// ── Output row ───────────────────────────────────────────────────────────────

#[derive(Debug, Clone, Serialize)]
pub struct WindowPiRow {
    pub chrom: String,
    pub bin_start: u64,
    pub bin_end: u64,
    pub n_variants: u64,
    pub n_monomorphic: u64,
    pub pi: f64,
}

impl WindowPiRow {
    /// Render as the tab-separated line vcftools emits.
    pub fn to_text(&self) -> String {
        format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            self.chrom,
            self.bin_start,
            self.bin_end,
            self.n_variants,
            self.n_monomorphic,
            format_g(self.pi),
        )
    }
}

/// Format a float like vcftools `%g` (6 significant figures, no trailing zeros/dot).
pub fn format_g(x: f64) -> String {
    if x == 0.0 {
        return "0".to_string();
    }
    let mag = x.abs().log10().floor() as i32;
    let dec = (5 - mag).max(0) as usize;
    let s = format!("{x:.dec$}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    s.to_string()
}

// ── Window accumulator ───────────────────────────────────────────────────────

struct WindowAcc {
    window_size: u64,
    /// Current chromosome being accumulated.
    cur_chrom: String,
    /// Completed rows (all chromosomes, in order).
    rows: Vec<WindowPiRow>,
    /// Partially-accumulated bins for cur_chrom: bin_index → (pi_sum, n_variants).
    bins: HashMap<u64, (f64, u64)>,
    /// Maximum bin index seen for cur_chrom (to emit zero-variant bins).
    max_bin: Option<u64>,
}

impl WindowAcc {
    fn new(window_size: u64) -> Self {
        Self {
            window_size,
            cur_chrom: String::new(),
            rows: Vec::new(),
            bins: HashMap::new(),
            max_bin: None,
        }
    }

    fn bin_of(&self, pos: u64) -> u64 {
        // vcftools bins: bin 0 covers positions 1..window_size (1-based POS)
        (pos - 1) / self.window_size
    }

    fn flush_chrom(&mut self) {
        if self.cur_chrom.is_empty() {
            return;
        }
        let Some(max_bin) = self.max_bin else {
            return;
        };
        for k in 0..=max_bin {
            let (pi_sum, n_var) = self.bins.get(&k).copied().unwrap_or((0.0, 0));
            let bin_start = k * self.window_size + 1;
            let bin_end = (k + 1) * self.window_size;
            let n_mono = self.window_size - n_var;
            let pi = pi_sum / self.window_size as f64;
            self.rows.push(WindowPiRow {
                chrom: self.cur_chrom.clone(),
                bin_start,
                bin_end,
                n_variants: n_var,
                n_monomorphic: n_mono,
                pi,
            });
        }
        self.bins.clear();
        self.max_bin = None;
    }

    fn push_site(&mut self, chrom: &str, pos: u64, pi_site: f64) {
        if chrom != self.cur_chrom {
            self.flush_chrom();
            self.cur_chrom = chrom.to_string();
        }
        let k = self.bin_of(pos);
        let e = self.bins.entry(k).or_insert((0.0, 0));
        e.0 += pi_site;
        e.1 += 1;
        self.max_bin = Some(self.max_bin.map_or(k, |m: u64| m.max(k)));
    }

    fn finish(mut self) -> Vec<WindowPiRow> {
        self.flush_chrom();
        self.rows
    }
}

// ── Per-site π (Tajima's estimator) ─────────────────────────────────────────

/// Compute π for one biallelic site from a genotype string slice.
/// Returns `None` if no called genotypes exist.
///
/// `π = (2n / (2n-1)) * 2 * p * (1-p)` where p = ALT frequency, 2n = haplotype count.
pub fn site_pi(gt_fields: &[&str]) -> Option<f64> {
    let mut n_ref: u64 = 0;
    let mut n_alt: u64 = 0;
    for field in gt_fields {
        let gt = if let Some(c) = field.find(':') {
            &field[..c]
        } else {
            field
        };
        for a in gt.split(['/', '|']) {
            match a {
                "0" => n_ref += 1,
                "." => {}
                _ if a.parse::<u64>().is_ok_and(|v| v > 0) => n_alt += 1,
                _ => {}
            }
        }
    }
    let two_n = n_ref + n_alt;
    if two_n < 2 {
        return None;
    }
    let p = n_alt as f64 / two_n as f64;
    let q = 1.0 - p;
    // Tajima's estimator with finite-sample correction
    Some((two_n as f64 / (two_n - 1) as f64) * 2.0 * p * q)
}

// ── VCF scanner ─────────────────────────────────────────────────────────────

fn open_reader(path: &Path) -> Result<Box<dyn Read>> {
    let file = std::fs::File::open(path).map_err(|e| {
        RsomicsError::Io(std::io::Error::new(
            e.kind(),
            format!("cannot open {}: {e}", path.display()),
        ))
    })?;
    let is_gz = path
        .extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"));
    Ok(if is_gz {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    })
}

const FIRST_SAMPLE: usize = 9;
const COL_CHROM: usize = 0;
const COL_POS: usize = 1;
const COL_ALT: usize = 4;

pub fn compute_window_pi(path: &Path, window_size: u64) -> Result<Vec<WindowPiRow>> {
    let reader = open_reader(path)?;
    let mut lines = BufReader::new(reader).lines();
    let mut acc = WindowAcc::new(window_size);
    let mut found_chrom = false;

    for line in lines.by_ref() {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if line.starts_with('#') {
            found_chrom = true;
            continue;
        }
        if !found_chrom {
            return Err(RsomicsError::InvalidInput(
                "VCF missing #CHROM header line".into(),
            ));
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= FIRST_SAMPLE {
            continue;
        }
        // Only biallelic sites (no comma in ALT, non-missing ALT)
        let alt = cols[COL_ALT];
        if alt == "." || alt.contains(',') {
            continue;
        }
        let pos: u64 = match cols[COL_POS].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let chrom = cols[COL_CHROM];
        let gt_fields = &cols[FIRST_SAMPLE..];
        // Only count polymorphic sites in N_VARIANTS; monomorphic (pi=0) are
        // not counted by vcftools even though they contribute zero to window PI.
        if let Some(pi) = site_pi(gt_fields) {
            if pi > 0.0 {
                acc.push_site(chrom, pos, pi);
            }
        }
    }

    if !found_chrom {
        return Err(RsomicsError::InvalidInput(
            "VCF missing #CHROM header line".into(),
        ));
    }
    Ok(acc.finish())
}

pub fn header() -> &'static str {
    "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g_zero() {
        assert_eq!(format_g(0.0), "0");
    }

    #[test]
    fn format_g_simple() {
        assert_eq!(format_g(0.0012), "0.0012");
    }

    #[test]
    fn format_g_repeating() {
        // 1/3000 = 0.000333333...
        assert_eq!(format_g(1.0 / 3000.0), "0.000333333");
    }

    #[test]
    fn site_pi_monomorphic_ref() {
        // All REF: p=0, pi=0
        let fields = ["0/0", "0/0", "0/0"];
        let pi = site_pi(&fields).unwrap();
        assert_eq!(pi, 0.0);
    }

    #[test]
    fn site_pi_all_het() {
        // S1=0/0, S2=0/1, S3=1/1 → REF=3, ALT=3, 2n=6, p=0.5
        // pi = (6/5) * 2 * 0.5 * 0.5 = 0.6
        let fields = ["0/0", "0/1", "1/1"];
        let pi = site_pi(&fields).unwrap();
        let expected = (6.0 / 5.0) * 2.0 * 0.5 * 0.5;
        assert!((pi - expected).abs() < 1e-15);
    }

    #[test]
    fn site_pi_rare_allele() {
        // S1=0/0, S2=0/0, S3=0/1 → REF=5, ALT=1, 2n=6, p=1/6
        // pi = (6/5) * 2 * (5/6) * (1/6) = (6/5) * 10/36 = 1/3
        let fields = ["0/0", "0/0", "0/1"];
        let pi = site_pi(&fields).unwrap();
        let expected = (6.0 / 5.0) * 2.0 * (5.0 / 6.0) * (1.0 / 6.0);
        assert!((pi - expected).abs() < 1e-15);
    }

    #[test]
    fn site_pi_missing_skipped() {
        let fields = ["./.", "0/1"];
        let pi = site_pi(&fields).unwrap();
        // Only S2 counted: REF=1, ALT=1, 2n=2
        // pi = (2/1) * 2 * 0.5 * 0.5 = 1.0 — NOT the n/(n-1) form with n=1 haploid pair
        let expected = (2.0 / 1.0) * 2.0 * 0.5 * 0.5;
        assert!((pi - expected).abs() < 1e-15);
    }

    #[test]
    fn site_pi_all_missing_returns_none() {
        let fields = ["./.", "./."];
        assert!(site_pi(&fields).is_none());
    }
}
