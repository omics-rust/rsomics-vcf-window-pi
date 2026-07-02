//! Nucleotide diversity (π) in fixed-width windows from VCF, matching
//! vcftools 0.1.17 `--window-pi` byte-for-byte.
//!
//! Each chromosome is split into non-overlapping bins of `window_size` bp
//! (1-based, bin `k` covers `k*W+1 ..= (k+1)*W`). For every bin that contains
//! at least one polymorphic site vcftools reports
//!
//! ```text
//! π = Σ_variant d_site / [ Σ_variant g·(g−1) + (W − N_variant) · T·(T−1) ]
//! ```
//!
//! where at a site `g` is the number of non-missing chromosomes and
//! `d = Σ_a c_a·(g − c_a)` is twice the number of differing chromosome pairs,
//! summed over the per-allele counts `c_a`. `T = 2·n_samples` is the full
//! chromosome count assumed for every invariant position, and a site counts as
//! a variant iff `d > 0`. Working in the doubled numerator/denominator keeps the
//! arithmetic integral: an out-of-range allele index (`≥ 1 + #ALT`) raises `g`
//! but joins no `c_a`, so it can leave `d` odd — vcftools charges such a
//! chromosome to the total pairwise count without ever pairing it to an allele.
//!
//! Only fully-diploid sites contribute: a site with any haploid genotype is
//! skipped, and any genotype of ploidy > 2 is a fatal error (vcftools:
//! "Polyploidy found, and not supported").

use std::collections::BTreeMap;
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

/// C `printf` `%g` with the default precision of 6 significant figures.
pub fn format_g(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    if x == 0.0 {
        return "0".to_string();
    }

    const SIG: i32 = 6;
    // Decompose via scientific notation rounded to SIG figures so that any
    // rounding that bumps the exponent (9.999995e-5 → 1e-4) is already applied.
    let sci = format!("{:.*e}", (SIG - 1) as usize, x);
    let (mantissa, exp_str) = sci.split_once('e').expect("scientific format has 'e'");
    let exp: i32 = exp_str.parse().expect("exponent is an integer");

    if !(-4..SIG).contains(&exp) {
        let mantissa = mantissa.trim_end_matches('0').trim_end_matches('.');
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{mantissa}e{sign}{:02}", exp.abs())
    } else {
        let decimals = (SIG - 1 - exp).max(0) as usize;
        let fixed = format!("{x:.decimals$}");
        if fixed.contains('.') {
            fixed
                .trim_end_matches('0')
                .trim_end_matches('.')
                .to_string()
        } else {
            fixed
        }
    }
}

// ── Per-site diploid statistics ──────────────────────────────────────────────

/// Doubled pairwise-difference and comparison counts for one site.
struct SiteStat {
    /// `Σ_a c_a·(g − c_a)` — twice the number of differing chromosome pairs.
    diffs: u128,
    /// `g·(g − 1)` — twice `C(g, 2)`, `g` the non-missing chromosome count.
    comparisons: u128,
}

/// Parse allele statistics for one site.
///
/// `n_alleles = 1 + #ALT`. A genotype allele index at or beyond `n_alleles` is
/// out of range: vcftools counts that chromosome into the total `g` but into no
/// per-allele count, so it lifts the pairwise total without pairing to a base.
///
/// Returns `Ok(None)` when the site must be skipped (not fully diploid), and a
/// fatal error when any genotype has ploidy > 2 (vcftools rejects polyploidy).
fn site_stat(
    gt_fields: &[&str],
    gt_index: usize,
    n_alleles: u64,
    chrom: &str,
    pos: u64,
) -> Result<Option<SiteStat>> {
    let mut counts: BTreeMap<u64, u128> = BTreeMap::new();
    let mut g: u128 = 0;
    let mut fully_diploid = true;

    for field in gt_fields {
        let gt = field.split(':').nth(gt_index).unwrap_or("");
        let mut ploidy = 0u32;
        for allele in gt.split(['/', '|']) {
            ploidy += 1;
            if allele == "." {
                continue;
            }
            let idx: u64 = allele.parse().map_err(|_| {
                RsomicsError::InvalidInput(format!(
                    "malformed genotype allele '{allele}' at {chrom}:{pos}"
                ))
            })?;
            g += 1;
            if idx < n_alleles {
                *counts.entry(idx).or_insert(0) += 1;
            }
        }
        if ploidy > 2 {
            return Err(RsomicsError::InvalidInput(format!(
                "Polyploidy found, and not supported: {chrom}:{pos}"
            )));
        }
        if ploidy != 2 {
            fully_diploid = false;
        }
    }

    if !fully_diploid {
        return Ok(None);
    }

    let diffs = counts.values().map(|&c| c * (g - c)).sum();
    let comparisons = g * g.saturating_sub(1);
    Ok(Some(SiteStat { diffs, comparisons }))
}

// ── Window accumulator ───────────────────────────────────────────────────────

#[derive(Default)]
struct BinAcc {
    diffs: u128,
    comparisons: u128,
    n_variants: u64,
}

struct WindowAcc {
    window_size: u64,
    /// `T·(T − 1)` for a fully-called invariant position, T = 2·n_samples.
    invariant_comparisons: u128,
    cur_chrom: String,
    bins: BTreeMap<u64, BinAcc>,
    rows: Vec<WindowPiRow>,
}

impl WindowAcc {
    fn new(window_size: u64, n_samples: u64) -> Self {
        let t = 2 * n_samples;
        Self {
            window_size,
            invariant_comparisons: t as u128 * (t.saturating_sub(1)) as u128,
            cur_chrom: String::new(),
            bins: BTreeMap::new(),
            rows: Vec::new(),
        }
    }

    fn flush_chrom(&mut self) {
        if self.cur_chrom.is_empty() {
            return;
        }
        let w = self.window_size;
        for (&k, bin) in &self.bins {
            let baseline = (w - bin.n_variants) as u128 * self.invariant_comparisons;
            let denom = bin.comparisons + baseline;
            let pi = bin.diffs as f64 / denom as f64;
            self.rows.push(WindowPiRow {
                chrom: self.cur_chrom.clone(),
                bin_start: k * w + 1,
                bin_end: (k + 1) * w,
                n_variants: bin.n_variants,
                n_monomorphic: w - bin.n_variants,
                pi,
            });
        }
        self.bins.clear();
    }

    fn push_variant(&mut self, chrom: &str, pos: u64, stat: &SiteStat) {
        if chrom != self.cur_chrom {
            self.flush_chrom();
            self.cur_chrom = chrom.to_string();
        }
        let bin = self.bins.entry((pos - 1) / self.window_size).or_default();
        bin.diffs += stat.diffs;
        bin.comparisons += stat.comparisons;
        bin.n_variants += 1;
    }

    fn finish(mut self) -> Vec<WindowPiRow> {
        self.flush_chrom();
        self.rows
    }
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
const COL_FORMAT: usize = 8;

pub fn compute_window_pi(path: &Path, window_size: u64) -> Result<Vec<WindowPiRow>> {
    compute_window_pi_from_reader(open_reader(path)?, window_size)
}

pub fn compute_window_pi_from_reader<R: Read>(
    reader: R,
    window_size: u64,
) -> Result<Vec<WindowPiRow>> {
    let lines = BufReader::new(reader).lines();

    let mut acc: Option<WindowAcc> = None;

    for line in lines {
        let line = line?;
        let line = line.trim_end_matches('\r');
        if line.starts_with("##") {
            continue;
        }
        if let Some(rest) = line.strip_prefix('#') {
            let n_samples = rest.split('\t').count().saturating_sub(FIRST_SAMPLE) as u64;
            if n_samples == 0 {
                return Err(RsomicsError::InvalidInput(
                    "Require Genotypes in VCF file in order to output Nucleotide Diversity \
                     Statistics."
                        .into(),
                ));
            }
            acc = Some(WindowAcc::new(window_size, n_samples));
            continue;
        }
        let acc = acc
            .as_mut()
            .ok_or_else(|| RsomicsError::InvalidInput("VCF missing #CHROM header line".into()))?;

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= FIRST_SAMPLE {
            continue;
        }
        let pos: u64 = match cols[COL_POS].parse() {
            Ok(v) => v,
            Err(_) => continue,
        };
        let chrom = cols[COL_CHROM];
        let alt = cols[COL_ALT];
        let n_alleles = if alt == "." {
            1
        } else {
            1 + alt.split(',').count() as u64
        };
        let gt_index = cols[COL_FORMAT]
            .split(':')
            .position(|f| f == "GT")
            .unwrap_or(0);

        if let Some(stat) = site_stat(&cols[FIRST_SAMPLE..], gt_index, n_alleles, chrom, pos)? {
            // vcftools bins a site by (POS−1)/W; POS=0 gives a negative index it discards.
            if pos != 0 && stat.diffs > 0 {
                acc.push_variant(chrom, pos, &stat);
            }
        }
    }

    acc.map(WindowAcc::finish)
        .ok_or_else(|| RsomicsError::InvalidInput("VCF missing #CHROM header line".into()))
}

pub fn header() -> &'static str {
    "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n"
}

#[cfg(test)]
mod tests {
    use super::*;

    fn g(x: f64) -> String {
        format_g(x)
    }

    #[test]
    fn format_g_matches_printf() {
        assert_eq!(g(0.0), "0");
        assert_eq!(g(0.0012), "0.0012");
        assert_eq!(g(1.0 / 3000.0), "0.000333333");
        assert_eq!(g(0.0006), "0.0006");
        assert_eq!(g(0.00025), "0.00025");
        assert_eq!(g(0.000166806), "0.000166806");
        assert_eq!(g(6.67289e-05), "6.67289e-05");
        assert_eq!(g(3.57488e-05), "3.57488e-05");
        assert_eq!(g(1.51665e-05), "1.51665e-05");
        assert_eq!(g(0.000833333), "0.000833333");
        assert_eq!(g(1.0), "1");
        assert_eq!(g(0.5), "0.5");
        assert_eq!(g(1.0 / 6.0), "0.166667");
        assert_eq!(g(0.0001), "0.0001");
        assert_eq!(g(9.999e-05), "9.999e-05");
        assert_eq!(g(1.23456789e-07), "1.23457e-07");
    }

    fn stat(fields: &[&str], n_alleles: u64) -> Option<SiteStat> {
        site_stat(fields, 0, n_alleles, "chr1", 1).unwrap()
    }

    #[test]
    fn diploid_het_counts() {
        let s = stat(&["0/0", "0/1", "1/1"], 2).unwrap();
        // REF=3, ALT=3 → d = 3·3 + 3·3 = 18, comparisons = 6·5 = 30.
        assert_eq!(s.diffs, 18);
        assert_eq!(s.comparisons, 30);
    }

    #[test]
    fn multiallelic_included() {
        // 0/1, 1/2 → counts 0:1, 1:2, 2:1; d = 1·3 + 2·2 + 1·3 = 10.
        let s = stat(&["0/1", "1/2"], 3).unwrap();
        assert_eq!(s.diffs, 10);
        assert_eq!(s.comparisons, 12);
    }

    #[test]
    fn half_call_counts_present_allele() {
        // 0/. , 1/1 → REF=1, ALT=2, d = 1·2 + 2·1 = 4, comparisons = 3·2 = 6.
        let s = stat(&["0/.", "1/1"], 2).unwrap();
        assert_eq!(s.diffs, 4);
        assert_eq!(s.comparisons, 6);
    }

    #[test]
    fn monomorphic_has_zero_diffs() {
        let s = stat(&["0/0", "0/0"], 2).unwrap();
        assert_eq!(s.diffs, 0);
    }

    #[test]
    fn out_of_range_lifts_total_only() {
        // ALT=G → n_alleles = 2. 1/2: index 2 out of range, counts only into g.
        // counts {1:1}, g = 2, d = 1·(2−1) = 1, comparisons = 2·1 = 2.
        let s = stat(&["1/2"], 2).unwrap();
        assert_eq!(s.diffs, 1);
        assert_eq!(s.comparisons, 2);
    }

    #[test]
    fn all_out_of_range_is_monomorphic() {
        // 2/2 with n_alleles = 2 → both copies out of range, no pairs differ.
        let s = stat(&["2/2"], 2).unwrap();
        assert_eq!(s.diffs, 0);
        assert_eq!(s.comparisons, 2);
    }

    #[test]
    fn out_of_range_with_missing() {
        // 1/2, ./. with ALT=G → g = 2 (one sample), counts {1:1}, d = 1.
        let s = stat(&["1/2", "./."], 2).unwrap();
        assert_eq!(s.diffs, 1);
        assert_eq!(s.comparisons, 2);
    }

    #[test]
    fn haploid_site_skipped() {
        assert!(stat(&["0", "1", "1"], 2).is_none());
        assert!(stat(&["0/1", "0"], 2).is_none());
        assert!(stat(&["0/1", "."], 2).is_none());
    }

    #[test]
    fn polyploid_is_error() {
        assert!(site_stat(&["0/1", "0/1/1"], 0, 2, "chr1", 100).is_err());
    }

    #[test]
    fn gt_index_respected() {
        // FORMAT = GT:DP with GT at index 0. REF=1, ALT=3 → d = 1·3 + 3·1 = 6.
        let s = stat(&["0/1:10", "1/1:9"], 2).unwrap();
        assert_eq!(s.diffs, 6);
    }
}
