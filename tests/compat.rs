//! Value-exact compatibility with vcftools 0.1.17 `--window-pi`.
//!
//! Every EXPECTED block below is frozen output captured once from
//! vcftools 0.1.17 (`vcftools --vcf IN --window-pi W`). vcftools is never
//! invoked at test time — the goldens are hardcoded string constants and the
//! computation runs entirely in-process over an in-memory reader.

use std::io::Cursor;

use rsomics_vcf_window_pi::{compute_window_pi_from_reader, header};

fn run(vcf: &str, window_size: u64) -> String {
    let rows = compute_window_pi_from_reader(Cursor::new(vcf), window_size).unwrap();
    let mut out = header().to_string();
    for row in &rows {
        out.push_str(&row.to_text());
    }
    out
}

/// `(name, vcf, window, expected)` — expected is verbatim vcftools 0.1.17 output.
struct Case {
    name: &'static str,
    vcf: &'static str,
    window: u64,
    expected: &'static str,
}

const CASES: &[Case] = &[
    // Multi-chromosome, no missing data.
    Case {
        name: "basic_multichrom",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n\
chr1\t1100\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\n\
chr2\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.0012\n\
chr1\t1001\t2000\t1\t999\t0.000333333\n\
chr2\t1\t1000\t1\t999\t0.0006\n",
    },
    // Missing genotypes: denominator weighted by 2*n_samples (missing included).
    Case {
        name: "missing_genotypes",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t50\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t./.\t./.\t0/0\n\
chr1\t400\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/1\t./.\t1/1\n\
chr1\t5000\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\t./.\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.000393349\n\
chr1\t4001\t5000\t1\t999\t0.000178654\n",
    },
    // Half-calls (0/.) and phased genotypes (0|1).
    Case {
        name: "halfcall_phased",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t10\t.\tA\tT\t60\tPASS\t.\tGT\t0|1\t0/.\t1|1\n\
chr1\t20\t.\tA\tT\t60\tPASS\t.\tGT\t0|0\t0|1\t./1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.000800534\n",
    },
    // Multiallelic sites are INCLUDED (not skipped).
    Case {
        name: "multiallelic",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT,G\t60\tPASS\t.\tGT\t0/1\t1/2\t2/2\n\
chr1\t200\t.\tA\tC,G,T\t60\tPASS\t.\tGT\t0/3\t1/2\t0/0\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.00153333\n",
    },
    // Gap windows with zero variants are NOT emitted; PI uses %g scientific < 1e-4.
    Case {
        name: "gap_windows_scientific",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\n\
chr1\t9500\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/0\n\
chr1\t25000\t.\tA\tT\t60\tPASS\t.\tGT\t1/1\t0/1\n",
        window: 10000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t10000\t2\t9998\t0.0001\n\
chr1\t20001\t30000\t1\t9999\t5e-05\n",
    },
    // Monomorphic sites (all-REF, all-ALT) are absorbed into the invariant baseline.
    Case {
        name: "monomorphic_absorbed",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\n\
chr1\t200\t.\tA\tT\t60\tPASS\t.\tGT\t1/1\t1/1\t1/1\t1/1\n\
chr1\t300\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.00025\n",
    },
    // All-missing site among variants: contributes nothing, not a variant.
    Case {
        name: "all_missing_site",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t./.\t./.\t./.\n\
chr1\t200\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/1\t1/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.000533333\n",
    },
    // GT not the first FORMAT field.
    Case {
        name: "gt_not_first",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tDP:GT\t10:0/1\t9:1/1\t8:0/0\n\
chr1\t200\t.\tA\tT\t60\tPASS\t.\tDP:GT\t5:0/0\t7:0/1\t6:0/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t2\t998\t0.00113333\n",
    },
    // Two chromosomes, gaps, missing, small windows, scientific PI.
    Case {
        name: "two_chrom_gaps",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chrX\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/0\t./.\n\
chrX\t12000\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\t0/1\n\
chrY\t300\t.\tA\tT\t60\tPASS\t.\tGT\t1/1\t0/1\t0/0\n",
        window: 5000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chrX\t1\t5000\t1\t4999\t4.00048e-05\n\
chrX\t10001\t15000\t1\t4999\t0.000106667\n\
chrY\t1\t5000\t1\t4999\t0.00012\n",
    },
    // Out-of-range allele indices (idx ≥ 1 + #ALT): each out-of-range copy
    // lifts the total chromosome count but joins no per-allele count.
    Case {
        name: "out_of_range_single_sample",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t5\t.\tA\tG\t.\t.\t.\tGT\t1/2\n\
chr1\t8\t.\tA\tG,T\t.\t.\t.\tGT\t2/3\n\
chr1\t1500\t.\tA\tG\t.\t.\t.\tGT\t0/1\n",
        window: 10,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t10\t2\t8\t0.1\n\
chr1\t1491\t1500\t1\t9\t0.1\n",
    },
    // Out-of-range mixed with in-range and missing genotypes; ALT='.' → 1 allele.
    Case {
        name: "out_of_range_mixed_multisample",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/2\t2/2\n\
chr1\t300\t.\tA\tG,T\t.\t.\t.\tGT\t0/3\t1/2\t./.\n\
chr1\t900\t.\tA\t.\t.\t.\t.\tGT\t0/1\t0/0\t0/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t3\t997\t0.0010006\n",
    },
    // A wholly out-of-range site has zero differing pairs → absorbed as
    // monomorphic, not counted as a variant.
    Case {
        name: "out_of_range_absorbed_monomorphic",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t2/2\t3/3\n\
chr1\t400\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.0005\n",
    },
    // POS=0 (telomere convention): vcftools bins a site by (POS−1)/W, so POS=0
    // lands on a negative window index it discards — the record produces no row.
    Case {
        name: "pos_zero_dropped",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t0\t.\tA\tT\t60\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/0\n\
chr1\t1100\t.\tC\tT\t60\tPASS\t.\tGT\t0/0\t0/0\t0/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.0006\n\
chr1\t1001\t2000\t1\t999\t0.000333333\n",
    },
    // POS=0 as the sole record on a chromosome: that chromosome emits nothing.
    Case {
        name: "pos_zero_first_on_chrom",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t0\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t1/1\n\
chr2\t500\t.\tA\tG\t60\tPASS\t.\tGT\t1/1\t0/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr2\t1\t1000\t1\t999\t0.0005\n",
    },
    // Single sample.
    Case {
        name: "single_sample",
        vcf: "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n\
chr1\t700\t.\tA\tT\t60\tPASS\t.\tGT\t1/1\n\
chr1\t1500\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\n",
        window: 1000,
        expected: "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n\
chr1\t1\t1000\t1\t999\t0.001\n\
chr1\t1001\t2000\t1\t999\t0.001\n",
    },
];

#[test]
fn window_pi_matches_goldens() {
    for c in CASES {
        assert_eq!(
            run(c.vcf, c.window),
            c.expected,
            "case `{}` diverged",
            c.name
        );
    }
}

/// vcftools aborts with "Polyploidy found, and not supported"; we must fail too.
#[test]
fn polyploid_is_rejected() {
    let vcf = "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\t0/1/1\n";
    assert!(compute_window_pi_from_reader(Cursor::new(vcf), 1000).is_err());
}

/// A VCF with no genotype columns (sites-only) makes vcftools abort with
/// "Require Genotypes in VCF file..." and exit 1; we must fail the same way.
#[test]
fn sites_only_requires_genotypes() {
    let sites_only = "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\n";
    let err = compute_window_pi_from_reader(Cursor::new(sites_only), 1000).unwrap_err();
    assert!(err.to_string().contains("Require Genotypes"));

    // FORMAT column present but still no sample columns is the same abort.
    let format_no_samples = "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\n";
    assert!(compute_window_pi_from_reader(Cursor::new(format_no_samples), 1000).is_err());
}

/// Haploid genotypes make a site "not fully diploid" — vcftools skips it,
/// yielding no windows here.
#[test]
fn haploid_site_skipped() {
    let vcf = "##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0\t1\t1\n";
    assert_eq!(
        run(vcf, 1000),
        "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI\n"
    );
}
