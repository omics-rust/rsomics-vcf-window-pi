# rsomics-vcf-window-pi

Nucleotide diversity (π) in fixed-width windows from VCF — byte-identical to `vcftools --window-pi`.

```
rsomics-vcf-window-pi [OPTIONS] <VCF>
```

## Usage

```
Options:
  --window-size <BP>  Window size in bp  [default: 10000]
  --json              Emit JSON envelope
```

Output columns: `CHROM BIN_START BIN_END N_VARIANTS N_MONOMORPHIC PI`

Only biallelic, polymorphic sites contribute to π. Monomorphic sites (all REF or
all ALT) count toward N_MONOMORPHIC but not N_VARIANTS, matching vcftools behavior.

## Algorithm

For each polymorphic biallelic site, the per-site π uses Tajima's estimator with
finite-sample correction:

```
π_site = (2n / (2n−1)) × 2p(1−p)
```

where 2n is the total haplotype count (2 × samples with called genotypes), and p is
the ALT allele frequency. Window π = Σ(π_site) / window_size.

## Performance

Measured on mini_m2 (aarch64-apple-darwin), 100k-variant VCF (5 samples, 3 chromosomes),
`--window-size 10000`, hyperfine --warmup 3 --runs 20:

| | wall time (mean ± σ) |
|---|---|
| rsomics-vcf-window-pi 0.1.0 | 42.7 ms ± 0.3 ms |
| vcftools 0.1.17 | 135.1 ms ± 5.3 ms |

**3.16× faster.** Output is byte-identical on the same input.

## Install

```
cargo install rsomics-vcf-window-pi
```

## Origin

This crate is an independent Rust reimplementation of `vcftools --window-pi` based on:

- The vcftools 0.1.17 man page and black-box behavior observation
- Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by
  DNA polymorphism. _Genetics_ 123(3), 585–595 (for the π estimator with finite-sample correction)

License: MIT OR Apache-2.0.
Upstream credit: [vcftools](https://vcftools.github.io/) (LGPLv3).
