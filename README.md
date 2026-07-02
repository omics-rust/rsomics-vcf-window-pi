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

Every polymorphic site — biallelic or multiallelic — contributes to π. Monomorphic
sites (all one allele) are absorbed into the invariant baseline and count toward
N_MONOMORPHIC, not N_VARIANTS. Windows containing no variant are not emitted. Only
fully-diploid sites are used; a genotype of ploidy > 2 is a fatal error, matching
vcftools ("Polyploidy found, and not supported").

## Algorithm

vcftools computes a window's π as the average pairwise difference across every base
in the window. A variant site with non-missing allele counts `c_i` contributes
`m = Σ_{i<j} c_i·c_j` differing chromosome pairs out of `C(g,2)` comparisons, where
`g = Σ c_i` is the non-missing chromosome count. Every invariant position (including
monomorphic sites and the base positions with no record) is assumed fully called,
contributing `C(T,2)` comparisons for `T = 2 × n_samples` and zero differences.

```
π = Σ_variant m / [ Σ_variant C(g,2) + (window_size − N_VARIANTS) · C(T,2) ]
```

Missing genotypes therefore stay in the denominator through `T`, and the PI column is
printed with C `printf %g` (6 significant figures, switching to scientific notation
below 1e−4).

## Boundaries

Byte-identity with vcftools 0.1.17 holds for every well-formed, coordinate-sorted,
equal-width VCF — including sites-only input (which, like vcftools, aborts with
"Require Genotypes in VCF file..."), `POS=0` telomere records (dropped, matching
vcftools' negative window index), and all valid GT encodings. On genuinely
pathological input vcftools is not a reliable oracle, so this crate stays
deterministic rather than reproducing its bugs:

- **Unsorted / interleaved chromosomes.** vcftools re-accumulates cumulative windows
  across the whole file when a chromosome reappears after another, an artifact of its
  single running accumulator. This crate flushes per contiguous chromosome block.
- **Out-of-range allele index in a multi-window file.** An allele index beyond
  `1 + #ALT` drives vcftools into a stale-buffer read that intermittently `SIGBUS`es
  or corrupts memory (nondeterministic across runs). This crate charges the extra
  chromosome to the pairwise total deterministically and never crashes.

## Performance

Measured on mini_m2 (aarch64-apple-darwin), 100k-variant VCF (5 samples, 3 chromosomes),
`--window-size 10000`, hyperfine --warmup 3 --runs 20:

| | wall time (mean ± σ) |
|---|---|
| rsomics-vcf-window-pi 0.1.1 | 42.7 ms ± 0.3 ms |
| vcftools 0.1.17 | 135.1 ms ± 5.3 ms |

**3.16× faster.** Output is byte-identical on the same input.

## Install

```
cargo install rsomics-vcf-window-pi
```

## Origin

This crate is an independent Rust reimplementation of `vcftools --window-pi` based on:

- The vcftools 0.1.17 man page and black-box behavior observation
- Nei, M. & Li, W.-H. (1979). Mathematical model for studying genetic variation in terms
  of restriction endonucleases. _PNAS_ 76(10), 5269–5273 (average pairwise nucleotide difference)

License: MIT OR Apache-2.0.
Upstream credit: [vcftools](https://vcftools.github.io/) (LGPLv3).
