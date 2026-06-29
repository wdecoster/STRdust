# Changelog

All notable changes to STRdust are documented here.

## [0.20.0]

### Changed (breaking)

- **Haploid genotypes are now reported as a single allele value.** For loci on
  chromosomes listed under `--haploid`, the `GT` field is now a single value
  (e.g. `1`, `0`, or `.` when missing) instead of the previous diploid
  representation (`1/1`, `0/0`, `./.`). This follows the VCF specification, which
  asks for one allele value at haploid loci (e.g. male non-PAR X, Y, mitochondria).
- Per-allele `FORMAT`/`INFO` fields (`RB`, `FRB`, `MRL`, `SUP`, `SC`, `STDEV`)
  likewise carry a **single value** at haploid loci instead of a duplicated pair.
  Their header `Number` changed from `2` to `.` to allow mixed ploidy in one file.
- STRdust now prints a warning to stderr when `--haploid` is used, noting the
  changed output so downstream tooling is not silently surprised.

### Notes

- `--haploid` remains whole-chromosome. Loci in pseudoautosomal regions (PAR1/PAR2),
  which are diploid even on sex chromosomes, are **not** special-cased. This has no
  practical effect on current pathogenic STR catalogs, none of whose loci fall in a PAR.
- **Downstream impact:** mainstream tools (bcftools, GATK, vcftools, plink) handle
  haploid/mixed-ploidy genotypes. Custom parsers that assume a diploid `GT` (e.g.
  always splitting on `/` or `|` into two alleles) may need updating.
