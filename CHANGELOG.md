# Changelog

All notable changes to STRdust are documented here.

## [0.21.0]

### Changed (breaking)

- **`RB` and `MRL` are one base lower than before.** Both were computed against
  `end - start`, but the annotated repeat (1-based inclusive start, BED end) spans
  `end - start + 1` reference bases - the same length as the `REF` sequence in the
  record. Every `RB` and `MRL` value STRdust has produced was therefore one too high,
  and `RB` did not equal `FRB - len(REF)` as its header description implies. Both now
  derive from the `REF` sequence itself. `DBSCAN_RB` carried the same offset and is
  fixed with them. `FRB` (the raw consensus length) is unaffected.
  **Downstream impact:** any threshold on `RB` or `MRL`, and any comparison of values
  across STRdust versions, shifts by one base.

### Added

- `--mapq` (default 10) sets the minimum mapping quality of a read to be used, and is
  applied consistently in both the batched and non-batched paths. Previously only reads
  with mapping quality 0 were dropped, and only in the non-batched path, so the batched
  path - which is what runs for every locus with coverage - applied no filter at all.
  Pass `--mapq 0` to keep ambiguously mapped reads, which can matter in segmental
  duplications.

### Notes

- Reads that merely overlap a locus, rather than spanning it, are still used by the
  `sensitive` path: re-aligning them to the repeat-compressed reference is part of how
  it recovers alleles the original alignment placed badly. `--mode fast` cannot use them
  and drops them, falling back to `sensitive` when too few reads are left.

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
