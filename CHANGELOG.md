# Changelog

All notable changes to STRdust are documented here.

## [1.0.0]

Benchmarked for the first time against a truth set (GIAB HG002 tandem repeats over the
adotto v1.2.1 catalog, ~30x ONT), rather than against STRdust's own other mode. Exact
allele-length concordance on loci never used for tuning rose from **53.0% to 85.0%**. Most
of that is one upstream bug; the rest is defaults that had never been measured.

### Changed (breaking)

- **`--mode` now defaults to `fast`** (was `sensitive`). `fast` is both quicker and more
  accurate, including at long expansions — the case it was expected to lose. On 5,072 loci
  carrying an expansion over 200 bp it scores 71.9% exact against sensitive's 69.1%, with
  fewer no-calls; overall it is ~5x cheaper. `sensitive` re-aligns each read to a
  repeat-compressed reference and rebuilds the allele by concatenating insertions near the
  junction, which folds in stray flanking sequence and cannot subtract deletions; `fast`
  derives the allele from a reference span and gets both right for free.
- **Allele lengths change.** The POA consensus was running long on a quarter of alleles
  (see *Fixed* below), so `RB`, `FRB` and `MRL` shift for many loci — almost always shorter,
  and closer to the truth. Any threshold on those fields should be re-checked.
- **`--unphased` now takes the clustering strategy and `--phasing` is gone.**
  `--unphased --phasing dbscan` becomes `--unphased dbscan`; values are `ward`, `dbscan` and
  `both`. The value is required: with two positional arguments an optional-value flag would
  be ambiguous. On ordinary samples `ward` is the right choice — see #32.
- **`--minlen` defaults to 3** (was 1). The drop to 1 in 2025-11 was an unintended
  regression worth about 11 points of concordance; insertions shorter than a few bases at the
  repeat junction are read noise far more often than sequence.
- **Removed `--poa-gap-open`, `--poa-match`, `--poa-mismatch`.** Tuning the POA scoring
  gained up to 11 points on its own, but −0.1 once the consensus endpoint is fixed: both were
  addressing the same defect, and one addresses it at the source.

### Added

- **`--priority {expanded,sensitive,balanced,precise}`**, default `balanced`, replacing a
  reference/variant rule nobody chose. The old rule stacked a 0.9 length-ratio prefilter, a
  `len/20` floor division and a strictly-less-than comparison, landing at 52.1% recall /
  83.0% precision — on no efficient frontier. Measured alternatives:

  | value | recall | precision | F1 |
  |---|---|---|---|
  | `sensitive` | 99.0% | 54.2% | 70.1 |
  | `balanced` | 85.5% | 58.7% | 69.6 |
  | `precise` | 39.9% | 91.8% | 55.7 |

  `sensitive`, `balanced` and `precise` change **only the emitted `GT`** — `RB`, `FRB` and
  `MRL` are byte-identical under all three. `expanded` differs in kind: it reports a locus as
  reference *without genotyping it* when every read looks near-reference, which is ~14%
  faster and improves overall concordance, at the cost of resolution below ~3 bases
  (`1-10bp` exact falls about 6 points). Long expansions are unaffected. At skipped loci
  `RB` and `MRL` are `0` — a real finding rather than a placeholder, since the check
  established every read matches the reference length — while `SUP` and `SC` are `.`. Note
  this often *adds* length information: many such loci would otherwise have been no-calls
  with `RB` of `.`.

- **Reference calls now use the same `FORMAT` as every other record.** Loci called
  homozygous reference by the fast CIGAR check previously emitted `GT:SUP`, dropping `RB`,
  `FRB` and `MRL` even though all three were known. A `FORMAT` that varies between records
  in one file breaks any reader that indexes the sample column positionally, and the omission
  discarded a genuine measurement. They now report `RB=0`, `FRB=len(REF)`, `MRL=0`.

### Fixed

- **The POA consensus ran to the deepest sink rather than the best-supported path.**
  `bio::alignment::poa::Aligner::consensus()` picks its endpoint with `max_by_key` over a
  cumulative score to which every edge contributes at least 1, so the score rises strictly
  along every edge and the argmax is *always* a sink. A single read extending one base past
  the others created a deeper sink and won, with no vote taken. The excess was 29:1 one-sided
  and got **worse with more reads** (14.9% at 0-4 reads, 32.6% at 15-19) — a consensus cannot
  do that. Corrected by a traversal that refuses to run out along ends carried by a minority
  of reads. Worth **+26.0 points** of exact concordance on held-out loci. Reported upstream.
- **`--unphased` panicked with more than one thread** (`RefCell already borrowed`). The
  thread-local BAM reader borrow was held across the caller's closure, and rayon could
  schedule the nested clustering task onto the same thread. Load-dependent; phased input was
  never affected.
- **Read downsampling was unseeded**, so two identical invocations could report different
  alleles at loci above the read cap. Now seeded, matching the consensus downsampler.
- **`Batch::new` took the batch end from the last repeat** rather than the furthest, so a
  repeat nested inside an earlier, longer one could truncate the fetch region. No such pairs
  exist in the adotto catalog; it would affect BEDs with nested or overlapping intervals.
- **`VCFRecord::single_read` still used `end - start`**, the last survivor of the 0.21.0
  off-by-one, making single-read calls one base too long and underflowing on contractions.

### Notes

- All measurements come from one sample, one chemistry and one catalog. The orderings should
  generalise; the absolute numbers should not be quoted as universal.
- The investigation behind this release was kept in `KNOWN_ISSUES.md`, which is retired at
  1.0.0. Everything still outstanding became an issue (#32, #34-#38); the document itself
  remains in history at
  [`3fba652a171d01b89cb723cfa02f3cf6494b6693`](https://github.com/wdecoster/STRdust/blob/3fba652a171d01b89cb723cfa02f3cf6494b6693/KNOWN_ISSUES.md),
  including the conclusions that were overturned along the way.

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
