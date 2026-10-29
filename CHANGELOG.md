# Changelog

All notable changes to STRdust are documented here.

## [1.3.0]

### Added

- **The VCF ID column carries the catalog's locus name.** STRchive BED files, including the
  catalog `--pathogenic` downloads (name in column 5), and TRGT catalogs (`ID=` in column 4)
  now give e.g. `HD_HTT` or `CANVAS_RFC1` instead of `.`, so calls no longer have to be mapped
  back to genes by coordinate. Other BED files are unchanged. Taken from #2 by fellen31.
- **`MOT`: the repeat units of each allele, detected from its sequence.** A new per-allele
  FORMAT field listing up to three units joined by `/`, strongest first: `AAAGA` at RFC1
  (the benign AAAAG unit), `CAG/CAA` at TBP (its CAA interruptions), and a motif switch such
  as `AAGGG/AAAAG` or `CCTT/CT`, as expansions at GOLGA8A can carry. Found as the period
  (1-100) with the best tandem stretch, so it works inside padded catalog intervals, and
  written in the reference's rotation so alleles at a locus compare directly. On 10,000
  adotto intervals on chr2, 24% of called loci have no tandem repeat (`.`) and 32% of the
  rest list more than one unit. It costs no measurable time, and every other field is
  unchanged. The FORMAT column gains `MOT` before `PS`, so parsers that index the sample
  column by position rather than by the FORMAT keys need updating.

## [1.2.0]

### Added

- **Expansions too long for the reads are flagged, not called reference (#34).** When no
  read spans an expansion, its reads end in a large clip at the repeat and used to be
  discarded: phased input then no-called the whole locus, and `--unphased` reported a
  confident homozygous reference (0/0) from the other haplotype alone. Reads clipped by at
  least 100 bp at the repeat are now counted and reported as `UNSPANNED_EXPANSION=<n>` once
  there are at least `--support` of them. The haplotype they belong to is reported
  unresolved and the other one is still called: `0|.` for phased input, `0/.` unphased. On
  a synthetic ladder at the HTT CAG repeat (+100 bp to +8 kb, 10 and 20 kb reads), every
  expansion with spanning reads is sized within 0.1%, and the one without (+8 kb, 10 kb
  reads) is now `0|.` / `0/.` with the flag instead of a no-call or a false 0/0. With
  `--unphased`, only clips more than 200 bp longer than the called allele count, since
  reads not spanning the called allele itself clip too. On 5,072 HG002 long expansions:
  31 phased no-calls become half-calls, each called allele matching the truth; 16 unphased
  calls become half-calls, 13 of which had missed a second allele (up to 31 kb, e.g. 0/0
  where the truth is 0/22062). Loci without the flag are byte-identical to before.

### Changed

- **`--unphased ward` is about 10x faster, with identical output.** Ward's pairwise read
  distances now use Myers' bit-parallel edit distance instead of the cell-by-cell dynamic
  programme: the same exact distance, computed 64 cells at a time. Unphased runs of 10,000
  loci went from 214 to 17 CPU seconds, and from 3,177 to 379 on 5,072 long expansions.

- **Length outliers are removed from tighter read clusters.** Before the consensus, reads
  more than two standard deviations from the cluster's mean length are dropped; this used to
  happen only when the standard deviation was at least 5 bp and now applies from 2 bp. On
  10,000 held-out HG002 loci: exact allele length +0.3 points, within 1 bp +0.5, with no
  change in no-calls and no change on 5,072 long expansions. Testing it also settled #36: the
  filter does not bias expansions short. Errors at long loci lean long, and switching the
  filter off makes them worse.

## [1.1.3]

### Fixed

- **`--threads 1` was not respected with `--unphased`.** The single-threaded path never
  configured rayon, so the parallel distance matrix used to split reads into haplotypes ran
  on rayon's default pool of one thread per core: on a test set `--unphased -t 1` used 4.3
  cores, and on a large server it could use them all. The thread pool is now built for every
  thread count. Phased input was not affected, in either mode.

## [1.1.2]

### Fixed

- **Builds with newer clang.** `hts-sys` 2.2.1 moves to bindgen 0.72, which generates the
  htslib bindings correctly against the libclang on ubuntu 26.04 and in current conda
  environments, so the CI pin to ubuntu-24.04 is dropped. The code is unchanged from 1.1.1.

## [1.1.1]

### Fixed

- **Release binaries.** v1.1.0 was published without them: GitHub moved `ubuntu-latest` to
  26.04, where the bindgen in `hts-sys` cannot generate the htslib bindings. The Linux
  workflows are pinned to ubuntu-24.04. The code is unchanged from 1.1.0.
- **Dependency advisories** RUSTSEC-2026-0338 to -0342, -0346 to -0349 and -0383, in
  `rustls`, `quinn-proto` and `webpki-root-certs`. All are reached only through `reqwest`,
  which STRdust uses to download the `--pathogenic` catalog.

## [1.1.0]

### Changed

- **The allele consensus is built with [poasta](https://github.com/broadinstitute/poasta)
  instead of rust-bio's POA.** On HG002 truth-set loci never used for tuning (two samples of
  10,000), exact allele-length concordance rose by 1.2-1.4 points in `fast` mode and 2.8-3.3
  in `sensitive`, within 1 bp by 0.4-0.7, and genotype concordance by 0.6-1.5; the gains in
  exact and within 1 bp concordance are significant in paired tests on both samples. A
  further 50,000 loci in `fast` mode confirm it: exact +1.3, within 1 bp +0.3 and genotype
  +0.9 points, improving or level in every size stratum. On 4,569 loci with an expansion
  over 200 bp, exact concordance rose by 1.0 (`fast`) and 1.9 (`sensitive`) points, with
  within 1 bp unchanged. `fast` mode is 3-4x quicker on random loci and 10x on long
  expansions.
- **`SC` is now a cost**: the cost of aligning the consensus back to its read graph, where 0
  is a perfect fit and lower is better. It used to be rust-bio's alignment score, roughly
  three times the allele length.
- **Allele lengths change** at about 10% of called loci in `fast` mode and 13% in
  `sensitive`, by 1-2 bases at 85% of them, so thresholds on `RB` and `FRB` should be
  re-checked (`MRL` does not depend on the consensus). Equally supported paths through a
  homopolymer are now resolved by the cluster's median read length (hidden `--poa-ties`):
  following the longer path, poasta's own choice, made alleles too long, and following the
  better-covered one made them too short. Near long expansions, where reads run short of the
  true allele, the median is not consulted.

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
