# Known issues and open findings

Written 2026-09-08, from an investigation on 2026-08-28/29 that started as a review of
`--mode fast` and turned into a hunt for why STRdust scored poorly in the Genome Biology
Nanopore tandem-repeat benchmark (`s13059-026-04210-y_reference.pdf`, which benchmarked
v0.16.0).

All measurements below are on a 30x ONT sample, 2000 catalog loci on chr1:20–21.4 Mb,
unless stated otherwise. Current version: 0.21.0.

---

## 1. Summary

Nine defects were found. Five are fixed and merged, one is fixed on an open PR, and three
open issues describe work that changes behaviour and needs a decision.

| # | What | State |
|---|---|---|
| #19 | Multi-mapped reads counted more than once in `SUP` | fixed, merged (#26) |
| #20 | Batched path fed non-spanning, MAPQ-0 reads to the aligner | fixed, merged (#27, `--mapq`) |
| #21 | `--haploid` silently dropped HP-tagged reads in the batched path | fixed, merged (#26) |
| #22 | `RB`/`MRL`/`DBSCAN_RB` one base too high on every allele | fixed, merged (#27) |
| — | Three VCF conformance issues (bare `;`, `\|` vs `/`, undeclared version) | fixed, merged (#29) |
| #30 | QUICKREF called variant loci homozygous reference (coordinate off-by-one) | **PR #30 open — correct, but not a sufficient fix** |
| #31 | No tolerance for read noise; over-calls alternate alleles at short repeats | **open — the big one** |
| #24 | ALT can describe bases outside `[POS, END]` | open, documented, deliberate |
| #23 | Neither mode has ever been benchmarked against a truth set | open |

**Read §13 and §14 before acting on §2–§9.** Those early sections were written from
measurements taken against STRdust's own output; §13 and §14 are the first measurements
against a truth set, and they overturn three of the conclusions below. Where they disagree,
the truth-set measurement wins. In particular:

| earlier claim | what the benchmark showed |
|---|---|
| #31 (over-calls at short repeats) is the big one | Under-calling at long loci is the larger half; the threshold fails at **both** ends (§13.3) |
| QUICKREF is a crude heuristic worth little | The **tolerant** one fires at 32.9% with 96.8% precision — better than full genotyping (§14.2) |
| `--mode fast` trades accuracy for speed | `fast` is **more** accurate than `sensitive` on every combined metric, at half the CPU (§14.1) |
| #24's flank folding is cosmetic | Candidate root cause of a 29% insertion bias in the consensus (§13.5) |

The single most important finding is that **the POA consensus is the root cause** (§14.1):
reading lengths from the CIGAR gives an 8.2% over-call rate where the consensus gives 26.6%.
Everything else — the similarity threshold, QUICKREF's criterion — is downstream of that.

---

## 2. The benchmark result, and what explains it

The paper reports, for STRdust v0.16.0:

- overall concordance 72.56% (R10), rising to **94.93% when an off-by-motif deviation is permitted**;
- the only tool below 90% concordance at `0/0` genotypes — and `0/0` loci are **83.6% of all loci**;
- "STRdust tends to call more genotypes with an alternate allele";
- Mendelian consistency at `0/0` of 99.5% vs 99.8–99.99% for every other tool, "driven by a
  higher rate of motif-length-off calls at this class";
- below 40% consistency at mononucleotide repeats.

Four independent contributors were identified, in rough order of expected impact:

1. **No noise tolerance (#31)** — STRdust calls 83% of loci non-reference where the paper
   reports 83.6% invariant. Directly explains the `0/0` concordance, the alt-allele bias,
   the motif-length-off Mendelian errors and the homopolymer result.
2. **QUICKREF false negatives (#30)** — ~11% of loci silently emitted as `0/0` with
   `SUP=.,.`, so no filter could catch them. Hits recall hardest at *small* expansions,
   which are most of any catalog.
3. **`RB` off-by-one (#22)** — every allele length one base too high. Any benchmark scoring
   exact allele-length match marked STRdust wrong at loci it otherwise got right.
4. **Non-parsimonious records** — `bcftools norm` realigns 1557/2000 records. If the
   benchmark normalised before comparing, positions moved and records failed to match.

(2) and (3) are fixed; (4) is a deliberate design choice shared with TRGT; (1) is open.

---

## 3. Open: #31 — no tolerance for read noise

> **SUPERSEDED in its diagnosis — see §29.** This section's measurement stands: STRdust did
> call far too many loci non-reference. Its *explanation* was wrong. Two thirds of the
> over-calling was `rust-bio`'s POA consensus running to the deepest sink (§29.1), and most
> of the remainder was `is_similar_to_ref` rounding real 1-2 base alleles to reference
> (§29.6). False non-reference calls at truth-hom-ref loci fell from 25.2% to 7.9% by fixing
> those two, without touching the tolerance this section proposes. Kept because the
> measurement and the reasoning that led from it to the wrong culprit are the record of how
> the diagnosis was found.

### Measurement

With QUICKREF disabled so every locus is genotyped properly, STRdust calls
**1652 of 1999 loci (83%) non-reference**. The paper reports 83.6% of loci are invariant in
the R10 HPRC assemblies. Over half of the non-reference calls are noise-scale:

| largest \|RB\| per locus | count | cumulative |
|---|---|---|
| 0–1 bp | 262 | 15.9% |
| 2 bp | 263 | 31.8% |
| 3 bp | 341 | 52.4% |
| ≥ 4 bp | 786 | 100.0% |

Median REF length at these loci is 11 bp.

### Mechanism

`is_similar_to_ref` (`src/vcf.rs:458`) decides whether an allele counts as reference using an
edit-distance threshold of `reference.len() / 20` — 5%, floor division (`src/vcf.rs:478`).
The code comment states the consequence outright:

```rust
// For 12bp: threshold = 0, so only exact matches pass
```

At a median REF length of 11 bp, **most loci require the POA consensus to match the
reference exactly** to be called `0/0`. A POA consensus over ONT reads at a short STR
essentially never does, so a one-base stutter becomes `0/1`.

QUICKREF is the same intolerance in cruder form: it requires *every* overlapping read to
differ by exactly zero. Sweeping a tolerance there does not rescue it, because a few bp is
also the size of most genuine variation at these loci:

| tolerance | QUICKREF fires | of which wrong |
|---|---|---|
| 0 | 1/1999 | 0 |
| 3 | 26/1999 | 9 |
| 5 | 153/1999 | 99 |

Even at tolerance 10 it catches only 48% of the loci STRdust itself calls `0/0`.

### What can be done

The most useful hint in the paper is the 72.56% → 94.93% jump when an off-by-motif deviation
is permitted: **the natural unit of tolerance is one motif, not a percentage of interval
length.** STRdust currently does not know the motif at all — the catalog BED carries it in
column 4 and STRdust ignores it, and `src/motif.rs` is an 18-line stub with
`unimplemented!()`.

### The threshold is not where the length error is

`is_similar_to_ref` feeds only `determine_genotypes` (`src/vcf.rs:403-455`) — the GT field.
It does **not** touch the consensus, so `RB`, `FRB`, `MRL` and the ALT sequence are
unaffected by it. That has two consequences worth settling before any work starts:

- **Raising the threshold fixes the genotype class and leaves the reported length wrong.**
  A locus whose truth is `0/0` and whose consensus is 3 bp long still reports `RB=3`; only
  the GT flips from `0/1` to `0/0`. Any length-based benchmark — the paper's, and
  `misc/benchmark_truth.py`, which scores `RB` — would show **no improvement at all**.
- **The paper's own description points at the consensus, not the threshold.** The Mendelian
  errors at `0/0` are "driven by a higher rate of motif-length-off calls", i.e. the
  consensus length is off by roughly one repeat unit. That is stutter in the POA consensus,
  not a classification threshold.

So the more useful fix is one step earlier: **make the allele length itself right** — snap
the consensus to the mode/median of the spanning-read length distribution, or fit a stutter
model as Medaka Tandem (error model estimated from the reads) and LongTR (HipSTR's stutter
model) do; those are the two tools at the top of the benchmark. Get the length right and the
genotype follows for free, because `is_similar_to_ref` then sees an exact match. Fixing the
threshold alone buys a better-looking GT column and nothing a benchmark can see.

Proposed direction, in order:

1. **Parse the motif from catalog column 4** and carry it on `RepeatInterval`. Fall back to
   inferring the motif from the reference sequence when the column is absent. Needed by
   everything below, and by the benchmark's off-by-motif class.
2. **Correct the allele length**, not just its classification — snap the consensus to the
   mode/median of the per-read length distribution at the locus, behind a CLI knob so it can
   be swept. This is the change a length-based benchmark can actually see.
3. **Then revisit `reference.len() / 20`** with a motif-aware threshold — roughly
   `max(motif_len - 1, reference.len() / 20)`. After 2 this should matter much less; if it
   still does, that is informative in itself.
4. **Decide the question the genotyper is asking.** Both QUICKREF and `is_similar_to_ref`
   currently ask "is this identical to reference?" in different ways. The question that
   matters is whether the *read length distribution* is centred on the reference. A median-
   or mode-based statement is the natural form of that, and would also fix QUICKREF (§4).
5. **Re-run the benchmark** (#23) before and after, stratified by motif size — homopolymers
   are where the current threshold is most degenerate and where the paper reports <40%.

Do not tune QUICKREF against full genotyping. The two disagree because **both** are
noise-intolerant, in opposite directions — QUICKREF says "any read disagrees → not
reference", the genotyper says "consensus isn't identical → not reference". Matching one to
the other just picks which broken reference to agree with. Treat it as one problem.

---

## 4. On PR #30 (open) — the QUICKREF coordinate bug

### 4.1 Is #30 a sufficient fix? No.

It fixes the *false positives* — QUICKREF wrongly calling variant loci `0/0` — and does
that completely (212 wrong calls → 0). It does nothing about the *false negatives*, which
is the failure mode you suspected: QUICKREF almost never fires at loci that genuinely are
homozygous reference. After the fix it fires on **1 locus in 2000 (0.1%)**, where the paper
reports **83.6% of loci are invariant**. That is three orders of magnitude too low.

The fix in fact makes firing *rarer*, for two reasons: the coordinate correction widens the
inspected window to the repeat's true bounds, and the ±15 padding widens it further, from
11 bp to ~41 bp — three times as much room for indel noise to appear in some read. The
padding cut false `0/0` from 2 to 0, but also cut firing.

Two consequences worth planning around:

- **The criterion, not the coordinates, is the reason QUICKREF does not fire** (§4.4). No
  coordinate or tolerance setting rescues it.
- **#30 in isolation could lower a benchmark score.** The false `0/0` calls it removes were
  accidentally *masking* the over-calling of #31: at ~11% of loci QUICKREF was emitting
  `0/0` where full genotyping would have emitted a spurious `0/1`, and at the majority of
  those loci the truth is `0/0`. Removing a wrong mechanism that was accidentally right
  some of the time is still the correct change, but do not read a small drop as a
  regression — score #30 and #31 together, not in sequence.

### 4.2 Mechanism

`calculate_all_length_diff_from_cigar` compared a **1-based** `start` against a **0-based**
`reference_position`, with strict `<` at both ends:

```rust
Cigar::Ins(len) if start < reference_position && reference_position < end => ...
```

The inspected window therefore began two bases into the repeat and stopped one short of its
end. Aligners overwhelmingly place a repeat's indel on the repeat's **first base**: at
chr1:20861973 (an 11 bp interval) 26 of 46 spanning reads carry an indel at exactly offset 0,
and none were counted. `has_variation` stayed false and the locus was emitted `0/0` with
`SUP=.,.`.

Both modes are affected — QUICKREF runs before either. The batched path passed the interval
unpadded and was fully exposed; the non-batched path pads ±15, which masked the off-by-one
and is why it survived.

### 4.3 Effect of the fix

| | QUICKREF fires | false `0/0` |
|---|---|---|
| before | 231 / 1999 | **212 (92%)** |
| coordinates fixed | 10 / 1999 | 2 |
| + ±15 padding in the batched path | **1 / 1999** | **0** |

Genotypes match a QUICKREF-disabled run at 1998/1999 loci. On the chr7 test data QUICKREF
goes from 427/800 to 35/800. Cost: sensitive 300 loci 69 s → 91 s CPU; fast 2000 loci
7.1 s → 9.1 s — the speed QUICKREF was buying was largely fraudulent.

### 4.4 Still open inside QUICKREF, not addressed by #30

1. **No minimum read count in the batched path.** The non-batched path requires ≥5 reads
   checked; the batched one has no guard, so a locus covered by one read, or only by reads
   clipping its edge, can still be called `0/0`.
2. **`SUP=.,.`** — a QUICKREF record carries no read support, so a `0/0` backed by 40 reads
   is indistinguishable from one backed by 1. Nothing downstream can filter it.
3. **Different criteria between the paths** — batched requires `diff == 0` exactly,
   non-batched allows `|diff| ≤ 3`. A read with an equal-sized insertion and deletion inside
   the interval nets to zero and looks reference.

### 4.5 The criterion itself is wrong — this is why it does not fire

QUICKREF fires only if *every* overlapping read shows exactly zero net length difference —
a conjunction over reads, not a statement about the genotype. At the median locus only **25%
of reads** are clean; requiring all ~32 to be simultaneously clean is `0.25^32`.

| fraction of reads clean | loci |
|---|---|
| 100% | 1 (0.1%) |
| ≥90% | 9 (0.5%) |
| ≥75% | 44 (2.2%) |
| ≥50% | 248 (12.4%) |

It therefore measures read quality, not genotype, and gets **worse with depth** — exactly
backwards for a genotype test:

| depth | every read clean | median read diff == 0 |
|---|---|---|
| 16–25x | 0.3% | 43.5% |
| 26–35x | 0.0% | 53.1% |
| 36x+ | 0.0% | 49.1% |

Replacing the conjunction with a consensus statement lands the firing rate where it should be:

| criterion | fires |
|---|---|
| every read exactly 0 (current) | 1/2000 (0.1%) |
| ≥80% of reads exactly 0 | 23/2000 (1.1%) |
| **median read diff == 0** | **1004/2000 (50%)** |
| **modal read diff == 0** | **1195/2000 (60%)** |

50–60% is the same order as the paper's 83.6% invariant rate, rather than three orders of
magnitude below it. Note that a median-based QUICKREF would fire on 1004 loci of which full
genotyping calls 767 non-reference — that is not a 76% error rate, it is #31 showing up from
the other side.

### 4.6 Also: the careful QUICKREF is dead code

Over 2000 loci the **non-batched QUICKREF fired 0 times, and the non-batched path was
entered 0 times.** It is reached only when a target has no usable reads, in which case it
reports no coverage rather than QUICKREF. So the permissive, ≥5-read, `|diff| ≤ 3`
implementation has never run; the strict one is the only one that ever executes.

### 4.7 Decision needed

On a dense STR catalog a *correct* QUICKREF fires on 1 locus in 2000 and earns nothing. It
may still pay off on panels with more genuinely reference loci. The options are: merge #30
and leave it near-inert; merge #30 and then rebuild the criterion as a median/mode statement
(preferably as part of #31, since it is the same question); or drop QUICKREF entirely.

---

## 5. Open: #24 — ALT can describe bases outside `[POS, END]`

`POS`, `END` and `REF` come from the catalog interval only; nothing derives them from the
extracted sequence, so the same BED entry yields identical coordinates in every sample
regardless of that sample's variants. Only `ALT`, `RB`, `FRB` and `MRL` move.

Both modes deliberately fold in indels the aligner placed near but outside the interval
(±30 bp in `sensitive` via `parse_cs`, `--fast-flank` in `fast`). The record is therefore not
a faithful "replace REF at these coordinates with ALT" statement: `bcftools consensus` would
reconstruct a haplotype of the right *length* with up to ~30 bp of bases out of place. At
chr1:20463511 (a 10 bp interval) both modes report a ~920 bp ALT.

**Recommendation: document, do not chase.** Widening `POS`/`END` to cover the flank would be
spec-exact but costs the stable-coordinates property. Shrinking the flank costs real
sensitivity — `--fast-flank 0` missed a 226 bp expansion entirely at one measured locus.
STRdust's consumers read `RB`/`FRB`/`MRL` rather than reconstructing haplotypes.

Middle path if wanted: a **`FLANKINDEL`** INFO flag on the affected records, so anyone doing
reconstruction can skip them. Cheap, and honest about which records are approximate.

### VCF state, verified with bcftools

Good: 2000 records parse with no warnings; every emitted INFO/FORMAT key is declared;
`bcftools norm -c w` reports **zero REF mismatches**; no duplicate ALT alleles and none
identical to REF; `END == POS + len(REF) - 1` on every record.

Fixed in #29: a bare trailing semicolon in INFO on 231/2000 records (`END={end};{somatic}`
with an empty `somatic`, hitting every QUICKREF and missing-genotype record); `|` vs `/`
separator misuse; and the missing `##fileformat=VCFv4.2` declaration, plus a README note
describing the ALT. Consequence documented at the time: QUICKREF records now use `/` even in
a phased run, since they never fetch reads and so have no phase set.

Remaining, not a bug: `bcftools norm` realigns 1557/2000 records because STRdust emits the
full repeat sequence rather than a parsimonious left-aligned allele. That is a SHOULD, not a
MUST, and is what TRGT does too — but see §2(4) for its benchmark consequence.

---

## 6. Open: #23 — no accuracy benchmark

`--mode fast` has only ever been compared against the `sensitive` path, which measures
agreement and silently assumes the alignment path is right. `sensitive` has never been
benchmarked against a truth set either.

Current numbers, 2000 loci, single thread:

| | sensitive | fast |
|---|---|---|
| CPU | 1016 s | 13 s |
| peak RSS | 643 MB | 289 MB |

`MRL` agrees within 2 bp at 95.7% of alleles (median difference 1 bp); 14 of 1760 loci differ
by more than 10 bp, nearly all in one region where a ~900 bp insertion sits between four
overlapping catalog entries and the modes attribute it to different ones.

Needed:

- a truth set — GIAB HG002 TR benchmark VCF, or assembly-derived allele lengths from the
  HG002 trio assemblies;
- per-locus allele-length concordance for both modes against truth, stratified by allele
  length, motif size and locus complexity;
- specifically the case `fast` should lose on: expansions long enough that reads clip rather
  than span, where `fast` drops the read and falls back to `sensitive`. On the sample tested
  the fallback never triggered, so it is untested against real data;
- re-scoring with an off-by-one-motif tolerance, to reproduce the paper's 72.56/94.93 split
  and confirm #31 is the cause.

`misc/benchmark_truth.py` (merged in #28) scores both modes against a truth set and is the
starting point.

---

## 7. Fixed, for reference

- **#19** — (closed by #26) `find_insertions` (`src/genotype.rs`) iterated every mapping minimap2 returned,
  so a multi-mapping read contributed several entries to `insertions`, each counting as a
  separate read for `SUP`, the `--support` threshold, clustering and the consensus. Observed
  `SUP=25,25` at a locus with 37 spanning reads (17/17/3).
- **#20** — `call.rs::process_batch` required only *overlap* and applied no MAPQ filter,
  while `parse_bam::get_overlapping_reads` requires a read to *span* the locus and drops
  MAPQ 0. The permissive path is the one that actually runs. Closed by #27, which added
  `--mapq` and applies it in both paths.
- **#21** — (closed by #26) on a `--haploid` chromosome, `genotype_with_extracted_reads` reads pooled reads
  from `phase0`, but `process_batch` routed reads by HP tag whenever `--unphased` was unset,
  regardless of haploidy — so HP-tagged reads were silently dropped.
- **#22** — (closed by #27) `Allele::from_consensus` computed `ref_len = end - start` while REF is fetched
  inclusive on both ends (`end - start + 1` bases), so `RB` was consistently `FRB - len(REF) + 1`.
  `MRL` and `DBSCAN_RB` carried the same offset; `FRB` was unaffected.
- **#29** — the three VCF conformance issues described in §5.

---

## 8. The structural cause

#19, #20, #21 and the QUICKREF coordinate bug all have the same shape: **the batched path is
a less careful copy of the non-batched path, and it is the one that actually runs.** The
non-batched path is more careful in each case (spanning + MAPQ filter, pooled haploid reads,
±15 padding, ≥5-read guard, `|diff| ≤ 3` tolerance) and, per §4, is very nearly dead code.

That is a code-structure problem, not four coincidences. The two paths should share one
read-selection implementation and one reference-check implementation. Doing that would have
prevented every one of these bugs, and is worth doing before adding more logic to either.

---

## 9. Suggested order of work — SUPERSEDED

**This ordering predates the benchmark. Use §14.4 instead.** Kept because the reasoning
behind items 5 and 6 still holds, and because the difference between this list and §14.4 is
itself the record of what the truth set changed. Item 1 has been done (§13, §14); item 2 has
dropped from first to fourth priority; item 4's premise was wrong (§14.3).

1. **Run the benchmark (#23) before changing anything else** — it is the only thing that
   settles whether #31 is really the cause, and every measurement above is one sample scored
   against STRdust's own full genotyping, not against truth. See §10 for the configurations
   worth running in one pass.
2. **#31**: parse the motif from catalog column 4, replace the `len()/20` threshold with a
   motif-aware one behind a CLI knob. Highest expected benchmark impact.
3. Merge **#30** — correct on its own terms, but score it *together with* #31 rather than
   before it (§4.1).
4. Rebuild the **QUICKREF criterion** (§4.5) as part of 2 — median or mode of the read
   length distribution, plus a minimum read count and real `SUP` values — or remove QUICKREF
   entirely. Do not tune it against full genotyping; both sides are noise-intolerant in
   opposite directions.
5. **Unify the batched and non-batched paths** (§8) — one read selection, one reference check.
6. **#24**: add the `FLANKINDEL` flag if haplotype reconstruction ever matters downstream.

---

## 10. What to measure when the benchmark runs

The confounds above are separable, but only if the run covers the right configurations. In
one pass over the truth set:

| run | flags | isolates |
|---|---|---|
| A | current `main` (QUICKREF buggy) | the baseline the paper scored |
| B | `--alignment-all` (QUICKREF off) | full genotyping alone — the #31 signal, uncontaminated |
| C | `fix/quickref-coordinates` (#30) | whether the fix helps or hurts on truth (§4.1) |
| D | #30 + a consensus-length fix (§3) | the #31 fix — note a *threshold*-only fix moves nothing here |
| E | D + a median-based QUICKREF | whether QUICKREF can earn its place at all |

Score each with **and without** an off-by-one-motif tolerance. Reproducing the paper's
72.56% → 94.93% split on run A is the single most valuable confirmation available: if it
reproduces, #31 is the diagnosis and B/D quantify the fix; if it does not, the diagnosis
needs revisiting before any of this work is worth doing.

Stratify by motif size, and separate mononucleotides — the paper reports STRdust below 40%
there, and it is exactly where `reference.len() / 20` degenerates to an exact-match test.

**The script reports two axes, and a change can move one without the other.**

- **Length** (`exact_pct`, `within_1bp_pct`, …) scores `RB` against the truth allele
  lengths, per locus, taking the worse of the two alleles. Same axis as the paper's
  concordance, and the right primary metric. A consensus-length fix moves this.
- **Genotype class** (`gt_match_pct`, plus `gt=0/0` / `gt=0/1` / `gt=1/1` / `gt=1/2` rows
  and a `.genotype_confusion.tsv`) scores the GT field. A change to `is_similar_to_ref`
  moves *only* this, since that function never touches the consensus (§3).

The #31 signature to watch for is a `gt=0/0` row with high `within_1bp_pct` but low
`gt_match_pct`: lengths near-right, genotypes wrong. That is the shape the paper reports —
STRdust alone below 90% concordance at `0/0`.

It is also not a reimplementation of the paper (different truth set, positional rather than
minimum-deviation allele pairing, tolerance bands rather than Match/off-by-1/off-by-motif
classes, no Levenshtein metric). That is fine for tracking improvement over time — it just
means its absolute numbers should never be quoted against the paper's.

Also worth recording per run, since they are cheap and diagnostic:

- **fraction of loci called non-reference** (against 83.6% invariant in the paper's HPRC set);
- **QUICKREF firing rate and its false-`0/0` count against truth** — this is the number that
  says whether §4.5 is right, and the only version measured so far was scored against
  STRdust's own genotyping, which §4.5 argues is itself wrong at those loci;
- **`SUP` distribution at `0/0` calls**, which is currently `.,.` for QUICKREF records
  (§4.4) and therefore unfilterable.

Both modes should be run: `sensitive` has never been scored against a truth set either, and
the fast/sensitive agreement numbers in #23 assume the alignment path is right.

---

## 11. Reproducing the measurements

The scratchpad from the original session is gone; the inputs are not. Everything in §3–§4
came from:

```bash
REF=/home/wdecoster/database/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
C=/home/wdecoster/testdata/rr_NHB05_082.cram      # 30x ONT

# the 2000-locus catalog subset: chr1:20-25 Mb from the inquiSTR catalog, first 2000 entries
zcat /home/wdecoster/repos/inquiSTR/repeat_catalog_v1.hg38.1_to_1000bp_motifs.bed.gz \
  | awk -F'\t' '$1=="chr1" && $2>20000000 && $3<25000000' | cut -f1-3 > real.bed
head -2000 real.bed > real2000.bed                # spans chr1:20-21.4 Mb

./target/release/STRdust --mode fast -R real2000.bed $REF $C > qr_on.vcf
./target/release/STRdust --mode fast --alignment-all -R real2000.bed $REF $C > qr_off.vcf
```

`qr_off.vcf` (QUICKREF disabled) is the "full genotyping" reference used throughout —
note §4.5's caveat that it is not truth.

The firing-rate tables in §4.5 came from replicating
`calculate_all_length_diff_from_cigar` in Python over the padded window, rather than from
STRdust itself, so criteria could be swept without rebuilding:

```python
PAD = 15   # QUICKREF_PADDING in src/call.rs
def net_diff(r, lo, hi):          # lo, hi = 0-based padded interval
    d = 0; ref = r.reference_start
    for op, l in r.cigartuples:
        if op in (0, 7, 8): ref += l
        elif op == 2:
            if lo <= ref <= hi: d -= l
            ref += l
        elif op == 3: ref += l
        elif op in (1, 4):
            if lo <= ref <= hi: d += l
    return d
# per locus: diffs = [net_diff(r, s-PAD, e+PAD) for r in bam.fetch(c, s, e)
#                     if r.mapping_quality >= 10 and not r.is_secondary and not r.is_supplementary]
# current rule: all(d == 0 for d in diffs)   vs   median(diffs) == 0 / mode(diffs) == 0
```

The tolerance sweep in §3 used a throwaway branch that read the threshold from the
environment, so the whole sweep ran without recompiling per value:

```rust
// src/call.rs, in place of `if diff != 0 {`
let tol: i64 = std::env::var("QR_TOL").ok().and_then(|v| v.parse().ok()).unwrap_or(0);
if diff.abs() > tol {
```

---

## 12. Loose end

The last question of the original session was never answered: **the §4.5 measurement (the
conjunction-over-reads analysis, the depth table and the median/mode firing rates) was never
filed anywhere.** Issue #31 carries only the tolerance sweep. Either add it to #31 or open a
separate issue for the QUICKREF criterion — otherwise this document is its only record.

---

## 13. Benchmark results, 2026-09-08 (v0.16.0b, 50k loci)

First real measurement against a truth set. HG002, GIAB `HG002_GRCh38_TandemRepeats_v1.0.1`
truth VCF over `adotto_TRregions_v1.2.1` loci, 50,000 randomly sampled loci genome-wide,
`--mode sensitive` equivalent (v0.16 has one path), 36,272 s CPU.

`--confident-bed HG002_GRCh38_TandemRepeats_v1.0.bed.gz`, `--sample-loci 50000 --seed 1`,
`--threads 24`. The confident BED is what makes the `0/0` class trustworthy: a catalog locus
with no truth record counts as homozygous reference only when GIAB actually covers it. (The
BED is v1.0 against a v1.0.1 VCF — worth a glance, but the numbers below behave as if the
regions match.) Reuse `--seed 1` for every later run so the same 50,000 loci are scored and
the deltas are real.

### 13.1 The harness reproduces the paper

`RB` in v0.16 is one base too high on every allele (#22) — verified here directly, at
`RB - (FRB - len(REF)) == 1` on all 63,524 alleles. Correcting for it:

| | this benchmark (corrected) | paper (v0.16.0) |
|---|---|---|
| exact / Match | 62.8% | — |
| **within 1 bp** | **73.8%** | **72.56%** (Match + off-by-1) |
| **within 5 bp** | **96.5%** | **94.93%** (off-by-motif permitted) |

Different sample, different truth set, different catalog, and both headline numbers land
within ~1.5 points. The harness is measuring what the paper measured, and §2's diagnosis can
be trusted. **Note the summary TSV reports 33.3% exact / 66.5% within-1bp for a v0.16 run —
those columns are shifted by the `RB` bug and must be corrected before comparing to `main`.**

### 13.2 The dominant error is under-calling, not over-calling

| truth | called | loci | |
|---|---|---|---|
| 0/0 | 0/0 | 32,598 | 79.2% |
| 0/0 | non-ref | 6,790 | 16.5% |
| variant | **0/0** | **3,210** | **36.3% of variant loci** |
| variant | non-ref | 4,893 | |

Truth is 82.3% `0/0`, matching the paper's 83.6% invariant. But recall at variant loci is
poor: over a third are called homozygous reference.

**It is not QUICKREF.** Of the 3,210 false hom-refs, 2,676 (83%) are ordinary genotyped
calls and only 534 are QUICKREF records. Of the genotyped ones, 2,563 (96%) reported a
**non-zero `RB`** — STRdust measured a length difference and called `0/0` anyway. So this is
`is_similar_to_ref`, not read selection and not the consensus failing to see the allele.

### 13.3 The threshold scales with locus length; the noise does not

`is_similar_to_ref` forgives `len(REF)/20`. Every one of those 2,676 false hom-refs has all
alleles within `len(REF)/20`. Their median locus is 130 bp (all loci: 73 bp), where the
threshold forgives 6 bp — and their truth deviations are median 2 bp, p90 6 bp. Real
variation is being absorbed by a tolerance sized for the interval rather than for the error.

| locus length | loci | threshold | over-call | under-call |
|---|---|---|---|---|
| 20–50 bp | 4,182 | 2 bp | **43.8%** | 25.1% |
| 50–100 bp | 32,987 | 3 bp | 17.2% | 26.8% |
| 100–200 bp | 6,842 | 6 bp | 6.7% | 34.0% |
| 200–500 bp | 2,601 | 14 bp | 0.8% | 55.3% |
| >500 bp | 1,934 | 34 bp | 0.7% | **76.0%** |

This is one formula failing in both directions. §3 described only the short-locus end
(threshold 0, over-calls); the long-locus end is the larger effect on a realistic catalog.

(This table counts every scored locus, QUICKREF records included. `misc/compare_runs.py`
prints the same table restricted to *genotyped* loci, so its over-call rates are higher —
70.2% rather than 43.8% in the 20–50 bp bin for v0.16. Same effect, different denominator.)

### 13.4 But no threshold rule fixes it

Simulating the reference-vs-alt decision as a pure length rule over 32,104 genotyped loci
(edit distance's dominant term):

| rule | over-call | under-call | balanced accuracy |
|---|---|---|---|
| current `len/20` | 11.9% | 55.4% | 66.4% |
| constant 1 bp | 37.3% | 14.0% | 74.4% |
| **constant 2 bp** | 21.4% | 29.6% | **74.5%** |
| constant 3 bp | 13.0% | 38.5% | 74.3% |
| constant 5 bp | 4.9% | 58.6% | 68.3% |
| `min(4, len/20)` | 13.1% | 45.6% | 70.6% |

The best achievable by retuning the threshold is ~74.5%, against 66.4% today — worth about
8 points, and then it stops. The ceiling is set by the consensus, not the classifier.

### 13.5 The consensus is insertion-biased — this is the root cause

Consensus length error at 47,419 alleles of truth-homozygous-reference loci, where the
correct answer is exactly 0:

| deviation | share |
|---|---|
| 0 bp | 67.1% |
| +1 bp | 8.8% |
| +2 bp | 9.1% |
| +3 bp | 4.6% |
| +4 bp and beyond | 6.1% |
| any negative | 3.8% |

**29.1% too long against 3.8% too short.** That asymmetry is not random error, and the bump
at +2 exceeding +1 looks like whole repeat units being added — consistent with the paper's
"off-by-motif" finding. Two candidate mechanisms, both testable:

1. **Flank folding (#24).** Both paths deliberately fold in indels the aligner placed up to
   ±30 bp outside the interval (`parse_cs` in the sensitive path). At a locus that is truly
   reference, a nearby unrelated insertion is added to the allele. This was documented and
   deliberately not chased; it now looks like a primary accuracy problem, not a cosmetic
   spec issue. Testable by re-running with the flank reduced and re-measuring this table.
2. **POA stutter.** The consensus adds a repeat unit under ONT indel noise. Testable by
   comparing the consensus length against the mode of the spanning-read length distribution.

### 13.6 What this changes

- **#31 needs rewriting.** "Over-calls alternate alleles at short repeats" is one half of a
  symmetric failure, and the smaller half on a realistic catalog. The title should be that
  the reference/alternate decision is sized by interval length rather than by the error.
- **The fix order in §9 stands, but for a better reason.** Fixing the consensus length is
  not merely the change a benchmark can see (§3) — it is the only change that lifts the
  ceiling. Retuning the threshold is worth ~8 points as a stopgap; a capped threshold such
  as `min(4, len(REF)/20)` is a one-line change if something is wanted before the real fix.
- **#24 is promoted** from "document, do not chase" to a candidate root cause of the
  insertion bias, and the flank-reduction experiment in 13.5 should be run early — it is
  cheap and it discriminates between the two mechanisms.
- ~~**QUICKREF is a side issue at this scale.**~~ **Wrong — see §14.2.** This was based on
  v0.16's *false* hom-refs (534 loci) without looking at how often QUICKREF fired correctly
  (16,471 loci, 96.8% right). QUICKREF is one of the better-performing components, not a
  side issue.

---

## 14. Benchmark results, 2026-09-09 (`main`, run A)

Same 50,000 loci, same seed, `--modes sensitive,fast`. `main`'s `RB` needs no correction —
verified at `RB - (FRB - len(REF)) == 0` on 90,752 alleles, so #27 is doing its job and the
v0.16 numbers below are the corrected ones.

| run | exact | ≤1 bp | ≤5 bp | GT match | recall | precision | no-call | CPU |
|---|---|---|---|---|---|---|---|---|
| v0.16 (corrected) | 62.8% | 73.8% | 96.5% | 74.9% | 62.2% | 43.8% | 1,425 | 36,272 s |
| **main sensitive** | **52.2%** | **66.2%** | 96.2% | **68.0%** | 64.8% | **35.8%** | 3,582 | 40,215 s |
| **main fast** | **72.9%** | **85.1%** | **98.6%** | **81.9%** | 55.5% | **60.8%** | 3,546 | 21,561 s |

Two results, both large, and neither expected.

### 14.1 `--mode fast` is markedly *more accurate* than `sensitive`, not just faster

20 points of exact concordance, 19 points of length concordance, 14 points of genotype
agreement, and 25 points of precision — at roughly half the CPU. The decisive number is the
over-call rate among **genotyped** homozygous-reference loci, which strips out every other
difference between the paths:

| path | over-call rate at genotyped hom-ref loci |
|---|---|
| main sensitive (POA consensus) | 26.6% |
| main fast (length straight from the CIGAR) | **8.2%** |

Reading the length off the existing alignment is three times more accurate than building a
POA consensus and measuring that. This is §13.5's insertion bias identified by mechanism:
**the POA consensus is the root cause**, and `fast` avoids it entirely.

`fast` trades recall for precision (55.5% vs 64.8%), so it is not uniformly better on every
axis, but it wins decisively on both combined metrics. The obvious question this raises is
whether `sensitive` should remain the default, or whether the alignment path should be kept
only for the case it is genuinely needed — expansions long enough that reads clip rather
than span.

### 14.2 `main` sensitive looks like a regression against v0.16 — but is not

10 points of exact concordance and 7 points of genotype agreement worse than a release from
February. The cause is not the genotyper:

| run | QUICKREF fires | of those, truly 0/0 |
|---|---|---|
| v0.16 (tolerant, non-batched path) | **16,471 (32.9%)** | **96.8%** |
| main (strict, batched path) | 1,042 (2.1%) | 98.8% |

v0.16 short-circuited **a third of the benchmark** with an answer that was right 96.8% of
the time. `main` short-circuits 2.1%, so 31% more loci are routed into full genotyping,
where the over-call rate is ~27%. Per genotyped locus `main` is slightly *better* than
v0.16 (26.6% vs 29.0% over-call); it only looks worse because the good path stopped
catching anything.

### 14.3 This inverts §4 and §13.6

- **QUICKREF is not a side issue, and the "dead code" is the good implementation.** §4.6
  established that the tolerant non-batched QUICKREF never runs in `main`. It is now clear
  what that cost: a heuristic firing at 32.9% with 96.8% precision was replaced by one
  firing at 2.1%. The tolerant criterion (`|diff| ≤ 3`, first 25 reads, ≥5 required, padded
  coordinates) is not crude — at these loci it is **more accurate than full genotyping**,
  which manages ~73% precision on the same question.
- **§4.5's tolerance sweep was measured against the wrong reference.** It scored QUICKREF
  against STRdust's own full genotyping and concluded tolerance 3 was "9 of 26 wrong". Against
  truth, tolerance 3 is right 96.8% of the time across 16,471 loci. The sweep was measuring
  disagreement with a genotyper that is itself wrong ~27% of the time at these loci.
- **§13.6 said QUICKREF is a side issue at this scale. That was wrong**, and it was wrong
  because it was based on `main`'s firing rate rather than on what a correct QUICKREF can
  reach.
- **#30 remains correct but is now clearly insufficient.** Fixing the strict criterion's
  coordinates leaves it firing on ~2% of loci. The larger win is giving the batched path the
  tolerant criterion and the padding the non-batched path already had.

### 14.4 Revised priorities

1. **Run B (`--alignment-all`)** — with QUICKREF removed entirely, this separates the
   genotyper's own accuracy from QUICKREF's coverage, and gives the clean baseline for the
   consensus work. It is now the most informative single run available.
2. **Restore a tolerant QUICKREF in the batched path** (padding + `|diff| ≤ 3` + a minimum
   read count), on top of #30's coordinate fix. On this evidence that is worth more than any
   threshold change, and it is a small change to code that already exists.
3. **Fix the consensus** (§13.5) — `fast` shows the achievable gain: 26.6% → 8.2% over-call.
   Either estimate the length from the read distribution rather than the POA consensus, or
   reconsider which path is the default.
4. **Then** the `is_similar_to_ref` threshold (§13.3/13.4), which is worth ~8 points on top.
5. Report `SUP` on QUICKREF records — at 33% firing, `SUP=.,.` on a third of the output is
   no longer a cosmetic problem.

---

## 15. Where this stands, and what to do when run B lands

Written 2026-09-09 while run B is executing, so the thread is not lost.

### 15.1 Tooling now in the repo

- `misc/benchmark_truth.py` — scores one build against a truth set. Since this
  investigation: `--no-mode-flag` (drive a pre-v0.21.0 binary, which has no `--mode`),
  `--repeat-bed` (truth from a general call set: each repeat locus is one target, its truth
  alleles are the summed per-haplotype length change of the variants inside, and a covered
  locus with no variant is homozygous reference), `--sample-loci`/`--seed` (a representative
  random subset instead of the first N), a genotype-class axis (`gt_match_pct`, `gt=` strata,
  `<prefix>.genotype_confusion.tsv`), and a warning when `--confident-bed` is missing.
- `misc/quickref_paired.py` — follows every QUICKREF firing in one run into another,
  splitting its contribution into rescues and corrections (§16.2).
- `misc/benchmark_truth.py --truth-cache PATH` — build the truth set once instead of once
  per run. Worth it for smoke runs, noise for the 50k series (§18.4).
- `misc/compare_runs.py` — **new; reproduces every table in §13 and §14 in one command.**
  It measures each run's `RB` offset from its own VCF, so the #22 correction is applied
  automatically rather than remembered:

  ```bash
  python3 misc/compare_runs.py v0.16=v016_50k main=main_50k
  ```

  Neither script is covered by tests, and neither is run by CI.

### 15.2 The runs so far

All on 50,000 loci, `--seed 1`, so every run scores the identical loci.

| run | what | state |
|---|---|---|
| v0.16 | `v016_50k` | done, §13 |
| A | labelled `main`, sensitive + fast — **actually the #30 build**, §23 | done, §14 |
| B | `main` + `--alignment-all`, sensitive (`noqr_50k`) | done, §16 |
| C | `fix/quickref-coordinates` (#30), sensitive (`pr30_50k`) | valid, = A, §23 |
| C2 | `pr30_50k_again` — **the first run of actual `main`** | valid, §23 |
| D | the junction-inflation matrix: `minlen` × flank window, 4 runs on a shared 10k sample | designed, §20 |

### 15.3 What B is for, and what to expect

B removes QUICKREF entirely, so every locus is genotyped. It is the clean baseline for all
consensus work: after it, any change to the consensus can be measured without QUICKREF's
coverage moving underneath the result.

Expect B to be **slightly worse than A sensitive**, not better. QUICKREF covers 2.1% of loci
in `main` at 98.8% precision, and full genotyping is right ~73% of the time on the same
loci, so removing it should cost a few tenths of a point overall. **If B differs from A
sensitive by more than ~1 point, something other than QUICKREF is involved and is worth
chasing before continuing.**

First command when the files arrive:

```bash
python3 misc/compare_runs.py v0.16=v016_50k main=main_50k noqr=main_noqr_50k
```

Then read, in order: the concordance table (does B sit just below A sensitive?), the
consensus-bias table (B is the honest measurement of the POA consensus, uncontaminated by
QUICKREF), and the threshold sweep (B is the right run to size the `is_similar_to_ref`
change against, since every locus went through it).

### 15.4 Then, in order

1. **C: run PR #30** — the same command as A with `--modes sensitive` on the
   `fix/quickref-coordinates` branch. Settles §4.1's open question empirically. Expect it
   to land between A and B, and close to A.
2. **The tolerant-QUICKREF experiment** — the largest identified win (§14.3). Give the
   batched path the criterion the non-batched path already has: `|diff| <= 3`, a minimum
   read count, and the ±15 padding, on top of #30's coordinate fix. v0.16 reached 32.9%
   firing at 96.8% precision with exactly that criterion; `main` reaches 2.1%. Measure with
   a fourth run.
3. **The junction-inflation matrix (run D, §20)** — the flank-folding experiment of §13.5
   and the `--minlen` question of §19, run together because they are two knobs on one
   mechanism: how much stray inserted sequence near the junction is folded into the allele.
   Four runs on a shared 10,000-locus sample. If the +1/+2/+3 excess shrinks with either
   knob, folding is a root cause and #24 moves from "documented, deliberate" to a fix; if
   neither moves it, the bias is POA stutter and the fix is a read-distribution-based length
   estimate.
4. **Decide what `sensitive` is for.** `fast` is more accurate *and* cheaper on this data
   (§14.1). Before changing any default, check the case `fast` is expected to lose: loci
   whose expansion is long enough that reads clip rather than span. The `>200bp` stratum in
   A is the place to look — `fast` scored 57.1% exact there against `sensitive`'s 46.0%, but
   with 56 no-calls against 69, so neither path is convincing at long expansions and the
   stratum only holds 219 loci. A run enriched for long expansions would settle it.
5. **Only then** the `is_similar_to_ref` threshold (§13.4), worth ~8 points, and #31's
   motif work.

### 15.5 Open questions

- **Why is `fast` only 1.9× cheaper here** (21,561 s vs 40,215 s CPU) when #23 measured 78×
  on a dense chr1 catalog? Likely the clipped-read fallback into the alignment path firing
  often on adotto's long loci — which, if true, means `fast`'s accuracy advantage is being
  measured on a mixture of both paths, not on the CIGAR path alone. Worth checking how often
  the fallback triggers before drawing conclusions about either path.
- **Is the GIAB BED/VCF version mismatch harmless?** `HG002_GRCh38_TandemRepeats_v1.0.bed.gz`
  against a `v1.0.1` VCF (§13).
- **#31's title and body still describe the old diagnosis** and should be rewritten once B
  and C are in.

---

## 16. Benchmark results, 2026-09-09 (`main --alignment-all`, run B)

Same 50,000 loci, same seed, `--modes sensitive`, QUICKREF disabled outright so every
locus goes through full genotyping. Prefix `noqr_50k`. Reproduce with:

```bash
python3 misc/compare_runs.py v0.16=v016_50k main=main_50k noqr=noqr_50k
```

| run | scored | exact | ≤1 bp | ≤5 bp | GT match | recall | precision |
|---|---|---|---|---|---|---|---|
| v0.16 (corrected) | 48,575 | 62.8% | 73.8% | 96.5% | 74.9% | 62.2% | 43.8% |
| main fast | 46,454 | 72.9% | 85.1% | 98.6% | 81.9% | 55.5% | 60.8% |
| main sensitive | 46,418 | 52.2% | 66.2% | 96.2% | 68.0% | 64.8% | 35.8% |
| **B (noqr)** | 46,205 | **51.9%** | 65.9% | 96.2% | **67.6%** | 64.9% | 35.5% |

### 16.1 B is 0.3 points below A sensitive — §15.3's decision rule is satisfied

The rule was: expect B slightly worse than A sensitive, and if it differs by more than
~1 point, something other than QUICKREF is involved and must be chased first. Observed is
−0.3 exact / −0.4 GT. **Nothing other than QUICKREF is involved.** B is the clean baseline
for consensus work, and every §13/§14 conclusion about the genotyper stands unchanged.

### 16.2 The paired measurement on QUICKREF's own loci

A and B score the identical loci, so the 1,042 QUICKREF firings in A can be followed into B
one-to-one. This is the direct measurement §4.5 could not make, because it compared QUICKREF
against STRdust's own genotyper rather than against truth.

| on the 1,042 loci QUICKREF answered in A | loci | share |
|---|---|---|
| QUICKREF exact | 1,030 | **98.8%** |
| full genotyping exact (B) | 749 | **71.9%** |
| full genotyping called it, but wrong | 80 | 7.7% |
| full genotyping produced no call at all | 213 | 20.4% |
| *of those, loci QUICKREF got right — a true rescue* | *201* | *19.3%* |
| QUICKREF wrong (all but two of them truth 1/1) | 12 | 1.2% |

Two things follow that were not visible before:

- **QUICKREF's main contribution is rescue, not accuracy.** Of the 1,030 loci it gets right,
  201 are loci where full genotyping fails outright and only 80 are loci it would get wrong.
  Where full genotyping does produce a call it is 749/829 = 90.3% exact on this subset —
  still below 98.8%, but far better than its ~72% over the whole run. QUICKREF is firing
  preferentially on loci the alignment path finds hard.
- **The whole footprint is inside the truth 0/0 class.** The genotype confusion tables are
  identical between A and B outside it (0/1, 1/2 differ by ≤2 loci); 0/0→0/0 drops 28,138 →
  27,822 and 0/0→missing rises 3,199 → 3,400. The 12 wrong firings are the 1/1 class's
  `missing` going 280 → 290. QUICKREF is not silently answering non-reference loci.

This strengthens §14.3 rather than changing it: the tolerant criterion is worth restoring,
and a restored one should be measured on both axes — precision *and* how many otherwise
no-call loci it rescues.

### 16.3 The consensus bias is unchanged, and now measured without QUICKREF

| run | alleles at truth hom-ref loci | exact | too long | too short |
|---|---|---|---|---|
| v0.16 | 47,419 | 67.1% | 29.1% | 3.8% |
| main fast | 73,900 | **86.5%** | **11.4%** | 2.1% |
| main sensitive | 73,872 | 71.9% | 25.6% | 2.5% |
| **B (noqr)** | 75,530 | 72.4% | 25.1% | 2.5% |

B adds the 829 QUICKREF loci it managed to genotype and moves the number by half a point.
**QUICKREF was neither masking nor inflating the consensus bias** — ~25% of alleles at
hom-ref loci come out too long either way, against `fast`'s 11.4%. §13.5 and §14.1 hold on
the uncontaminated run, and the over/under-call table by locus length is identical to A's to
within 0.5 points in every stratum.

### 16.4 The threshold sweep, on the run that should size it

Every locus in B went through `is_similar_to_ref`, so this is the right run to size a change
against (45,376 genotyped loci in A vs 46,205 here).

| rule | over-call | under-call | balanced acc |
|---|---|---|---|
| current: `len/20` | 9.7% | 57.5% | 66.4% |
| **constant 1 bp** | 32.2% | 16.0% | **75.9%** |
| constant 2 bp | 17.7% | 32.4% | 75.0% |
| constant 3 bp | 10.2% | 41.7% | 74.0% |
| constant 5 bp | 3.5% | 61.1% | 67.7% |
| `min(4, len/20)` | 10.6% | 48.2% | 70.6% |

Unchanged from A: ~9.5 points available, constant 1 bp best, and the gain comes from cutting
under-calling at the cost of over-calling. Worth noting against §14.1: the same sweep on
`fast` reaches **85.0%** at constant 1 bp versus 75.9% here. The threshold is worth more on
the more accurate length estimate, which is an argument for fixing the consensus *first* and
sizing the threshold afterwards, not the other way round.

### 16.5 What is next, unchanged from §15.4

B was the gate, and it passed cleanly, so the order in §15.4 stands: **C (PR #30)**, then
the **tolerant-QUICKREF experiment**, then the **flank-folding test** at
`src/genotype.rs:912`. Nothing in B argues for reordering them.

---

## 17. Run C is staged — how to run it and what to do with the result

> **Executed; see §21-§23.** The run was mislabelled and the verdict on #30 inverted twice before the binaries were identified by content.

Written 2026-09-11, at the same point in the cycle as §15: the run is queued, not started.

### 17.1 What is staged

`fix/quickref-coordinates` (PR #30, one commit `6fd8eea`, `MERGEABLE`, **not behind `main`**)
built as a stripped static musl binary, so C differs from A by exactly that commit and
nothing else:

```
~/Downloads/fix_strdust/run_C/STRdust-pr30-musl   md5 0e26dafe248696b780657d4a59794a09
~/Downloads/fix_strdust/run_C/benchmark_truth.py  (the working-tree version, as used for A and B)
```

Built with `make musl` on the branch; `--version` reports 0.21.0, same as `main`, so the
binary name is the only thing distinguishing them — do not rename it.

### 17.2 The command

Copy `run_C/` to the machine that holds the CRAM and reference, then, from the directory
holding the data:

```bash
python3 benchmark_truth.py \
    --binary ./STRdust-pr30-musl \
    --modes sensitive \
    --fasta GRCh38.fa --bam ont_downsampled.cram \
    --truth-vcf HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz \
    --repeat-bed adotto_TRregions_v1.2.1.bed.gz \
    --confident-bed HG002_GRCh38_TandemRepeats_v1.0.bed.gz \
    --sample-loci 50000 --seed 1 --threads 24 \
    --truth-cache truth_50k_seed1.json --out-prefix pr30_50k
```

Identical to A except for the binary, `--modes sensitive` (no `fast`, which #30 does not
touch), the prefix, and `--truth-cache` (§18). **`--sample-loci 50000 --seed 1` must not
change** — it is what makes all four runs score the same loci. Copy `pr30_50k.*` back to
`~/Downloads/fix_strdust/`, and keep `truth_50k_seed1.json` on the machine: every later run
in the series reuses it.

**Check before trusting the result.** All three completed runs produced a byte-identical
target list, so C must too:

```bash
md5sum pr30_50k.targets.bed   # 50a099114bd2a67311e9b12d3313ac53
```

A different sum means the truth set moved and C is not comparable to A, B or v0.16 — stop
and find out why rather than reading the concordance numbers.

**Cost.** ~40,000 CPU s of STRdust, plus ~44 s of harness. Every CPU figure in §13–§16 is
`RUSAGE_CHILDREN`, i.e. the STRdust subprocess *only* — the harness's own truth-building
never appeared in them, which is why the wall time of earlier runs ran ahead of the quoted
cost. Until 2026-09-11 that phase was quadratic and cost ~2.2 h on top of every run in
§13–§16; see §18. `--truth-cache` is in the command above because it costs nothing, not
because it saves much: at ~44 s it is now noise against the run itself.

### 17.3 First command when the files land

```bash
python3 misc/compare_runs.py v0.16=v016_50k main=main_50k noqr=noqr_50k pr30=pr30_50k
```

### 17.4 What to expect, and the decision rule

C should land **between B and A, close to A** — A 52.2% exact, B 51.9%, so the window is
narrow and the headline concordance number is *not* where the answer is. #30 fixes which
loci QUICKREF is allowed to fire on; it does not make it fire more. Read these instead, in
order:

1. **The QUICKREF footprint table.** `fires` is the number that matters. If #30 leaves it
   near 2.1%, §4.1's answer is confirmed empirically: the coordinate fix is correct but
   insufficient, and the tolerant-criterion experiment (§14.3, §15.4 step 2) is the real
   work. If it rises materially, #30 does more than expected and the priority order changes.
2. **Precision within `fires`.** A is 98.8% on 1,042 loci. #30 should raise it or hold it.
   A *drop* would mean the coordinate fix admits loci the old bug was accidentally excluding,
   which is worth understanding before merging.
3. **The paired follow-through, the §16.2 analysis re-run with `pr30` in place of `main`.**
   This is the measurement that actually settles #30, because it separates the two things
   QUICKREF does: rescuing loci full genotyping no-calls (201 in A) versus correcting loci it
   calls wrong (80 in A). #30 could plausibly improve the second and leave the first alone.
   Run `python3 misc/quickref_paired.py pr30_50k noqr_50k` — the script is in the repo and
   reproduces §16.2 exactly when given `main_50k noqr_50k`.

Merging PR #30 does not depend on C — it fixes a real coordinate bug either way (§4.2). C
decides *how much of §4 is left afterwards*, which is what sets the next experiment.

### 17.5 State of the tree at handoff

`KNOWN_ISSUES.md`, `misc/compare_runs.py` and `misc/quickref_paired.py` still untracked,
`misc/benchmark_truth.py` still modified, on `main`. None of it is committed and none of it is covered by CI. The branch was
checked out only to build; the repo was returned to `main` afterwards.


---

## 18. Harness: the truth-set phase was quadratic

Found 2026-09-11, from the observation that `[benchmark] reading truth set` ran far longer
than the quoted CPU cost of a run.

### 18.1 What was wrong

`inside_confident` walked the confident-region list from the start of the chromosome for
every query. That is harmless for its original caller, `read_truth`, which walks the truth
VCF and asks about a few thousand records. `read_truth_over_repeats` — the `--repeat-bed`
path every run in §13–§17 uses — asks once per *catalog* locus instead: ~1.8M queries against
a comparable number of regions, with the scan growing as the locus moves along the
chromosome.

Measured on a synthetic chromosome at realistic scale (140k regions, 140k loci):

| locus index | cost of one query |
|---|---|
| 0 | 0.009 ms |
| 20,000 | 1.27 ms |
| 60,000 | 4.93 ms |
| 120,000 | 8.20 ms |

which integrates to **~2.2 hours genome-wide**, entirely invisible in every reported CPU
number because those measure `RUSAGE_CHILDREN`. That is the whole of the "reading truth set
takes very long" observation this came from.

### 18.2 The fix

`load_confident_regions` now also stores, per chromosome, the largest region end seen up to
each index. A locus is contained in some region iff the largest end among regions starting
at or before it reaches past the locus — two bisected lookups instead of a walk. The same
synthetic benchmark: **~2.2 h → ~1 s**, a ~10,000× speed-up on one chromosome.

This handles overlapping and nested regions correctly, which a plain "find the one region
containing it" bisect would not.

That ~1 s is the *function*, not the phase. Timed end to end at genome scale (1.78M catalog
loci, 1.5M truth records, synthetic but realistically shaped), the whole truth phase after
the fix is:

| step | time |
|---|---|
| `load_confident_regions` | <0.1 s |
| `read_truth_over_repeats` (VCF walk + catalog + sampling) | 43.8 s |
| **total truth phase** | **~44 s** |

So the phase is now dominated by the plain Python walk over the truth VCF, and there is no
second quadratic hiding behind the first.

### 18.3 Why the earlier results are still comparable

The change is a speed-up, not a semantic change, and that was verified rather than assumed:

- **120,000 randomized queries** against deliberately overlapping and nested region lists,
  old implementation vs new: **0 mismatches**. Empty and missing-chromosome cases match too.
- **End to end on synthetic inputs**, old vs new: the truth loci, the target BED and the
  skip counters are **identical objects**, not merely similar totals.
- **In the field**: v0.16, A and B all wrote the same `targets.bed`
  (`50a099114bd2a67311e9b12d3313ac53`), so C reproducing that sum confirms it on the real
  data. This is the check in §17.2.

`TRUTH_CACHE_VERSION` stays at 1 for exactly this reason: the content of the truth set did
not move, so caches and completed runs remain valid.

### 18.4 `--truth-cache`

The truth set is a pure function of the truth VCF, the repeat catalog, the confident BED and
the sampling — all identical across a benchmark series, and all rebuilt from scratch on every
run until now. `--truth-cache PATH` writes it once and reuses it, keyed on those inputs
(path, size and mtime rather than a content hash, since hashing the VCF would cost more than
the work saved), plus `TRUTH_CACHE_VERSION`. A stale or unreadable cache is reported and
rebuilt, never silently trusted, and the write is atomic so an interrupted run cannot leave a
truncated one behind.

**It is worth much less than it looks, now that the phase is ~44 s.** Across the five queued
50,000-locus runs it saves ~3 minutes against runs that each cost hours — roughly 0.3% of
wall time, i.e. noise. It was written before the phase was timed, when the saving looked
like hours per run.

Where it does still earn its place is **iteration**: on a smoke run (`--sample-loci 2000`,
one chromosome) STRdust finishes in about a minute, so a 44 s truth phase is close to half
the wall clock, and the cache takes it to 0.1 s. Cache write is 0.4 s for a 3.5 MB file, so
it costs nothing to leave on.

Keep it for smoke runs and re-runs; do not expect it to matter for the 50k series. If it
ever needs maintenance, deleting it is the right call — the quadratic fix is what actually
solved the problem.

### 18.5 Not done

Neither the fix nor the cache is covered by a test in the repo — the checks above were run
ad hoc in a scratch directory, against synthetic inputs, and are gone. `misc/` has no tests
at all and is not touched by CI, which is the same gap §15.1 notes for `compare_runs.py`.

---

## 19. `--minlen` means one more than it says

Assessed 2026-09-11. Not filed as an issue yet.

### 19.1 The mismatch

`minlen` has exactly one use, `src/genotype.rs:940` and its `else if` at `:942`:

```rust
if cap[0][1..].len() > minlen && interval_around_junction.contains(&ref_pos) {
    result.push_str(&cap[0][1..]);
```

Strictly greater, so `--minlen N` keeps insertions of **N+1 bp and longer**, while the help
calls it "minimal length of insertion/deletion operation [default: 1]". At the default,
single-base insertions are dropped, not kept.

### 19.2 The `>` is original; the default moved under it

`> minlen` dates from the first commit of `genotype.rs`. What changed is the default, from 5
to 1, in `02d9540` (2025-11-14) — whose subject is "add min_haplotype_fraction parameter to
CLI and related functions". A one-line change, in an unrelated commit, with no test and no
mention in the message.

That matters because the off-by-one was cosmetic at 5 (≥6 rather than ≥5) and is not at 1:
it is now exactly the boundary between keeping and dropping single-base insertions.

### 19.3 Why not simply loosen it

Both obvious fixes — default 0, or `>` → `>=` — make 1 bp insertions count. On run B, the
called allele length at truth homozygous-reference loci (correct answer: 0, n=75,530):

| error | alleles | share | share of all over-calls |
|---|---|---|---|
| −1 | 1,232 | 1.6% | |
| **0** | **54,680** | **72.4%** | |
| **+1** | **5,711** | **7.6%** | **30.1%** |
| +2 | 6,231 | 8.3% | 62.9% cumulative |
| +3 | 3,057 | 4.0% | 79.0% cumulative |

**A single base is the single most common wrong answer.** The mechanism points both ways,
which is why this needs measuring rather than arguing:

- `result.push_str` **concatenates every** qualifying insertion inside
  `flanking ± 30`, so a stray 1 bp insertion in the flank of a read that already carries the
  repeat insertion adds a base to the allele. Loosening feeds the +1 pile directly. This is
  the same mechanism as the flank-folding experiment (§15.4 step 3, #24): `minlen` and the
  ±30 window are the two knobs on the same inflation.
- But `parse_cs` returns `None` when nothing qualifies, and `find_insertions` then drops the
  read **entirely** — it never reaches `insertions.len() < args.support`. A read whose only
  junction insertion is 1 bp is a read carrying a near-complete deletion of the repeat, and
  those currently vanish, costing support and biasing against short alleles. Loosening
  rescues them.

Net effect unknown. Given §16.3 and the table above, the prior is that loosening makes the
dominant error worse, but that is a prior, not a measurement.

### 19.4 What to do: measure first, then change one thing

There are two questions here and they have different answers.

**The semantics are a documentation bug and can be fixed now, behaviour-preserving:** `>=`
with the default set to 2. `--minlen N` would then mean "N and longer" as it claims, while
the effective threshold stays exactly where it is, so v0.16, A, B and C remain comparable.
Changing the *effective* threshold mid-series would silently invalidate the comparison,
which is the one thing this benchmark cannot afford.

**Whether 1 bp insertions should count is a tuning question, and is measured as run D
(§20), not argued.** §19.3 gives a prior, not a result. The measurement is free: `--minlen 0`
on the *current* binary is exactly the loosened behaviour (`len > 0` keeps ≥1), so it needs
no code change and no rebuild.

Order: run D, then set the default from what D shows, then make the `>=` change once — so
the default and the comparison operator move in a single deliberate commit rather than the
default drifting under the operator a second time (§19.2).

### 19.5 Two snags to fix alongside

- The help says "insertion/**deletion** operation", but `minlen` never filters deletions: the
  `'-'` branch only advances `ref_pos`. Harmless — the repeat is excised from the reference,
  so deletions carry no allele length — but the documentation is wrong.
- `src/main.rs:194` uses `args.minlen != 1` to detect "the user set `--minlen`" for the
  `--mode fast` warning. That hardcodes the current default and must move with it; clap's
  value-source API is the robust version.

---

## 20. Run D: the junction-inflation matrix

> **Superseded by §28.** The matrix was replaced by a single-binary sweep once the knobs became CLI flags, and its `--junction-window 0` arm turned out to test a position one base from the junction (§29, and the note at `src/main.rs`).

Designed 2026-09-11, not yet run. This replaces "change 30 to 10 and see" (§15.4 step 3) and
answers §19 by measurement instead of by prior.

### 20.1 Why one experiment and not two

`--minlen` (§19) and the `flanking ± 30` fold window (§13.5, #24) are two knobs on the same
mechanism. `parse_cs` walks a read's CS tag and **concatenates every insertion** that is both
longer than `minlen` and positioned inside the window:

```rust
if cap[0][1..].len() > minlen && interval_around_junction.contains(&ref_pos) {
    result.push_str(&cap[0][1..]);
```

`minlen` sets how *small* a piece may be folded in, the window sets how *far* from the
junction it may come from. Running them separately would confound them; run as a matrix, each
arm isolates one.

### 20.2 What the spectrum already says

`misc/compare_runs.py` now prints the error spectrum base by base, which the aggregate
"too long" number hid:

| run | −1 | +0 | **+1** | **+2** | +3 | +1 as share of over-calls |
|---|---|---|---|---|---|---|
| main sensitive | 1.7% | 71.9% | **7.7%** | **8.4%** | 4.1% | 30.1% |
| noqr (B) | 1.6% | 72.4% | **7.6%** | **8.2%** | 4.0% | 30.1% |
| main fast | 1.5% | 86.5% | **5.4%** | **1.9%** | 2.2% | 47.0% |

**In `sensitive`, +2 is more common than +1; in `fast` it is a third as common.** The CIGAR
path has a +1 tail and essentially no +2 mode; the alignment path has a distinct +2 mode on
top. That is a specific, mechanical thing to explain, and folding a short piece of flank is a
candidate for it. D should account for that +2 mode or show it is untouched by either knob.

### 20.3 The arms

Five runs, same loci, same binary except where noted:

| arm | change | rebuild |
|---|---|---|
| D0 | baseline, current defaults (effective ≥2 bp, window ±30) | no |
| D1 | `--minlen 0` → keeps ≥1 bp, i.e. the loosening §19 asks about | no |
| D2 | `--minlen 5` → keeps ≥6 bp, the pre-2025-11-14 default (§19.2) | no |
| D3 | window ±10: `src/genotype.rs:912`, `30` → `10` | yes |
| D4 | window ±0: same line, `30` → `0` | yes |

D1 and D2 need no code change at all: `> minlen` means `--minlen 0` already *is* the
loosened behaviour, and `--minlen 5` already *is* the old default. Only the window arms need
a rebuild, one each.

D2 matters as much as D1. §19.3's prior is that stray junction sequence inflates alleles; if
that is right, the informative arm is the one that folds in *less*, not more — and it also
tells us whether the undiscussed 5 → 1 change in `02d9540` was a regression.

### 20.4 The run

10,000 loci rather than 50,000: every question D asks is a *difference between arms*, so the
arms need to match each other, not the 50k series. That makes each arm ~1/5 of a full run,
and the matrix roughly the cost of one.

```bash
for arm in "D0 " "D1 --strdust-arg=--minlen --strdust-arg=0" "D2 --strdust-arg=--minlen --strdust-arg=5"; do
  set -- $arm
  name=$1; shift
  python3 benchmark_truth.py \
      --binary ./STRdust-linux-musl --modes sensitive \
      --fasta GRCh38.fa --bam ont_downsampled.cram \
      --truth-vcf HG002_GRCh38_TandemRepeats_v1.0.1.vcf.gz \
      --repeat-bed adotto_TRregions_v1.2.1.bed.gz \
      --confident-bed HG002_GRCh38_TandemRepeats_v1.0.bed.gz \
      --sample-loci 10000 --seed 1 --threads 24 \
      --strdust-arg=--alignment-all "$@" \
      --truth-cache truth_10k_seed1.json --out-prefix ${name}_10k
done
```

`--alignment-all` on every arm: QUICKREF short-circuits loci that neither knob can affect, so
leaving it on would dilute the signal for no benefit. `--truth-cache` genuinely pays here —
five arms sharing one truth set, unlike the 50k series (§18.4). D3 and D4 are the same
command with the rebuilt binary and their own prefix.

Read with:

```bash
python3 misc/compare_runs.py D0=D0_10k D1=D1_10k D2=D2_10k D3=D3_10k D4=D4_10k
```

### 20.5 Decision rules, written before the data

- **`--minlen`.** If D1 raises the +1 mass or lowers exact-at-hom-ref against D0, loosening is
  wrong: keep the effective threshold and fix only the semantics (`>=` with default 2, §19.4).
  If D1 is neutral or better, the read-rescue effect dominates and the default should produce
  ≥1. If **D2** lowers the +1/+2 mass materially, the threshold itself should go up and the
  5 → 1 change was a regression — the strongest available outcome, since it is a one-line
  revert.
- **Always read no-calls next to it.** Tightening `minlen` makes `parse_cs` return `None`
  more often, and `find_insertions` then drops the read entirely, so precision bought this way
  is paid for in support and no-calls. An arm that improves the spectrum while raising
  `no_call` has not obviously won; compare `loci_scored` and recall across arms before
  concluding.
- **The window.** If D3/D4 collapse the +2 mode specifically, folding is its cause, #24 moves
  from "documented, deliberate" to a fix, and the window becomes a tuning parameter worth
  exposing. If the spectrum barely moves, the bias is POA stutter and the fix is a
  read-distribution length estimate rather than either knob (§14.1 shows `fast` already
  achieves 86.5% by not building a consensus at all).
- **Null result is informative.** If no arm moves the spectrum, both knobs are exonerated and
  the consensus itself is the only remaining suspect, which promotes §15.4 step 4.

### 20.6 Why D does not disturb the main series

D is on its own 10,000-locus sample and answers only within-matrix questions, so it can run
before, after or beside C without touching the v0.16/A/B/C comparison. Nothing in D changes
committed code: D1 and D2 are CLI flags, D3 and D4 are throwaway builds. **No default should
change until D has been read** — that is the whole point of queueing it rather than acting on
§19.3's prior.

---

## 21. Run C landed and shows nothing — and that is probably an artefact

2026-09-11. **Do not read C as a result until the binary is confirmed** (§21.4).

### 21.1 What C shows

The §17.2 gate passed: `pr30_50k.targets.bed` is
`50a099114bd2a67311e9b12d3313ac53`, identical to v0.16, A and B, so the truth set is right
and the §18 harness change did not move it on real data.

Everything else is identical to `main` sensitive **to the decimal**:

| run | scored | exact | GT | QUICKREF fires | precision |
|---|---|---|---|---|---|
| main sensitive | 46,418 | 52.2% | 68.0% | 1,042 (2.1%) | 98.8% |
| pr30 | 46,418 | 52.2% | 68.0% | 1,042 (2.1%) | 98.8% |

Not merely the same count — **the same 1,042 loci**, set-identical, zero in either direction.
Across the whole VCF only **3 records of 50,000** differ.

### 21.2 Those 3 records are noise, and that is itself a result

PR #30's only behavioural change is the window QUICKREF uses. If the QUICKREF decisions are
identical at every locus — they are — then nothing the fix touches can have changed, so the
3 differing records cannot come from the fix. They are run-to-run variation: a `1|0` → `1|1`
at `chr6:160615264`, `SUP` `2,4` → `2,3` at a locus no-called in both, and `RB`
`1164,1163` → `1161,1163` with `STDEV` `4,7` → `3,7` at `chr19:8768202`. Consensus and
read-sampling jitter under `--threads 24`.

**This sets the noise floor for the whole benchmark series: ~3 loci in 50,000, 0.006%.**
Differences below that are not findings. It also means the raw VCFs cannot be diffed
directly — record order varies between runs, so a plain `diff` reported 93,308 changed lines
where a comparison keyed on `(chrom, pos)` finds 3.

### 21.3 Why the null result is not believable

`fix/quickref-coordinates` does two things, and the second is not small: `src/call.rs` now
passes `target.start - 15` / `target.end + 15` into the CIGAR check, where `main` passes the
bare interval. Any indel within 15 bp of a locus now marks it as varying, so QUICKREF must
fire **less**. Tested directly, with both musl binaries built from this tree and run on
`test_data/small-test-phased.bam` over 10,144 tiled chr7 loci:

| binary | QUICKREF fires |
|---|---|
| `main` | 253 |
| `fix/quickref-coordinates` | **140** |

A 45% drop. For the real run to show a 0% drop, not one of ~30 reads at any of 1,042 loci
could have a net length change anywhere in two 15 bp flanks — implausible for ONT data even
if TR flanks are cleaner than tiled windows. Scaled down tenfold it would still be ~45 loci,
not zero.

### 21.4 The check to run

The binary shipped in `run_C/` is byte-identical to a fresh build of the branch
(`0e26dafe248696b780657d4a59794a09`, confirmed by rebuilding and stripping), and `main`'s
musl build is `5e5f1b8d6515a5ed7b6028d9218dd9e4`. The VCF header records `argv[0]`, so
`##command=STRdust ./STRdust-pr30-musl` proves the *filename*, not the build.

```bash
md5sum ./STRdust-pr30-musl        # expect 0e26dafe248696b780657d4a59794a09
```

on the machine that ran C. If it differs, C ran the wrong build and must be re-run; the rest
of the series is unaffected, since A and B predate the staged binary entirely.

### 21.5 Resolved: the binary on the remote was correct, the run still was not

`md5sum ./STRdust-pr30-musl` on the remote returned `0e26dafe248696b780657d4a59794a09`, the
staged build. So the *file* was right at the time of checking. The output nevertheless
matches `main` at every locus, and a re-run (§22) with a different binary moves 1,384 loci,
so **C did not execute the patched code**. The likeliest account is that the file was
replaced after the run rather than before it. C is void either way; §22 supersedes it.

### 21.6 Superseded: what this section got wrong

Then #30 is not merely insufficient (§4.1) but **inert on this data**, which is a stronger
statement than §17.4 anticipated and would make §15.4 step 2 — the tolerant QUICKREF — the
only remaining QUICKREF work. In that case the local chr7 result above needs explaining
too, since the same binary pair behaves differently there; the likeliest reason would be
that adotto's TR flanks are much cleaner than tiled windows, which is testable by re-running
the chr7 tiles against a TR-only subset with real coverage.

---

## 22. Run C2 (`pr30_50k_again`): #30 has a large, positive effect

2026-09-11. **The numbers are solid; which binary produced them is not yet confirmed
(§22.4).**

### 22.1 The result

Same 50,000 loci, same gate passed (`targets.bed` = `50a099114bd2a67311e9b12d3313ac53`).

| run | scored | exact | GT | QUICKREF fires | QUICKREF precision |
|---|---|---|---|---|---|
| main sensitive | 46,418 | 52.2% | 68.0% | 1,042 (2.1%) | 98.8% |
| C (`pr30_50k`) | 46,418 | 52.2% | 68.0% | 1,042 (2.1%) | 98.8% |
| **C2 (`pr30_50k_again`)** | **46,566** | **53.4%** | **68.9%** | **2,278 (4.6%)** | **98.1%** |

**+1.2 points of exact concordance and +0.9 of genotype agreement**, from more than doubling
QUICKREF's coverage at essentially unchanged precision. §17.4 predicted `fires` would stay
near 2.1% and said a material rise would mean "#30 does more than expected"; it did.

### 22.2 Why firing goes *up*, which §21.3 got backwards

`calculate_all_length_diff_from_cigar` returns a **net** length difference. Widening the
window therefore cuts both ways: it catches indels just outside the interval (firing goes
down), but it also lets a flank indel **cancel** one inside the repeat (firing goes up).
Verified directly on the chr7 test data with both local builds: `main` fires at 253 loci,
the patched build at 140 — and **7 of those 140 fire only under the patch**, which is the
cancellation happening.

C2 against `main` shows the same bidirectional signature at scale: **1,310 loci fire only in
C2, 74 only in `main`.** That is what the patched code looks like. C, by contrast, is
set-identical to `main` — the signature of code that never ran.

Whether the cancellation is *desirable* is a separate question §22.5 raises.

### 22.3 What the extra coverage buys, measured

`misc/quickref_paired.py pr30_50k_again noqr_50k` follows all 2,278 firings into run B,
where the same loci went through full genotyping:

| on the 2,278 loci QUICKREF answered | loci | share |
|---|---|---|
| QUICKREF exact | 2,235 | **98.1%** |
| exact in B (full genotyping) | 1,343 | **59.0%** |
| called in B, wrong | 574 | 25.2% |
| no call in B | 361 | 15.8% |
| QUICKREF wrong as well | 34 | 1.5% |

**348 loci rescued from no-call and 574 corrected.** Full genotyping manages 59.0% on loci
QUICKREF gets 98.1% right — and only 70.1% even where it produces a call at all, against
90.3% on the narrower set C fires on (§16.2). The loci #30 newly admits are *harder* for the
aligner than the ones QUICKREF already had, which is exactly where a cheap correct answer is
worth most.

This is the first direct evidence for §14.3's claim that a better QUICKREF is the largest
available win, and it is only a third of the way to v0.16's 32.9% firing.

### 22.4 The open question: which binary ran

Three runs, and the labels do not add up:

- `main_50k` (Sep 9, `./STRdust`) → 1,042 fires.
- `pr30_50k` (`./STRdust-pr30-musl`, md5 confirmed `0e26dafe…` on the remote) → 1,042,
  set-identical to `main`.
- `pr30_50k_again` (`./STRdust`, freshly uploaded) → 2,278.

A freshly built `main` cannot differ from the Sep 9 `main` — the source has not moved since
`e067e0f` (Aug 29). So C2's binary is almost certainly the patched build, and C's was not,
despite the filename and the later md5. **One command settles it**, on the remote:

```bash
md5sum ./STRdust        # 0e26dafe248696b780657d4a59794a09 = patched, 5e5f1b8d6515a5ed7b6028d9218dd9e4 = main
```

If it is `5e5f1b8d…`, then C2 ran `main` and the Sep 9 `./STRdust` was *not* built from
`main` — which would put run A's provenance in question and matter far more than #30.

**Contributing cause, for next time:** the repo was left on `main` after the staged binary
was built (§17.5), so `make musl` in the working tree rebuilds `main`, not the PR. Stage
binaries under content-addressed names and record the md5 in the run's own notes before
running, rather than relying on the filename.

### 22.5 What to do once identity is confirmed

1. **Merge #30.** It is a correct fix independent of the benchmark (§4.2), and now has
   +1.2 points behind it.
2. **Check the cancellation case before trusting the mechanism.** Firing rose partly because
   flank indels cancel repeat indels. Precision only fell 98.8% → 98.1%, so it is not doing
   much harm at this scale, but "net difference over a padded window" is a weaker test than
   "no indel anywhere in the padded window". The 34 wrong firings are the place to look, and
   `|diff| <= 3` plus a minimum read count (§15.4 step 2) would replace the cancellation
   accident with a deliberate tolerance.
3. **Re-run the C2 command on the merged `main`** so the series has a clean post-#30
   baseline before the tolerant-QUICKREF work starts.

---

## 23. Resolved: the labels were wrong, and #30 is a net regression

2026-09-11. **This supersedes §21 and §22, both of which drew the wrong conclusion from the
wrong label. Nothing is void; one run was misattributed.**

### 23.1 The binaries, settled by content

| binary | md5 | QUICKREF on chr7 tiles | QUICKREF on the 50k |
|---|---|---|---|
| `main` (`e067e0f`) | `5e5f1b8d6515a5ed7b6028d9218dd9e4` | 253 | **2,278** |
| `fix/quickref-coordinates` (#30) | `0e26dafe248696b780657d4a59794a09` | 140 | **1,042** |

`main` fires **more**, the patched build fires **fewer** — consistently, locally and on the
real data, at a similar ratio (1.8× and 2.2×). With that direction fixed, every run falls
into place:

| run | binary named | fires | what it actually was |
|---|---|---|---|
| A (`main_50k`) | `./STRdust` | 1,042 | **the #30 build**, not `main` |
| C (`pr30_50k`) | `./STRdust-pr30-musl` (md5 confirmed) | 1,042 | the #30 build ✓ |
| C2 (`pr30_50k_again`) | `./STRdust` (md5 `5e5f1b8d…`) | 2,278 | **`main`** — the first true `main` run |

A and C being set-identical was never a contradiction: they ran the same code. The
contradiction was the label on A.

### 23.2 #30's real effect, and it is negative

The clean head-to-head is C (#30) against C2 (`main`), same 50,000 loci, both md5-confirmed:

| run | scored | exact | GT | fires | QUICKREF precision |
|---|---|---|---|---|---|
| C2 = `main` | 46,566 | **53.4%** | **68.9%** | 2,278 (4.6%) | 98.1% |
| C = #30 | 46,418 | 52.2% | 68.0% | 1,042 (2.1%) | 98.8% |

**#30 costs 1.2 points of exact concordance and 0.9 of genotype agreement.** What it does,
locus by locus:

- **removes 1,310 firings — of which 1,279 (97.6%) were correct**
- adds 74 firings, all 74 correct
- the 1,310 loci it sends back to full genotyping are answered exactly **48.8%** of the
  time there, with 162 becoming no-calls

So it discards ~1,279 correct cheap answers to suppress ~31 wrong ones, and the alignment
path gets fewer than half of them right afterwards. The precision gain is real but tiny
(98.1% → 98.8%) and bought far too expensively.

### 23.3 Why — and what to do instead

The ±15 padding is not wrong in principle; combined with a **strict** `diff != 0` it is
simply too blunt. At a TR locus there is nearly always some indel noise within 15 bp, so
padding plus zero tolerance disqualifies almost everything. The fix is not to drop the
padding but to pair it with the tolerance the non-batched path already has:
`|diff| <= 3` and a minimum read count (§15.4 step 2, §14.3). v0.16 reached 32.9% firing at
96.8% precision with exactly that combination; `main` reaches 4.6% and #30 drops it to 2.1%.

**Do not merge #30 as it stands.** It is still a correct coordinate fix (§4.2) and the
boundary change in `parse_bam.rs` should survive, but the `QUICKREF_PADDING` in `call.rs`
must land together with a tolerance, not before it. Worth splitting the PR to measure the
two halves separately.

### 23.4 Corrections to earlier sections

- **§14.2's "main fires 1,042 (2.1%)" is wrong** — that was the #30 build. `main` fires
  2,278 (4.6%). The comparison against v0.16's 32.9% still holds qualitatively (7× rather
  than 16×), and the §14.2 explanation of *why* `main` looks worse than v0.16 is unaffected.
- **§14.1 (`fast` beats `sensitive`) is unaffected**: both A modes ran the same binary, so
  the comparison is internally valid. The absolute numbers are the #30 build's, not `main`'s.
- **§16 (run B) is unaffected**: `--alignment-all` disables QUICKREF, and #30 touches nothing
  else, so B measures the same genotyper either way.
- **§16.2's paired analysis stands**, relabelled: it is the #30 build's 1,042 firings.
- **§21 is wrong throughout** — C was not void and the binary was not swapped.
- **§22 is wrong in its central claim**: C2's +1.2 points is `main` beating #30, not #30
  beating `main`. §22.3's measurement is right but describes `main`'s QUICKREF, not #30's.

### 23.5 The process fix, implemented

The root cause was that a run was identifiable only by a filename. `misc/benchmark_truth.py`
now **fingerprints the binary it executes** — md5, size and resolved path — prints it at
startup and records `binary_md5` as a column in `<prefix>.summary.tsv`. A series that
compares builds is only interpretable if every run states which bytes it ran; the VCF
header's `argv[0]` proves the filename and nothing more.

Verified against the two builds above: it reports `5e5f1b8d…` and `0e26dafe…` correctly, and
degrades to a warning plus a blank column when the binary cannot be read.

**Also:** stop leaving the working tree on `main` after building a staged binary (§17.5 did),
because `make musl` then rebuilds `main` under whatever name is expected. Every future run in
this series should have its `binary_md5` checked against the intended build *before* its
numbers are read.

---

## 24. Run D is staged: QUICKREF padding × tolerance, one binary, four arms

> **Executed; see §28 and §29.4.** The tolerant-QUICKREF arm it was built around was measured, found to buy coverage by flattening 1-3 bp alleles, and the knobs were retired.

Staged 2026-09-11 on branch `feat/quickref-tolerance` (off `main`, CI green: fmt, clippy,
87 unit + 5 integration tests). This replaces §20's separate flank-folding matrix as the next
run; §20 still stands for the `minlen` question.

### 24.1 What changed in the code

The batched QUICKREF criterion is now three CLI knobs instead of a constant, so the whole
matrix runs from **one binary** and no arm needs a rebuild — which is what made §23's
mix-up possible in the first place.

- `--quickref-padding` (default 0) — how far outside the interval a read's indels count.
- `--quickref-tolerance` (default 0) — net length difference still counted as
  reference-like. **Padding and tolerance belong together**: widening the window with a
  strict `diff != 0` rejects nearly every read, which is exactly why #30 lost 1.2 points
  (§23.2).
- `--quickref-min-reads` (default 0) — reads that must have been inspected before a locus
  may be called reference without aligning. `TargetInfo` now counts them.

`src/parse_bam.rs` carries #30's boundary fix unchanged, and is always on. Defaults keep
everything else as `main` has it, so the only default-behaviour change is that boundary fix.

### 24.2 Validated against the builds it must reproduce

On the chr7 test data, 10,144 tiled loci, against the two known binaries:

| configuration | QUICKREF fires |
|---|---|
| `main` (`5e5f1b8d…`) | 253 |
| knobs, defaults — boundary fix only | 237 |
| knobs, `--quickref-padding 15` | **140** |
| #30 (`0e26dafe…`) | **140** |
| knobs, padding 15 + tolerance 3 + min-reads 5 | **2,604** |

`--quickref-padding 15` reproduces #30 **exactly**, which is the check that the knobs are
faithful rather than merely similar. The boundary fix alone costs 16 firings. The tolerant
configuration reaches 25.7% of loci here, against v0.16's 32.9% on the real benchmark —
the first sign that §14.3's "largest available win" is reachable from the batched path.

Local firing counts carry no precision information; that is what the 50k run is for.

### 24.3 The md5 gate, so §23 cannot recur

`benchmark_truth.py` gained `--expect-binary-md5`. It refuses to start, with exit status 1,
unless the binary's content is exactly the intended build:

```
binary md5 mismatch, refusing to run:
  expected 5e5f1b8d6515a5ed7b6028d9218dd9e4
  found    0e21900e9c92ad9bd4464be51dfe0d03  (…/knobs.bin)
```

Together with the `binary_md5` column already in `summary.tsv` (§23.5), a run now states
which bytes produced it *and* cannot silently be the wrong ones. Tested in both directions.

### 24.4 Staged for transfer

`~/Downloads/fix_strdust/run_D/` — copy the directory next to the data and run `RUN.sh`:

Built from `a435f9a` with the tree clean, so the binary maps to a commit rather than to a
working state — the provenance failure of §23 in miniature. `make musl` on
`feat/quickref-tolerance` reproduces this md5 exactly.

| file | md5 |
|---|---|
| `STRdust-quickref-knobs` | `fda021812c9a72128d9eaebe63ef9a8a` |
| `benchmark_truth.py` | `fe8be95a41f00597010d8813181c281c` |
| `compare_runs.py` | `43773176cd5f0c9d0fcc523d691ca86b` |
| `quickref_paired.py` | `275b1672baffc1e985d57f5007d23f8d` |
| `RUN.sh` | `51ed288c6ae4a4835aa79557c15697f2` |

The arms, all on the same 50,000 loci and `--seed 1`, all gated on the binary md5:

| arm | flags | what it isolates |
|---|---|---|
| D1 `d1_bound_50k` | defaults | the boundary fix alone |
| D2 `d2_pr30_50k` | `--quickref-padding 15` | #30 as it stands |
| D3 `d3_tol3_50k` | padding 15, tolerance 3, min-reads 5 | the rework |
| D4 `d4_tol1_50k` | padding 15, tolerance 1, min-reads 5 | a tighter tolerance |

`RUN.sh` finishes by printing the five-way comparison against C2 (`pr30_50k_again`, the
`main` baseline) and the paired rescue/correction breakdown of each arm against run B.

### 24.5 What to read, and the decision rules

1. **D2 must reproduce #30.** On the 50k sample that is literal: 1,042 fires. On the local
   10,000-locus sample (§24.7) it cannot be, because the loci differ — there the check is
   that **D2 against the local `main` arm reproduces C against C2 in direction and rough
   magnitude**: firing roughly halving, exact concordance down about a point. The knobs
   already reproduce #30 exactly on the chr7 test data (§24.2, 140 both ways), so a local
   mismatch would point at the sample, not the code. Either way this is the first number to
   check: if D2 does not behave like #30, nothing else in D is interpretable.
2. **D1 against C2** prices the coordinate fix on its own. Expect a small loss (it removed
   16 of 253 locally); if it costs materially more than that, the boundary fix is not as
   free as §23.3 assumed and should be reconsidered rather than merged by default.
3. **D3 is the candidate.** It wins if it beats C2's 53.4% exact *and* holds QUICKREF
   precision near 98%. v0.16 managed 32.9% firing at 96.8%; anything approaching that at
   98% would be worth more than every threshold change in §13 combined.
4. **Watch precision, not just coverage.** The failure mode is a tolerance that fires on
   genuinely variant loci. `quickref_paired.py` splits each arm into rescues and
   corrections; an arm that raises firing while its "QUICKREF wrong" count grows faster than
   its rescues is buying coverage with errors.
5. **D4 against D3** says whether tolerance 3 is too loose. If D4 is close to D3 on exact
   but cleaner on precision, prefer it.

Nothing here is merged. `feat/quickref-tolerance` is local only, and #30 stays unmerged
pending D (§23.3).

### 24.6 Replication across seeds

`--sample-loci N --seed S` draws a random subset of the catalog, so two seeds are two
independent draws of loci from the same data. Every conclusion from run D that leads to a
default change should be confirmed on a second seed before the change is made.

This is cheap, and it is the only control available for the thing that has bitten this
investigation repeatedly: a number that looks decisive on one sample. Run D's arms differ
from each other by fractions of a point in places, and §21 already established a run-to-run
noise floor of ~0.006% of loci from threading alone — sampling noise across 10,000 loci is
considerably larger than that.

Rule: **a difference that does not reproduce at a second seed is not a result.** Two seeds
agreeing on direction and roughly on magnitude is enough; they will not agree exactly, and
should not be expected to.

Outputs are tagged `n<loci>s<seed>` so batches never collide.

### 24.7 Running it locally

The benchmark data now lives on the laptop (`~/testdata`), which takes the server and the
file-shuttling out of the loop. 12 cores against the server's 24, so a 50,000-locus arm is
roughly an hour; the matrix runs at **10,000 loci** instead, because every question it asks
is a difference *between arms* and the arms only need to match each other.

`~/testdata/bench/` holds both binaries, the three scripts, and:

- `preflight.sh` — inputs, indexes, both binary md5s, and an actual read fetch from the CRAM.
- `run_matrix.sh [n_loci] [seed] [threads]` — the five arms, each md5-gated, skipping arms
  already finished so an interrupted run resumes. Ends by printing the comparison and the
  `binary_md5` each arm recorded.
- `go.sh` — index (atomically, via a temp name) → preflight → 300-locus smoke → the matrix.
  `go.sh replicate` runs the same matrix at seed 2 per §24.6.

The local baseline arm is `main` itself (`5e5f1b8d…`) rather than C2, so the matrix is
self-contained and does not depend on the 50k series or its provenance.

---

## 25. The truth set attributed whole spanning deletions to tiny repeats

Found 2026-09-11 while characterising the length-error distribution of run D.

### 25.1 What was wrong

`read_truth_over_repeats` attributes a variant to a repeat when the two *overlap*, then sums
the variant's whole per-haplotype length change into the locus. A deletion reaching far past
the repeat therefore donated its entire length:

```
chr22:22893684  ref_len=61   truth allele = -145777   (2390x the locus length)
chr15:84180336  ref_len=116  truth allele = -78923    (680x)
chr3:146669771  ref_len=46   truth allele = -4902     (107x)
```

A 61 bp repeat cannot contract by 145 kb. No genotyper can reproduce these, and STRdust
mostly **no-called** them — the correct response, scored as failure.

### 25.2 Small overall, decisive in one stratum

19 alleles across 15 of 10,000 loci — 0.15%, invisible in the headline numbers. But the
stratum is assigned by truth allele size, so every artefact lands in the largest one:

| stratum | loci with an impossible truth allele |
|---|---|
| 1-10bp, 11-50bp, 51-200bp, reference | 0% |
| **>200bp** | **15 of 38 — 39.5%** |

Dropping them from `main`'s `>200bp` stratum: 38 loci → 23, **14 no-calls → 1**, exact
41.7% → 45.5%, ≤5 bp 83.3% → 90.9%. Nearly every no-call in the stratum was a locus whose
truth was nonsense.

**This invalidates §15.4 step 4's reading**, which took the `>200bp` stratum as evidence that
"neither path is convincing at long expansions" and proposed a run enriched for long
expansions to decide `fast` versus `sensitive`. That stratum was ~40% artefact and its
recall was an artefact almost entirely. The question is still open; the evidence offered for
it was not evidence.

### 25.3 The fix

`clamp_to_locus` restricts a deletion's contribution to the bases it actually removes from
the locus — the overlap of the deleted interval with the repeat. Insertions are not clamped:
the inserted bases enter at a position inside the locus, so their whole length belongs to it.
Deleted bases start one base after the VCF position, per the anchor-base convention.

Verified on ten boundary cases (spanning, half-covering, starting inside and running past,
ending exactly at the locus start, one base in, fully contained, insertion, SNV) and end to
end: a 5 kb deletion spanning a 60 bp locus goes from `[-5000, 0]` to `[-60, 0]` while a
contained 10 bp deletion and a 30 bp insertion at neighbouring loci are untouched. On inputs
with no spanning deletion the truth set is byte-identical.

`TRUTH_CACHE_VERSION` is bumped to 2, so every existing cache rebuilds rather than silently
serving pre-fix truth.

### 25.4 What this invalidates, and what it does not

The truth set has changed, so **runs scored before this fix are not comparable with runs
scored after it** on any stratified number. Within-run and within-matrix comparisons made
before the fix remain valid among themselves — every arm of run D seed 1 and seed 2 was
scored against the same truth.

The fix was deliberately *not* applied while seed 2 was running: `run_matrix.sh` starts a
fresh Python process per arm, so editing the staged script mid-matrix would have scored
later arms against different truth than earlier ones — the §23 failure in a new costume.
The staged copy in `~/testdata/bench/` is updated only once seed 2 finishes, and both seeds
are then re-run.

A truth set built from a general call set will always carry some artefacts of this kind;
the aim is not perfection but that a large variant overlapping a repeat should not be
reported as that repeat's genotype.

---

## 26. POA scoring, and the tuning flags' exit plan

### 26.1 The consensus scoring was never measured

`src/consensus.rs` built every consensus with `Scoring::new(-12, -6, match 3 / mismatch -4)`
under a comment saying the values were "empirically determined ... further testing on other
repeats would be good", and — pointedly — that the goal was to "make sure the consensus does
not get longer than the individual insertions".

It does get longer, on a quarter of homozygous-reference alleles, concentrated at +1 and +2
(§16.3, §20.2). The gap penalty is the obvious suspect: at -12 per gap base against +3 per
match, a single read's insertion can open a node of its own and be carried into the
consensus path. rust-bio's POA has no affine gaps and ignores `gap_extend` entirely
(rust-bio#677), so the model is linear and three numbers are the whole of it.

`--poa-gap-open`, `--poa-match` and `--poa-mismatch` now expose them, defaulting to the
current values so nothing changes by default. Penalties are given as positive numbers on the
command line and negated internally.

This is the knob that acts directly on the dominant error, unlike the QUICKREF knobs, which
move coverage and speed but not consensus accuracy.

### 26.2 These flags are temporary, and are hidden

> **Outcome, for the record (§29.4).** All six were resolved differently from the guesses
> below: `--poa-gap-open`, `--poa-match` and `--poa-mismatch` were **retired** rather than
> retuned, because they contribute −0.1 once the consensus endpoint is fixed; the three
> `--quickref-*` knobs were **retired** because no setting of them is one a user would
> choose; `--junction-window` is **kept hidden**, being a no-op on the default path (15 loci
> in 50,000); and `--poa-trim-fraction` and `--poa-medoid-seed` stayed hidden but had their
> **defaults changed** to the measured optimum. The user-facing knob that emerged was one
> nobody had predicted: `--priority` (§29.6).

All six tuning flags — three QUICKREF, three POA — are `hide = true`: absent from `--help`,
absent from the README, fully functional. They exist so a sweep needs no rebuild, which is
what lets runs be compared by binary content rather than by which branch someone had checked
out (§23). They are **not a stable interface.**

The exit plan, once the measurements settle:

| flag | likely disposition |
|---|---|
| `--quickref-tolerance` | promote: a genuine accuracy/resolution trade a user might want |
| `--quickref-padding` | fold into the tolerance decision, then remove |
| `--quickref-min-reads` | fold in, then remove |
| `--poa-gap-open` / `--poa-match` / `--poa-mismatch` | bake the winning values in and remove |
| `--junction-window` | bake in and remove, unless it turns out to be locus-dependent |

A flag that survives should do so because someone would reasonably set it, not because it
was convenient during tuning. Anything still hidden when this work closes should be deleted,
and the comment block above the group in `src/main.rs` says so at the definition site.

### 26.3 The POA arms to run

Cheap: no rebuild, and they only touch the sensitive path.

| arm | flags | hypothesis |
|---|---|---|
| P1 | `--poa-gap-open 20` | a dearer gap collapses stray insertions; +1/+2 mass should fall |
| P2 | `--poa-gap-open 30` | as above, further; watch for real alleles being collapsed too |
| P3 | `--poa-gap-open 8` | the opposite direction, to confirm the mechanism is the gap at all |
| P4 | `--poa-match 5` | changes the match/gap ratio without touching the gap |
| J1 | `--junction-window 10` | folds in less stray flank sequence (§20, #24) |
| J2 | `--junction-window 0` | folds in none at all |

`--junction-window` replaces the "edit `genotype.rs:947`, rebuild, re-run" recipe §20 called
for, so the whole consensus matrix now runs from a single binary. Default 30, hidden, same
exit plan as the rest.

Read the error spectrum (§20.2) first: if the +1/+2 mass moves with the gap penalty, the
consensus over-call is a scoring artefact and the fix is a better default. If it does not
move, POA scoring is exonerated and the remaining suspects are the junction fold window
(§20) and the consensus construction itself, which `--mode fast` sidesteps entirely by
reading length off the alignment (86.5% vs 71.9%, §14.1).

**Detection is not at risk in any POA arm** — POA runs only after a locus has been routed to
full genotyping, so it changes reported lengths, never whether an expansion is found (§25's
detection table is invariant across every arm run so far).

---

## 27. Independent parameter audit, 2026-09-11

An independent read-only pass over the codebase looking for anything tunable that could move
a called allele, deliberately not steered toward the knobs already exposed. It found more
than expected, including one confirmed source of non-reproducibility. Fixes in `cd5d7c9`.

### 27.1 Read downsampling was unseeded — confirmed, and quantified

`downsample_reads_inplace` (`src/parse_bam.rs`) used `rand::rng()`, the OS-seeded thread
RNG, while `src/consensus.rs` has seeded its own downsampling with a constant since it was
written, under a comment giving the reason: the subset kept is a performance measure, not a
genotyping decision. Two identical invocations could therefore report different allele
lengths at any locus above the read cap.

**Confirmed against the benchmark, not just argued:** of the three loci that differed between
two runs of identical code over 50,000 loci (§21.2), **two sit exactly at the per-haplotype
cap** — `SUP=20,30` and `SUP=30,29`. The cap truncates to exactly 30, so `SUP=30` is its
signature.

**But it is rarer than it sounds on this data.** Per-haplotype support here: median 10,
p99.9 = 23, maximum 30. Only **2 of 19,070 haplotype observations** reach the cap and none
reach the 60-read unphased cap, so the injected noise is ~2 loci per 10,000 — about 0.02%,
two orders of magnitude below the 0.3-1.1 point seed-to-seed sampling noise (§24.6). It was
therefore **not** worth re-baselining the completed matrices over; the fix is carried into
every subsequent run instead. On deeper data (the 90 GB CRAM, for instance) it scales with
how often the cap is hit and matters much more.

**Still unexplained:** the third differing locus had `SUP=2,4` against `2,3` — low coverage,
so downsampling cannot be the cause, and the *support count itself* changed. A second
non-determinism source exists, rarer than this one, and is not yet identified.

### 27.2 Two off-by-one bugs, both real, both inert here

- **`Batch::new` took the batch end from `repeats.last()`** (`src/batching.rs`) though repeats
  sort by *start*, so an interval nested inside an earlier, longer one is last while ending
  first. `create_batches` already tracked the correct maximum and discarded it. The fetch
  region could stop inside a long repeat and starve it of the reads covering its tail.
  **Zero nested pairs exist in the adotto catalog** (0 of 4,952 batches), so nothing measured
  here changes; it would bite a user whose BED has overlapping or nested regions.
- **`VCFRecord::single_read` still used `end - start`** (`src/vcf.rs`), the last survivor of
  #22, making single-read calls one base too long, and underflowing on a contracted allele.
  Reachable only with `--unphased` and exactly one read.

Both verified inert on this data rather than assumed: 400 loci produce a byte-identical VCF
with and without them, so the seeding change remains solely attributable.

### 27.3 The findings worth acting on next

Ranked by the audit, with the mechanism that makes each matter:

1. **`parse_cs` sums insertions and never subtracts deletions** (`src/genotype.rs`). Every
   insertion inside the junction window is *concatenated*; the `'-'` branch only advances
   `ref_pos`. So **`--junction-window` is one-sided**: widening can only lengthen alleles,
   narrowing can only shorten them. This is a structural explanation for the 25.9%-too-long
   against 2.3%-too-short asymmetry, and for why `--mode fast` is more accurate — it derives
   the allele from a reference span, so deletions shorten it automatically. **Prediction: the
   junction-window sweep should move length monotonically. If it shows an optimum, a second
   mechanism is cancelling the first, and that is worth more than the chosen value.**
2. **`remove_outliers` trims at ±2 std dev before the POA** (`src/consensus.rs`), with a
   `std_dev < 5` short circuit that disables it entirely at clean loci. Symmetric trimming on
   an asymmetric distribution: at an expanded locus it cuts the longest reads, biasing short.
   It is the only length filter on the always-live path, and it decides which reads the POA
   ever sees — so it should be swept **before** the POA scoring, to separate "the POA invents
   bases" from "the POA is fed a truncated distribution".
3. **`consensus()` applies `--support` *after* outlier removal**, so a 4-read haplotype that
   loses two reads to trimming becomes a no-call despite having had support. A product of two
   parameters rather than either alone, and a concrete contributor to the ~6% no-call rate.
4. **`is_primary`-only in `find_insertions` plus minimap2's `map_ont` defaults**
   (`zdrop=400`, `max_gap=5000`) — a large enough expansion can be split across alignments
   with the non-primary piece discarded, **silently undersizing** rather than no-calling.
   With `flanking = 5000` capping sizable expansions outright, this is the only item on the
   detection axis that no queued sweep touches. See §27.4.
5. **The POA is seeded by `seqs_bytes[0]`** — the first sampled read's indels become the
   backbone the others align onto. Worth testing whether the +1/+2 excess is seed-read
   dependent by sweeping `DOWNSAMPLE_SEED` alone, which is now meaningful because the
   downsampling is deterministic.

Unphased-only, so not on the current benchmark's path but high-consequence if `--unphased` is
used: `find_roots`' dissimilarity threshold of `5.0` (`src/phase_insertions.rs`) means a
genuine minority expansion supported by few reads is **discarded and the reference allele
split in two to replace it**, producing a confident `0/0`. `EXPANSION_OUTLIER` is the only
signal that this fired, and it is INFO-only.

### 27.4 The experiment the benchmark cannot do

A **synthetic expansion ladder**: reads carrying 100/500/1000/2000/4000/8000 bp insertions at
a single locus, genotyped with defaults. It answers the question no random sample of HG002
can — *at what size does sizing break, and does it fail loudly (no-call) or quietly
(undersized)?* Undersizing is the dangerous answer and is what item 4 predicts. Everything
else in this document moves averages; this one decides whether a pathogenic expansion is
reported at all.

---

## 28. The consensus sweep: the over-call is stray junction insertions

Run 2026-09-11/12, 10,000 loci, seed 1, corrected truth, one md5-gated binary
(`9015030012cd0b1aedd3608ae4331b67` from `ba05a42`), nine arms, no rebuilds.

### 28.1 Results

All arms below have **identical no-call counts (574)** unless stated, so these are not bought
by declining to answer.

| arm | exact | vs base | reference | 1-10bp | 11-50bp | verdict |
|---|---|---|---|---|---|---|
| `--minlen 0` | 34.2% | **-18.8** | 36.6% | 21.5% | 29.6% | rejected |
| `--poa-gap-open 8` | 51.5% | -1.5 | 54.8% | 34.4% | 45.7% | control, worse as predicted |
| `--poa-match 5` | 52.5% | -0.5 | 55.8% | 35.3% | 46.9% | control, worse as predicted |
| **base** (defaults) | 53.0% | — | 56.4% | 35.4% | 46.1% | current |
| `--poa-gap-open 20` | 56.0% | +3.0 | 59.6% | 37.7% | 47.3% | good |
| `--poa-gap-open 30` | 58.3% | +5.3 | 62.2% | 38.8% | 48.4% | good |
| `--junction-window 10` | 62.1% | **+9.1** | 66.1% | 41.3% | 53.9% | good |
| **`--minlen 5`** | **64.5%** | **+11.5** | 68.7% | 42.9% | 56.9% | **best** |
| `--junction-window 0` | *82.5%* | — | *87.1%* | *52.2%* | *65.2%* | **artefact, see 28.3** |

Every improving arm improves **every** stratum, including `1-10bp` — the small-allele
resolution axis. This is the opposite of the QUICKREF tolerance knob (§24), which bought the
aggregate by giving that axis away.

### 28.2 The mechanism, and that the knobs are not interchangeable

Three knobs improve, by two different routes, visible in the error spectrum at truth
homozygous-reference loci:

| arm | +1 bin | +2 and beyond | how |
|---|---|---|---|
| base | 8.0% | 8.4% / 3.7% / 2.5% | — |
| `--poa-gap-open 30` | **7.7%** | 7.3% / 3.1% / 1.9% | the consensus refuses to absorb insertions |
| `--minlen 5` | 8.3% | **4.2%** / 1.8% / 1.3% | stray insertions never reach it, filtered by size |
| `--junction-window 10` | 8.1% | **5.4%** / 2.2% / 1.6% | same, filtered by distance from the junction |

So the dominant error is **small stray insertions near the junction being folded into the
allele** — `parse_cs` concatenates every one inside the window (§27.3) — and it can be
attacked either by keeping them out or by making the POA reluctant to take them.

`minlen` and `junction-window` have near-identical spectra, so they are probably filtering
much the same insertions and may not compound with each other. The gap penalty works
differently and is the only knob that touches the **+1 bin at all** — which remains 30% of
all over-calls and is, after this sweep, the largest unexplained piece.

`--poa-match 5` is a weaker version of lowering the gap penalty (it makes a gap relatively
cheaper), confirming the effect is the gap-to-match **ratio**. It earns no place as a
separate flag.

### 28.3 `--junction-window 0` is survivorship, not a result

It scores 82.5% exact on **2,138 of 10,000 loci**, having no-called 7,862 including 6,391
reference loci. With a zero-width window only an insertion landing exactly on the junction
counts, so almost everything is rejected and the survivors are the easy cases.

This is the turnover predicted for the window: aligners genuinely misplace a repeat's
insertion by a few bases, so zero tolerance discards real signal. **It is also the clearest
argument for keeping no-call counts beside every concordance number** — the headline alone
reads as the best arm in the sweep by 18 points.

### 28.4 `--minlen` 5 -> 1 was a regression, now measured

The dose-response is clean and monotonic: `--minlen` 0 / 1 / 5 gives 34.2% / 53.0% / 64.5%
exact, improving every stratum. §19.2 found the default was changed from 5 to 1 in `02d9540`
(2025-11-14), a one-line change inside a commit titled "add min_haplotype_fraction parameter
to CLI and related functions", with no test and no mention in the message. **That change cost
about 11 points of exact concordance.**

It also settles §19's semantic question in the direction §19.4 guessed but could not show:
the strict `>` is load-bearing. Making `--minlen 1` mean "1 and longer" would land on the
`--minlen 0` column, i.e. -18.8 points. The fix is `>=` **with the default set to 2 or
higher**, never `>=` at the current default.

### 28.5 Not yet known

- **Where the optima are.** Every improving knob was still improving at the edge of its
  tested range. `--poa-gap-open` 45 and 70 are running; `--minlen` 3 and 10 are not yet run.
  A knob that only ever improves across its tested range is a knob whose range is too small.
- **Whether they compound.** `--minlen 5` + `--poa-gap-open 30` act on different parts of the
  spectrum and should add; `--minlen` + `--junction-window` probably overlap. Untested.
- **Seed 2.** Nothing here is replicated. The effects are 3-11 points against seed-to-seed
  noise of 0.3-1.1 (§24.6), so they should survive, but no default moves until they do.
- **Recall falls as precision rises** across every improving arm (65.1% -> 62.0% for
  `--minlen 5`). Worth understanding before choosing a value, not just tallying exact matches.
- **The +1 bin.** Untouched by everything except the gap penalty, and then barely.

---

## 29. The +1 bin was a bug in rust-bio's POA, not a parameter

2026-09-12/13. This section supersedes the tuning conclusions of §28: most of what those
knobs appeared to buy was a partial workaround for a single upstream defect.

### 29.1 The diagnosis

`bio::alignment::poa::Aligner::consensus()` chooses its endpoint with `max_by_key` over a
cumulative score to which **every edge contributes a weight of at least 1**. The score
therefore increases strictly along every edge, so the argmax is **provably always a sink**:
the consensus runs to the deepest point in the graph rather than the best supported one.
A global alignment forces a read longer than every graph path to emit its surplus as
insertion nodes, and surplus at either terminus creates a deeper sink that wins regardless of
how little supports it — one read in fifteen, with no vote taken.

Three independent lines of evidence, none of which any other hypothesis predicts:

| evidence | value |
|---|---|
| consensus longer than the median read, vs shorter | **29 : 1** |
| P(too long \| median read is exactly reference), 0-4 reads | 14.9% |
| … 10-14 reads | 26.9% |
| … 15-19 reads | **32.6%** |

**The defect gets worse with more reads.** A consensus cannot do that; an estimator tracking
the maximum can. And `MRL`, the median read length of the same cluster, is 91.9% exact where
the consensus is 71.2% — and is **invariant across every scoring and filtering arm tried**,
because the reads were never wrong.

It also explains the signature §28 could not: surplus in the *interior* forms a weight-1
bubble that the heaviest-edge choice rejects — that is the +2/+3/+4 mass, which raising
`--poa-gap-open` did drain — while surplus at a *terminus* wins on depth, which no score can
touch. Hence +1 sat at ~8% through every scoring change.

### 29.2 The fix, and why no fork was needed

`Aligner::graph()` is public and `POAGraph` is a plain petgraph graph, so the corrected
traversal lives in `src/consensus.rs` (`trimmed_consensus`): the same heaviest-edge walk,
with the ends corrected — a terminal base is dropped unless the edge reaching it carries
`--poa-trim-fraction` of the cluster's reads. No fork (a git dependency cannot be published
to crates.io) and no vendored module (~500 lines frozen at 4.0.0 to fix ~30).

**Worth reporting upstream**: the endpoint selection is wrong for any noisy-read use.

### 29.3 What it is worth

| | seed 1 (tuning) | seed 3 | seed 4 | seed 5 | seed 6 |
|---|---|---|---|---|---|
| base | 53.0% | 53.4% | 53.7% | 52.8% | 53.6% |
| trim 0.35 | 79.4% | 79.4% | 79.7% | — | — |
| best combination | 80.8% | — | — | 80.8% | 81.4% |
| **`fast` + trim 0.35** | — | — | — | **84.6%** | **85.0%** |

**+26.0 points on seeds never used for tuning**, against +26.4 where it was tuned — no
measurable selection effect despite comparing well over thirty arms. 0.35 is the optimum of
a six-point ladder (0.10/0.20/0.35/0.40/0.45/0.50) and turns over on both sides; above it
the −1 bin grows, which is real sequence being cut rather than one-read overhangs.

### 29.4 The tuning knobs were measuring the same defect

Once the endpoint is corrected, the knobs of §28 collapse:

| on top of trim | gain |
|---|---|
| `--poa-gap-open 30` | **−0.1** |
| `--minlen 3` | +1.4 |
| `--poa-medoid-seed` | +0.3 |
| `--junction-window 10` | +0.5 |

`--poa-gap-open`'s standalone +5.3 (and +11.4 at 70) was almost entirely an indirect
mitigation of the sink rule. It, `--poa-match` and `--poa-mismatch` were retired rather than
given new defaults.

### 29.5 `--mode fast` wins everywhere, including where it was expected to lose

On 5,072 loci that actually carry an expansion over 200 bp:

| arm | exact | ≤5bp | no-calls |
|---|---|---|---|
| sensitive | 46.1% | 88.3% | 594 |
| sensitive + trim + minlen | 69.1% | 91.6% | 606 |
| fast | 63.1% | 91.9% | 498 |
| **fast + trim** | **71.9%** | **92.6%** | 498 |

Detection on that set: 92.2% found for `fast` against 91.7% for sensitive, with `short`
(silently undersized) at 1.0% for both. So ~7% of real expansions are missed by either mode,
almost always as a visible no-call. `fast` also builds its ALT through the same POA, which is
why trimming improves it by 8 points too — its advantage comes from cleaner *per-read*
alleles, not from avoiding the consensus.

`--mode` therefore defaults to `fast` (breaking), and `--minlen` is a no-op on that path.

### 29.6 `--priority`, replacing a threshold nobody chose

The recall cost attributed to trimming was not trimming: of 98 loci that moved from variant
to reference, **94 still had a called length differing from the reference** — `is_similar_to_ref`
was rounding them away. Median REF length 77, so `floor(77/20)` with a strict `<` tolerated
2 edits. The length was correct in `RB` all along.

| `--priority` | recall | precision | F1 |
|---|---|---|---|
| `sensitive` | 99.0% | 54.2% | 70.1 |
| `balanced` (new default) | 85.5% | 58.7% | 69.6 |
| `precise` | 39.9% | 91.8% | 55.7 |
| *old default* | *52.1%* | *83.0%* | *64.0* |

The old rule was on no frontier. **The setting affects only the emitted `GT`; `RB`, `FRB` and
`MRL` are byte-identical under all three**, verified.

§13.4's "no threshold rule fixes the over-calling" was measured against a length distribution
distorted by the sink bug and does not carry forward.

### 29.7 Two bugs found by trying to benchmark untested paths

- **`--unphased --threads N` panicked** with `RefCell already borrowed`. `with_reader` held
  the thread-local borrow across the caller's closure; `call.rs` is already inside a
  `par_iter` when it takes a reader, and `phase_insertions` starts a nested one, which
  work-stealing can schedule onto the same thread. Load-dependent, invisible on phased input.
  Fixed in `aed861d`. **The unphased path is still unbenchmarked.**
- **Read downsampling was unseeded** (§27.1), confirmed as the cause of 2 of the 3
  differing loci in §21.2, though it fires on only ~0.02% of loci at this depth.

### 29.8 QUICKREF's rationale has inverted

It exists to save time on loci that are homozygous reference anyway. Measured on 50,000 loci
in `fast` mode:

| | CPU seconds | exact | scored |
|---|---|---|---|
| QUICKREF on | 6,628.8 | 85.2% | 46,579 |
| QUICKREF off | **6,319.8** | 85.1% | 46,241 |

**Turning it off is ~5% faster.** It fires on 4.1% of loci — 5.0% of the 82% that are truly
homozygous reference — at 99.0% precision, and on the fast path what it skips costs barely
more than the check itself. It now buys ~0.3 points and 338 rescued loci at a 5% CPU cost,
which is the opposite of its purpose. Whether it belongs to `--mode sensitive` only depends
on the sensitive-mode measurement, still running.

---

## 30. Still outstanding

A deliberate list, so the things not done are as visible as the things done.

### 30.1 The experiment that matters most, and has never been run

**A synthetic expansion ladder** (§27.4): reads carrying 100/500/1000/2000/4000/8000 bp
insertions at one locus, genotyped with defaults. It is the only planned experiment that
addresses *silently undersized* expansions, and no amount of HG002 benchmarking substitutes
for it, because the truth set contains what it contains.

The mechanism it tests: `find_insertions` keeps only `read.is_primary`, and minimap2's
`map_ont` defaults (`zdrop=400`, `max_gap=5000`) can split a sufficiently large insertion
across alignments — the non-primary piece is discarded and the allele comes back **short
rather than absent**. `flanking = 5000` caps it from the other side. Detection measured on
real data (92.2% found, 1.0% `short`) cannot distinguish "the truth set has few enormous
alleles" from "we size them wrongly", and a ladder can.

This is the one item where a missed pathogenic expansion is the failure mode, so it
outranks everything else on this list.

### 30.2 Measured but unexplained, or unmeasured

- **The unphased path has never been benchmarked.** It panicked under `--threads N`
  (§29.7); the crash is fixed, the accuracy is unknown. `find_roots`' dissimilarity
  threshold of 5.0 can discard a minority expansion and split the reference allele in two to
  replace it — a confident `0/0` where an expansion exists, with `EXPANSION_OUTLIER` as the
  only (INFO-only) signal.
- **`remove_outliers` at length-variable loci.** Ruled out for the +1 bin — it is inactive
  on 96.8% of homozygous-reference alleles because their cluster std dev is below 5 — but its
  symmetric ±2 SD trimming on an asymmetric distribution remains a candidate for
  under-calling at long loci, where it *is* active. Never tested.
- **A second non-determinism source.** Of the three loci differing between two runs of
  identical code (§21.2), two were the unseeded downsampler (§27.1). The third had
  `SUP=2,4` against `2,3` — low coverage, so downsampling cannot be the cause, and the
  support count itself changed. Unidentified.
- **The 2m18s arm.** One benchmark arm finished in 2m18s where comparable arms take 13+
  minutes, scoring a full 9,294 loci. Probably page cache on the CRAM after many passes, but
  that is a guess and the number was used.

### 30.3 Owed to others

- **Report the POA endpoint bug upstream.** `Aligner::consensus()` picks its endpoint by
  `max_by_key` over a strictly increasing score, so it is always a sink — wrong for any
  noisy-read use, not just ours. Alongside the existing rust-bio#677. Noted in
  `src/consensus.rs`, not filed.
- **#31 has no comment from us** although §29 largely explains it, and the issue still
  names the similarity threshold as the cause.
- **§12's loose end**: the §4.5 QUICKREF measurement is still recorded only here.

### 30.4 Process debt

- **`misc/` has no tests and CI does not touch it** (§18.5), and it now carries the harness,
  three analysis scripts and the truth-set construction whose bug (§25) silently distorted a
  whole stratum.
- **Everything is measured on one sample**, HG002 at ~30x ONT, one catalog, one chemistry.
  Seeds control locus sampling, not that. The *orderings* in this document should generalise;
  the absolute numbers should not be quoted as universal.
