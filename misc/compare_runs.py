#!/usr/bin/env python3
"""Compare benchmark_truth.py runs against each other, on the loci they share.

``benchmark_truth.py`` scores one build. This scores several against each other and prints
the tables that actually drove decisions in ``KNOWN_ISSUES.md`` §13-§14: length and genotype
concordance side by side, the QUICKREF footprint, the consensus length bias at truth
homozygous-reference loci, and a sweep of what the reference/alternate threshold could
achieve.

Runs are named ``label=prefix`` or just ``prefix``; ``prefix`` is what was passed to
``--out-prefix``, so ``<prefix>.per_locus.tsv`` is what gets read. A prefix with several
modes in it (``sensitive`` and ``fast``) is expanded into one run per mode.

The RB correction
-----------------
Before v0.21.0, ``RB`` was one base too high on every allele (#22), so its length columns
are not comparable with a later build's. Pass ``--shift label=1`` to subtract it, or leave
it to ``--auto-shift`` (the default), which reads the run's VCF if it sits next to the TSV
and measures ``RB - (FRB - len(REF))`` directly rather than guessing from a version number.

    python misc/compare_runs.py v0.16=v016_50k main=main_50k
"""

import argparse
import collections
import csv
import statistics
import sys
from pathlib import Path


def measure_shift(prefix, mode):
    """Read RB - (FRB - len(REF)) from the run's VCF; 0 when it cannot be determined."""
    vcf = Path(f"{prefix}.{mode}.vcf")
    if not vcf.exists():
        return 0
    seen = collections.Counter()
    with open(vcf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            rb, frb = sample.get("RB", "."), sample.get("FRB", ".")
            if "." in (rb, frb):
                continue
            for a, b in zip(rb.split(","), frb.split(",")):
                if a == "." or b == ".":
                    continue
                seen[int(a) - (int(b) - len(fields[3]))] += 1
            if sum(seen.values()) > 20000:
                break
    return seen.most_common(1)[0][0] if seen else 0


def load(prefix, shift, auto_shift):
    """{(mode, chrom, pos): row} for every mode present in the per-locus file."""
    path = Path(f"{prefix}.per_locus.tsv")
    if not path.exists():
        sys.exit(f"no such file: {path}")
    shifts = {}
    out = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            mode = row["mode"]
            if mode not in shifts:
                shifts[mode] = shift if shift is not None else (
                    measure_shift(prefix, mode) if auto_shift else 0
                )
            called = None
            if row["called"]:
                # a QUICKREF record reports zeros rather than a measured RB, so it carries
                # no off-by-one to undo
                delta = shifts[mode] if row["status"] == "called" else 0
                called = [int(x) - delta for x in row["called"].split(",")]
            truth = [int(x) for x in row["truth"].split(",")] if row["truth"] else None
            error = None
            if called and truth:
                pairs = min(len(called), len(truth))
                error = max(abs(called[i] - truth[i]) for i in range(pairs))
            out[(mode, row["chrom"], row["pos"])] = {
                "error": error,
                "called": called,
                "truth": truth,
                "truth_gt": row["truth_gt"],
                "called_gt": row["called_gt"],
                "status": row["status"],
                "ref_len": int(row["ref_len"]),
                "stratum": row["stratum"],
            }
    return out, shifts


def concordance(runs, shared):
    print("\nlength and genotype concordance, on loci common to every run")
    header = f"{'run':>22} {'scored':>7} {'exact':>7} {'<=1bp':>7} {'<=5bp':>7} {'GT':>7} {'recall':>7} {'prec':>7}"
    print(header)
    for name, rows in runs.items():
        sub = [rows[key] for key in shared]
        scored = [r for r in sub if r["error"] is not None]
        gt = [r for r in sub if r["truth_gt"] and r["called_gt"]]
        truth_var = [r for r in gt if r["truth_gt"] != "0/0"]
        called_var = [r for r in gt if r["called_gt"] != "0/0"]
        hits = sum(1 for r in truth_var if r["called_gt"] != "0/0")
        pct = lambda n, d: f"{100 * n / d:6.1f}%" if d else "      -"
        print(
            f"{name:>22} {len(scored):>7} "
            f"{pct(sum(1 for r in scored if r['error'] == 0), len(scored))} "
            f"{pct(sum(1 for r in scored if r['error'] <= 1), len(scored))} "
            f"{pct(sum(1 for r in scored if r['error'] <= 5), len(scored))} "
            f"{pct(sum(1 for r in gt if r['called_gt'] == r['truth_gt']), len(gt))} "
            f"{pct(hits, len(truth_var))} {pct(hits, len(called_var))}"
        )


def quickref(runs, shared):
    print("\nQUICKREF footprint, and over-calling among the loci it did not catch")
    print(f"{'run':>22} {'fires':>8} {'of all':>7} {'correct':>8} {'over-call if genotyped':>23}")
    for name, rows in runs.items():
        sub = [rows[key] for key in shared]
        fired = [r for r in sub if r["status"] == "quickref_reference"]
        right = sum(1 for r in fired if r["truth_gt"] == "0/0")
        genotyped = [r for r in sub if r["status"] == "called" and r["truth_gt"] and r["called_gt"]]
        hom = [r for r in genotyped if r["truth_gt"] == "0/0"]
        over = sum(1 for r in hom if r["called_gt"] != "0/0")
        print(
            f"{name:>22} {len(fired):>8} {100 * len(fired) / len(sub):>6.1f}% "
            f"{(100 * right / len(fired) if fired else 0):>7.1f}% "
            f"{(100 * over / len(hom) if hom else 0):>22.1f}%"
        )


def consensus_bias(runs, shared):
    print("\nconsensus length error at truth homozygous-reference loci (correct answer: 0)")
    print(f"{'run':>22} {'alleles':>8} {'exact':>7} {'too long':>9} {'too short':>10} {'median |err|':>13}")
    for name, rows in runs.items():
        deltas = []
        for key in shared:
            row = rows[key]
            if row["status"] != "called" or row["truth_gt"] != "0/0" or not row["called"]:
                continue
            deltas.extend(row["called"])
        if not deltas:
            continue
        exact = sum(1 for d in deltas if d == 0)
        long_ = sum(1 for d in deltas if d > 0)
        short = sum(1 for d in deltas if d < 0)
        n = len(deltas)
        print(
            f"{name:>22} {n:>8} {100 * exact / n:>6.1f}% {100 * long_ / n:>8.1f}% "
            f"{100 * short / n:>9.1f}% {statistics.median(abs(d) for d in deltas):>13}"
        )


def expansion_detection(runs, shared, thresholds=(50, 100, 200, 500)):
    """Did the caller *find* the expansion, regardless of whether it sized it correctly.

    For a pathogenic repeat the clinically relevant failure is missing an expansion, not
    reporting it a few bases short: a no-call and a 3 bp error are not the same kind of
    mistake, and exact-length concordance scores them as if the second were worse. This
    table asks the detection question instead, over loci whose truth carries a *positive*
    (expanded) allele of at least T:

      found      longest called allele >= T        - sized right and detected
      partial    longest called allele >= T/2      - undersized but unmistakably expanded
      short      called, but reported under T/2    - a false negative that looks like a call
      no-call    nothing reported                  - a false negative that is at least visible

    `found` and `partial` together are what matters for screening; `short` is the dangerous
    column, because a confidently wrong small number does not prompt a second look.
    """
    print("\nexpansion detection: loci whose truth has an expanded allele >= T")
    print(f"{'run':>22}{'T':>7}{'loci':>7}{'found':>9}{'+partial':>10}{'short':>8}{'no-call':>9}")
    for name, rows in runs.items():
        for threshold in thresholds:
            here = [rows[k] for k in shared
                    if rows[k]["truth"] and max(rows[k]["truth"]) >= threshold]
            if not here:
                continue
            found = partial = short = missing = 0
            for row in here:
                if not row["called"]:
                    missing += 1
                    continue
                longest = max(row["called"])
                if longest >= threshold:
                    found += 1
                elif longest >= threshold / 2:
                    partial += 1
                else:
                    short += 1
            n = len(here)
            print(f"{name:>22}{threshold:>7}{n:>7}{100 * found / n:>8.1f}%"
                  f"{100 * (found + partial) / n:>9.1f}%{100 * short / n:>7.1f}%"
                  f"{100 * missing / n:>8.1f}%")
        print()


STRATUM_ORDER = ["reference", "1-10bp", "11-50bp", "51-200bp", ">200bp"]


def by_truth_size(runs, shared):
    """Concordance split by how far the truth allele departs from the reference.

    The stratum comes from the *truth*, so it is the same set of loci for every run and the
    columns are directly comparable. This is the table that says whether a change bought
    accuracy at small allele differences by giving up long expansions -- which is what
    STRdust is for, and the thing the published benchmark found it best at. A knob that
    improves the aggregate while the >200bp column falls is not an improvement.
    """
    print("\nconcordance by truth allele size (exact / <=5bp / no-call, per stratum)")
    strata = [s for s in STRATUM_ORDER
              if any(rows[k]["stratum"] == s for rows in runs.values() for k in shared)]
    if not strata:
        return
    header = "".join(f"{s:>22}" for s in strata)
    print(f"{'run':>22}{header}")
    for name, rows in runs.items():
        cells = ""
        for stratum in strata:
            here = [rows[k] for k in shared if rows[k]["stratum"] == stratum]
            scored = [r for r in here if r["error"] is not None]
            if not here:
                cells += f"{'-':>22}"
                continue
            if not scored:
                cells += f"{'0 scored':>22}"
                continue
            exact = 100 * sum(1 for r in scored if r["error"] == 0) / len(scored)
            near = 100 * sum(1 for r in scored if r["error"] <= 5) / len(scored)
            cells += f"{exact:>7.1f}%{near:>7.1f}%{len(here) - len(scored):>7}"
        print(f"{name:>22}{cells}")
    print(f"{'':>22}" + "".join(f"{'n=' + str(sum(1 for k in shared if next(iter(runs.values()))[k]['stratum'] == s)):>22}"
                                for s in strata))


def error_spectrum(runs, shared, span=5):
    """How the consensus length error is distributed, base by base.

    The aggregate "too long" percentage hides where the mass sits. It sits at +1 and +2:
    a single base is the most common wrong answer, which is what makes the junction knobs
    (`--minlen`, the +/-30 fold window) worth measuring one base at a time -- see §19/§20.
    """
    print(f"\nconsensus error spectrum at truth hom-ref loci, -{span}..+{span} bp "
          f"(share of alleles; '>' is everything beyond)")
    header = "".join(f"{d:>+7d}" for d in range(-span, span + 1))
    print(f"{'run':>22} {'alleles':>8}   {'<':>6}{header}{'>':>7}    {'+1 as % of over':>16}")
    for name, rows in runs.items():
        deltas = []
        for key in shared:
            row = rows[key]
            if row["status"] != "called" or row["truth_gt"] != "0/0" or not row["called"]:
                continue
            deltas.extend(row["called"])
        if not deltas:
            continue
        n = len(deltas)
        counts = collections.Counter(deltas)
        below = sum(c for d, c in counts.items() if d < -span)
        above = sum(c for d, c in counts.items() if d > span)
        over = sum(c for d, c in counts.items() if d > 0)
        cells = "".join(f"{100 * counts.get(d, 0) / n:>6.1f}%" for d in range(-span, span + 1))
        share = f"{100 * counts.get(1, 0) / over:>15.1f}%" if over else f"{'-':>16}"
        print(f"{name:>22} {n:>8}   {100 * below / n:>5.1f}%{cells}{100 * above / n:>6.1f}%  {share}")


def by_length(runs, shared):
    print("\nover-call / under-call by locus length (the len(REF)/20 threshold in action)")
    bins = ((0, 20), (20, 50), (50, 100), (100, 200), (200, 500), (500, 10**9))
    for name, rows in runs.items():
        print(f"  {name}")
        print(f"    {'ref_len':>12} {'loci':>7} {'threshold':>10} {'over-call':>10} {'under-call':>11}")
        for low, high in bins:
            sub = [
                rows[key]
                for key in shared
                if low <= rows[key]["ref_len"] < high
                and rows[key]["status"] == "called"
                and rows[key]["truth_gt"]
                and rows[key]["called_gt"]
            ]
            if len(sub) < 50:
                continue
            hom = [r for r in sub if r["truth_gt"] == "0/0"]
            var = [r for r in sub if r["truth_gt"] != "0/0"]
            over = 100 * sum(1 for r in hom if r["called_gt"] != "0/0") / len(hom) if hom else 0
            under = 100 * sum(1 for r in var if r["called_gt"] == "0/0") / len(var) if var else 0
            median_len = statistics.median(r["ref_len"] for r in sub)
            label = f"{low}-{high if high < 10**9 else ''}"
            print(f"    {label:>12} {len(sub):>7} {int(median_len) // 20:>10} {over:>9.1f}% {under:>10.1f}%")


def threshold_sweep(runs, shared):
    print("\nwhat the reference/alternate decision could reach as a pure length rule")
    print("(is_similar_to_ref uses edit distance; |length difference| is its dominant term)")
    rules = [("current: len/20", lambda rl: rl // 20)]
    rules += [(f"constant {t} bp", lambda rl, t=t: t) for t in (1, 2, 3, 5)]
    rules += [("min(4, len/20)", lambda rl: min(4, rl // 20))]
    for name, rows in runs.items():
        data = [
            (rows[key]["ref_len"], rows[key]["called"], rows[key]["truth_gt"] != "0/0")
            for key in shared
            if rows[key]["status"] == "called" and rows[key]["called"] and rows[key]["truth_gt"]
        ]
        if not data:
            continue
        print(f"  {name}  ({len(data)} genotyped loci)")
        print(f"    {'rule':>18} {'over-call':>10} {'under-call':>11} {'balanced acc':>13}")
        for label, rule in rules:
            tp = fp = tn = fn = 0
            for ref_len, called, is_var in data:
                limit = rule(ref_len)
                called_var = any(abs(x) > limit for x in called)
                if is_var and called_var:
                    tp += 1
                elif is_var:
                    fn += 1
                elif called_var:
                    fp += 1
                else:
                    tn += 1
            if not (tp + fn) or not (fp + tn):
                continue
            balanced = 50 * (tp / (tp + fn) + tn / (tn + fp))
            print(
                f"    {label:>18} {100 * fp / (fp + tn):>9.1f}% {100 * fn / (tp + fn):>10.1f}% "
                f"{balanced:>12.1f}%"
            )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("runs", nargs="+", help="out-prefix of a run, optionally label=prefix")
    parser.add_argument(
        "--shift",
        action="append",
        default=[],
        help="label=N: subtract N from every RB in that run (#22 was +1 before v0.21.0). "
        "Overrides the measured value",
    )
    parser.add_argument(
        "--no-auto-shift",
        action="store_true",
        help="do not measure the RB offset from each run's VCF",
    )
    args = parser.parse_args()

    overrides = {}
    for item in args.shift:
        label, _, value = item.partition("=")
        overrides[label] = int(value)

    runs = {}
    for item in args.runs:
        label, _, prefix = item.partition("=")
        if not prefix:
            label, prefix = item, item
        rows, shifts = load(prefix, overrides.get(label), not args.no_auto_shift)
        modes = sorted({mode for mode, _chrom, _pos in rows})
        for mode in modes:
            name = f"{label} {mode}" if len(modes) > 1 else label
            runs[name] = {
                (chrom, pos): value
                for (row_mode, chrom, pos), value in rows.items()
                if row_mode == mode
            }
            note = f"  RB shift {shifts[mode]:+d}" if shifts[mode] else ""
            print(f"loaded {name}: {len(runs[name])} loci{note}", file=sys.stderr)

    shared = set.intersection(*(set(rows) for rows in runs.values()))
    print(f"\n{len(shared)} loci common to all {len(runs)} runs")

    concordance(runs, shared)
    quickref(runs, shared)
    by_truth_size(runs, shared)
    expansion_detection(runs, shared)
    consensus_bias(runs, shared)
    error_spectrum(runs, shared)
    by_length(runs, shared)
    threshold_sweep(runs, shared)


if __name__ == "__main__":
    main()
