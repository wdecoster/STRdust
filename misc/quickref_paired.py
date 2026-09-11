#!/usr/bin/env python3
"""Follow every QUICKREF firing from one run into another, locus by locus.

This is the measurement that separates the two things QUICKREF does: rescuing loci that
full genotyping no-calls, and correcting loci that full genotyping calls wrong. Both are
invisible in the aggregate concordance numbers, which is why ``KNOWN_ISSUES.md`` §4.5 drew
the wrong conclusion from a comparison against STRdust's own genotyper.

It needs two runs of ``benchmark_truth.py`` over the identical loci (same ``--sample-loci``
and ``--seed``): one whose ``per_locus.tsv`` contains ``quickref_reference`` rows, and a
reference run to follow those loci into -- normally the ``--alignment-all`` run, where the
same loci went through full genotyping instead.

    python misc/quickref_paired.py main_50k noqr_50k        # reproduces §16.2
    python misc/quickref_paired.py pr30_50k noqr_50k        # the same for PR #30

Both arguments are ``--out-prefix`` values; ``<prefix>.per_locus.tsv`` is what gets read.
"""

import argparse
import collections
import csv
import sys
from pathlib import Path

QUICKREF_STATUS = "quickref_reference"


def load(prefix, mode):
    """{(chrom, pos): row} for one mode of one run."""
    path = Path(f"{prefix}.per_locus.tsv")
    if not path.exists():
        sys.exit(f"no such file: {path}")
    rows = {}
    with open(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["mode"] != mode:
                continue
            rows[(row["chrom"], row["pos"])] = row
    if not rows:
        sys.exit(f"{path} holds no rows for mode {mode!r}")
    return rows


def is_exact(row):
    """True/False, or None when the run produced no call to compare."""
    if row is None or row["status"] == "no_call":
        return None
    try:
        return abs(float(row["max_abs_diff"])) == 0
    except (ValueError, KeyError):
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("quickref_run", help="out-prefix of the run whose QUICKREF firings to follow")
    parser.add_argument("full_run", help="out-prefix of the run to follow them into (normally --alignment-all)")
    parser.add_argument("--mode", default="sensitive", help="mode column to read in both runs (default: sensitive)")
    args = parser.parse_args()

    left = load(args.quickref_run, args.mode)
    right = load(args.full_run, args.mode)

    shared = set(left) & set(right)
    if len(shared) != len(left) or len(shared) != len(right):
        print(f"warning: {len(left)} vs {len(right)} loci, {len(shared)} shared -- "
              f"the runs did not score the identical loci, so the counts below are partial\n",
              file=sys.stderr)

    fired = [k for k in shared if left[k]["status"] == QUICKREF_STATUS]
    if not fired:
        sys.exit(f"{args.quickref_run} has no {QUICKREF_STATUS} rows -- was it run with --alignment-all?")

    quickref_right = sum(1 for k in fired if is_exact(left[k]))
    full_right = sum(1 for k in fired if is_exact(right[k]))
    full_wrong = sum(1 for k in fired if is_exact(right[k]) is False)
    full_nocall = sum(1 for k in fired if is_exact(right[k]) is None)
    both_wrong = sum(1 for k in fired if not is_exact(left[k]) and is_exact(right[k]) is not True)
    # a rescue is a locus QUICKREF got *right* that the other run failed to call at all;
    # the no-call total also holds loci QUICKREF itself got wrong, which rescue nothing
    rescued = sum(1 for k in fired if is_exact(left[k]) and is_exact(right[k]) is None)
    called = full_right + full_wrong

    n = len(fired)
    print(f"{args.quickref_run} vs {args.full_run}, mode {args.mode}: "
          f"{n} loci answered by QUICKREF\n")
    width = max(len(args.quickref_run), len(args.full_run)) + 26
    for label, count in [
        ("QUICKREF exact", quickref_right),
        (f"exact in {args.full_run}", full_right),
        (f"called in {args.full_run}, wrong", full_wrong),
        (f"no call in {args.full_run}", full_nocall),
        ("QUICKREF wrong as well", both_wrong),
    ]:
        print(f"  {label:<{width}} {count:6d}   {count / n:6.1%}")

    if called:
        print(f"\n  where {args.full_run} produced a call at all: "
              f"{full_right}/{called} = {full_right / called:.1%} exact")
    print(f"  QUICKREF's contribution: {rescued} loci rescued from no-call, "
          f"{full_wrong} corrected")

    truth_classes = collections.Counter(left[k]["truth_gt"] for k in fired)
    print("\n  truth genotype of the loci it fired on:")
    for gt, count in truth_classes.most_common():
        print(f"    {gt:<8} {count:6d}   {count / n:6.1%}")


if __name__ == "__main__":
    main()
