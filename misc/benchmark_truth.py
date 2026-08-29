#!/usr/bin/env python3
"""Benchmark STRdust allele lengths against a truth set, and compare its modes.

Everything STRdust reports about ``--mode fast`` so far is *agreement with*
``--mode sensitive``, which quietly assumes the alignment path is right. This script
measures both modes against an external truth set instead, so the two questions
"do the modes agree" and "which one is closer to the truth" stop being conflated.

Designed to be run where the data already lives (a cluster with HG002), not on a
laptop: it only needs Python's standard library to run STRdust and parse the VCFs.
pandas is imported lazily for the summary tables, and matplotlib only with ``--plot``.

Truth set
---------
Any VCF whose records describe tandem repeat alleles works, as long as REF and ALT
carry the sequences (not symbolic alleles) and the sample column has a GT. The
intended one is the GIAB HG002 tandem repeat benchmark, which ships a VCF plus a BED
of confident regions; pass the BED as ``--confident-bed`` so loci outside it are
excluded rather than silently scored. Grab the current release from the GIAB FTP
site (``ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/``) - the layout there
changes between releases, so the paths are deliberately not hardcoded here.

The truth VCF's own coordinates become the target BED handed to STRdust. That makes
the comparison apples-to-apples: both sides then use the same REF interval, so
STRdust's ``RB`` (consensus length minus REF length, since v0.21.0) is directly
comparable to ``len(ALT) - len(REF)`` from the truth record. Comparing against a
different repeat catalog would reintroduce the boundary mismatch this avoids.

Typical use
-----------
    python misc/benchmark_truth.py \\
        --binary ./target/release/STRdust \\
        --fasta GRCh38.fa --bam HG002.cram \\
        --truth-vcf HG002_TR_benchmark.vcf.gz \\
        --confident-bed HG002_TR_benchmark_regions.bed \\
        --threads 8 --out-prefix hg002_tr

    # quick smoke run on one chromosome before committing to the whole benchmark
    python misc/benchmark_truth.py ... --regions chr1 --max-loci 2000

Outputs
-------
    <out-prefix>.per_locus.tsv   one row per locus per mode: truth vs called lengths
    <out-prefix>.summary.tsv     concordance overall and by stratum
    <out-prefix>.<mode>.vcf      the raw STRdust output, kept for re-analysis
    <out-prefix>.png             (with --plot) called vs truth, per mode
"""

import argparse
import gzip
import statistics
import subprocess
import sys
import time
from pathlib import Path

MODES = ("sensitive", "fast")


def open_maybe_gzip(path):
    """Open a plain or gzipped text file."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_confident_regions(path):
    """Read a BED into {chrom: [(start, end), ...]} with each chromosome's list sorted."""
    if path is None:
        return None
    regions = {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            regions.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for spans in regions.values():
        spans.sort()
    return regions


def inside_confident(regions, chrom, start, end):
    """Whether [start, end) is fully contained in one confident region."""
    if regions is None:
        return True
    spans = regions.get(chrom)
    if not spans:
        return False
    # linear scan is fine: the caller walks the truth VCF in coordinate order and
    # benchmark region lists are small per chromosome
    for region_start, region_end in spans:
        if region_start <= start and end <= region_end:
            return True
        if region_start > end:
            break
    return False


def truth_alleles(ref, alts, genotype):
    """Allele lengths relative to the reference, for the alleles this sample carries.

    Returns None for a genotype that is missing, or that references an allele the
    record does not define. A symbolic ALT (``<...>``) has no sequence to measure,
    so those records are skipped as well.
    """
    if genotype in (".", "./.", ".|."):
        return None
    indices = []
    for part in genotype.replace("|", "/").split("/"):
        if part == ".":
            return None
        indices.append(int(part))
    sequences = [ref] + alts
    lengths = []
    for index in indices:
        if index >= len(sequences):
            return None
        allele = sequences[index]
        if allele.startswith("<") or allele == "*":
            return None
        lengths.append(len(allele) - len(ref))
    return sorted(lengths)


def read_truth(path, confident, regions_filter, max_loci):
    """Parse the truth VCF into {(chrom, pos): {"ref_len", "lengths"}} plus a BED list."""
    loci = {}
    bed = []
    skipped = {"symbolic_or_missing_gt": 0, "outside_confident": 0, "other_chrom": 0}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, _id, ref, alt = fields[0], int(fields[1]), fields[2], fields[3], fields[4]
            if regions_filter and chrom not in regions_filter:
                skipped["other_chrom"] += 1
                continue
            # VCF POS is 1-based; the record covers 0-based [pos - 1, pos - 1 + len(ref))
            start, end = pos - 1, pos - 1 + len(ref)
            if not inside_confident(confident, chrom, start, end):
                skipped["outside_confident"] += 1
                continue
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            lengths = truth_alleles(ref, alt.split(","), sample.get("GT", "."))
            if lengths is None:
                skipped["symbolic_or_missing_gt"] += 1
                continue
            loci[(chrom, pos)] = {"ref_len": len(ref), "lengths": lengths}
            bed.append((chrom, start, end))
            if max_loci and len(bed) >= max_loci:
                break
    return loci, bed, skipped


def write_bed(bed, path):
    with open(path, "w") as handle:
        for chrom, start, end in bed:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def run_strdust(binary, fasta, bam, bed, mode, threads, extra_args, out_vcf):
    """Run one mode and return its CPU seconds (user + system of the child process)."""
    command = [
        str(binary),
        "-R", str(bed),
        "--mode", mode,
        "--threads", str(threads),
        *extra_args,
        str(fasta),
        str(bam),
    ]
    print(f"[benchmark] {' '.join(command)}", file=sys.stderr)
    before = _child_cpu_seconds()
    wall_start = time.monotonic()
    with open(out_vcf, "w") as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.PIPE)
    wall = time.monotonic() - wall_start
    cpu = _child_cpu_seconds() - before
    if result.returncode != 0:
        sys.exit(
            f"STRdust failed in --mode {mode} (exit {result.returncode}):\n"
            + result.stderr.decode(errors="replace")[-4000:]
        )
    return {"cpu_seconds": cpu, "wall_seconds": wall}


def _child_cpu_seconds():
    """CPU seconds consumed by all reaped children so far."""
    import resource

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime


def read_strdust(path):
    """Parse a STRdust VCF into {(chrom, pos): {"rb", "ref_len", "quickref", "support"}}."""
    calls = {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, ref, info = fields[0], int(fields[1]), fields[3], fields[7]
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            rb = sample.get("RB", ".")
            lengths = []
            for value in rb.split(","):
                if value not in (".", ""):
                    lengths.append(int(value))
            calls[(chrom, pos)] = {
                "ref_len": len(ref),
                "lengths": sorted(lengths) if lengths else None,
                "quickref": "QUICKREF" in info,
                "support": sample.get("SUP", "."),
            }
    return calls


def pair_alleles(truth, called):
    """Best pairing of called to truth alleles: both are sorted, so pair positionally.

    A haploid truth genotype against a diploid call (or the reverse) is compared on
    the alleles that exist in both, which is the most that can be said without
    guessing which allele was dropped.
    """
    n = min(len(truth), len(called))
    return list(zip(truth[:n], called[:n]))


def stratum(truth_lengths, ref_len):
    """Bucket a locus by how far its longer allele departs from the reference."""
    longest = max(abs(x) for x in truth_lengths) if truth_lengths else 0
    if longest == 0:
        return "reference"
    for limit, label in ((10, "1-10bp"), (50, "11-50bp"), (200, "51-200bp")):
        if longest <= limit:
            return label
    return ">200bp"


def compare(truth_loci, calls, mode):
    """One row per locus per mode, with the per-allele differences."""
    rows = []
    for key, truth in truth_loci.items():
        call = calls.get(key)
        if call is None:
            rows.append(
                {
                    "chrom": key[0], "pos": key[1], "mode": mode,
                    "ref_len": truth["ref_len"], "stratum": stratum(truth["lengths"], truth["ref_len"]),
                    "truth": ",".join(map(str, truth["lengths"])), "called": "",
                    "max_abs_diff": None, "status": "not_reported", "support": ".",
                }
            )
            continue
        if call["lengths"] is None:
            status = "quickref_reference" if call["quickref"] else "no_call"
            # QUICKREF means STRdust called the locus homozygous reference without
            # aligning, so the called lengths are zeros rather than missing
            called = [0] * len(truth["lengths"]) if call["quickref"] else None
        else:
            status = "called"
            called = call["lengths"]
        if called is None:
            rows.append(
                {
                    "chrom": key[0], "pos": key[1], "mode": mode,
                    "ref_len": truth["ref_len"], "stratum": stratum(truth["lengths"], truth["ref_len"]),
                    "truth": ",".join(map(str, truth["lengths"])), "called": "",
                    "max_abs_diff": None, "status": status, "support": call["support"],
                }
            )
            continue
        pairs = pair_alleles(truth["lengths"], sorted(called))
        diffs = [called_len - truth_len for truth_len, called_len in pairs]
        rows.append(
            {
                "chrom": key[0], "pos": key[1], "mode": mode,
                "ref_len": truth["ref_len"], "stratum": stratum(truth["lengths"], truth["ref_len"]),
                "truth": ",".join(map(str, truth["lengths"])),
                "called": ",".join(map(str, sorted(called))),
                "max_abs_diff": max(abs(d) for d in diffs) if diffs else None,
                "status": status if diffs else "no_shared_allele",
                "support": call["support"],
            }
        )
    return rows


def summarise(rows, timings):
    """Concordance overall and per stratum, as a list of dict rows."""
    import collections

    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[(row["mode"], "all")].append(row)
        grouped[(row["mode"], row["stratum"])].append(row)

    summary = []
    for (mode, stratum_name), group in sorted(grouped.items()):
        scored = [r for r in group if r["max_abs_diff"] is not None]
        n = len(scored)
        if n == 0:
            continue
        diffs = sorted(abs(r["max_abs_diff"]) for r in scored)
        summary.append(
            {
                "mode": mode,
                "stratum": stratum_name,
                "loci_in_truth": len(group),
                "loci_scored": n,
                "not_reported": sum(1 for r in group if r["status"] == "not_reported"),
                "no_call": sum(1 for r in group if r["status"] == "no_call"),
                "exact_pct": round(100 * sum(1 for d in diffs if d == 0) / n, 2),
                "within_1bp_pct": round(100 * sum(1 for d in diffs if d <= 1) / n, 2),
                "within_5bp_pct": round(100 * sum(1 for d in diffs if d <= 5) / n, 2),
                "within_5pct_pct": round(
                    100
                    * sum(
                        1
                        for r in scored
                        if abs(r["max_abs_diff"])
                        <= max(5, 0.05 * max((abs(int(x)) for x in r["truth"].split(",")), default=0))
                    )
                    / n,
                    2,
                ),
                "median_abs_diff": statistics.median(diffs),
                "p90_abs_diff": diffs[int(0.9 * (n - 1))],
                "max_abs_diff": diffs[-1],
                # blank rather than a placeholder number when --from-vcfs skipped the run
                "cpu_seconds": (
                    round(timings[mode]["cpu_seconds"], 1) if mode in timings else None
                ),
            }
        )
    return summary


def write_tsv(rows, path, columns=None):
    if not rows:
        return
    columns = columns or list(rows[0].keys())
    with open(path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join("" if row.get(c) is None else str(row.get(c, "")) for c in columns) + "\n")


def plot(rows, path):
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[benchmark] matplotlib not available, skipping --plot", file=sys.stderr)
        return
    modes = sorted({r["mode"] for r in rows})
    fig, axes = plt.subplots(1, len(modes), figsize=(6 * len(modes), 5.5), squeeze=False)
    for ax, mode in zip(axes[0], modes):
        scored = [r for r in rows if r["mode"] == mode and r["max_abs_diff"] is not None]
        xs, ys = [], []
        for row in scored:
            for truth_len, called_len in zip(row["truth"].split(","), row["called"].split(",")):
                xs.append(int(truth_len))
                ys.append(int(called_len))
        ax.scatter(xs, ys, s=4, alpha=0.25, edgecolors="none")
        if xs:
            lo, hi = min(xs + ys), max(xs + ys)
            ax.plot([lo, hi], [lo, hi], linewidth=1, color="0.4")
        ax.set_title(f"--mode {mode}  (n={len(xs)} alleles)")
        ax.set_xlabel("truth allele length relative to reference (bp)")
        ax.set_ylabel("STRdust RB (bp)")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"[benchmark] wrote {path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--binary", default="./target/release/STRdust", help="STRdust binary")
    parser.add_argument("--fasta", required=True, help="reference genome (indexed)")
    parser.add_argument("--bam", required=True, help="BAM/CRAM for the truth-set sample")
    parser.add_argument("--truth-vcf", required=True, help="truth VCF with sequence-resolved TR alleles")
    parser.add_argument("--confident-bed", help="benchmark regions; loci not fully inside are skipped")
    parser.add_argument("--out-prefix", required=True, help="prefix for the output files")
    parser.add_argument("--modes", default=",".join(MODES), help=f"comma-separated, from {MODES}")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--regions", help="comma-separated chromosomes to restrict to (e.g. chr1)")
    parser.add_argument("--max-loci", type=int, help="stop after this many truth loci (smoke runs)")
    parser.add_argument(
        "--strdust-arg",
        action="append",
        default=[],
        help="extra argument passed through to STRdust, repeatable (e.g. --strdust-arg=--unphased)",
    )
    parser.add_argument("--from-vcfs", action="store_true", help="re-analyse existing output, do not rerun")
    parser.add_argument("--plot", action="store_true", help="write a called-vs-truth scatter per mode")
    args = parser.parse_args()

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    unknown = [m for m in modes if m not in MODES]
    if unknown:
        sys.exit(f"unknown mode(s) {unknown}; expected any of {MODES}")

    regions_filter = set(args.regions.split(",")) if args.regions else None
    confident = load_confident_regions(args.confident_bed)

    print("[benchmark] reading truth set", file=sys.stderr)
    truth_loci, bed, skipped = read_truth(args.truth_vcf, confident, regions_filter, args.max_loci)
    if not truth_loci:
        sys.exit("no usable truth loci: check --regions, --confident-bed and the truth VCF's sample column")
    print(
        f"[benchmark] {len(truth_loci)} truth loci; skipped "
        + ", ".join(f"{v} {k}" for k, v in skipped.items() if v),
        file=sys.stderr,
    )

    bed_path = Path(f"{args.out_prefix}.targets.bed")
    write_bed(bed, bed_path)

    timings = {}
    all_rows = []
    for mode in modes:
        out_vcf = Path(f"{args.out_prefix}.{mode}.vcf")
        if args.from_vcfs:
            if not out_vcf.exists():
                sys.exit(f"--from-vcfs given but {out_vcf} does not exist")
        else:
            timings[mode] = run_strdust(
                args.binary, args.fasta, args.bam, bed_path, mode, args.threads, args.strdust_arg, out_vcf
            )
            print(
                f"[benchmark] --mode {mode}: {timings[mode]['cpu_seconds']:.1f}s CPU, "
                f"{timings[mode]['wall_seconds']:.1f}s wall",
                file=sys.stderr,
            )
        all_rows.extend(compare(truth_loci, read_strdust(out_vcf), mode))

    per_locus = Path(f"{args.out_prefix}.per_locus.tsv")
    write_tsv(all_rows, per_locus)
    summary = summarise(all_rows, timings)
    summary_path = Path(f"{args.out_prefix}.summary.tsv")
    write_tsv(summary, summary_path)

    print(f"\n[benchmark] wrote {per_locus} and {summary_path}\n", file=sys.stderr)
    if summary:
        columns = list(summary[0].keys())
        def cell(value):
            return "" if value is None else str(value)

        widths = [max(len(c), max(len(cell(r[c])) for r in summary)) for c in columns]
        print("  ".join(c.ljust(w) for c, w in zip(columns, widths)))
        for row in summary:
            print("  ".join(cell(row[c]).ljust(w) for c, w in zip(columns, widths)))

    if args.plot:
        plot(all_rows, f"{args.out_prefix}.png")


if __name__ == "__main__":
    main()
