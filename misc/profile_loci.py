#!/usr/bin/env python3
"""Profile STRdust per-locus CPU cost across replicates and relate it to locus features.

STRdust emits a per-locus ``TIME=<seconds>s`` INFO field in ``--debug`` mode. Since
commit "perf: thread CPU time" that value is *thread CPU time* (CLOCK_THREAD_CPUTIME_ID),
not wall-clock, so it is largely immune to other processes stealing the core. Even so a
single run is noisy, so this script runs STRdust several times and aggregates per locus to
find the ones that are *consistently* slow, then correlates slowness with locus features
parsed from the same VCF (depth, allele length, stdev, cluster count, ...).

Typical use (run on the machine with the data, not the laptop):

    python misc/profile_loci.py \
        --binary ./target/release/STRdust \
        --fasta ref.fa --bam sample.cram --bed loci.bed \
        --replicates 5 --threads 1 \
        --out-prefix profile_run1

Outputs:
    <out-prefix>.per_locus.tsv   one row per locus: median/MAD/CV CPU time + features
    <out-prefix>.replicates.tsv  raw per-locus time for every replicate (long format)
    <out-prefix>.png             (optional, --plot) feature-vs-time scatter panels

Only the standard library is required to *run the replicates and parse*. pandas + scipy are
used for the aggregation/correlation report and are imported lazily, so a parse-only run
(``--from-vcfs``) still works in a bare environment.
"""

import argparse
import subprocess
import sys
import re
import statistics
from pathlib import Path

TIME_RE = re.compile(r"TIME=([0-9.]+)s")


def parse_info(info: str) -> dict:
    """Parse a VCF INFO column into a dict; flags map to True."""
    out = {}
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = True
    return out


def locus_features(line: str) -> dict | None:
    """Extract time + features from one VCF data line, or None if it has no TIME field."""
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 8:
        return None
    chrom, pos, _id, ref, alt, _qual, _filt, info = cols[:8]
    m = TIME_RE.search(info)
    if m is None:
        return None  # not a --debug line / no timing
    info_d = parse_info(info)

    # FORMAT/sample (present for genotyped loci, absent for some missing records)
    fmt = cols[8].split(":") if len(cols) > 8 else []
    sample = cols[9].split(":") if len(cols) > 9 else []
    sample_d = dict(zip(fmt, sample))

    def two_max(field, cast=int):
        """Max of a 2-value 'a,b' FORMAT field, ignoring '.'."""
        vals = [cast(x) for x in sample_d.get(field, "").split(",") if x not in (".", "")]
        return max(vals) if vals else None

    def two_sum(field, cast=int):
        vals = [cast(x) for x in sample_d.get(field, "").split(",") if x not in (".", "")]
        return sum(vals) if vals else None

    alt_alleles = [a for a in alt.split(",") if a != "."]
    outliers = info_d.get("OUTLIERS")

    return {
        "locus": f"{chrom}:{pos}",
        "chrom": chrom,
        "pos": int(pos),
        "time_s": float(m.group(1)),
        "ref_len": len(ref) if ref != "." else None,
        "n_alt": len(alt_alleles),
        "max_alt_len": max((len(a) for a in alt_alleles), default=0),
        # RB = length relative to ref, FRB = full repeat length, SUP = read support
        "max_rb": two_max("RB"),
        "max_frb": two_max("FRB"),
        "total_sup": two_sum("SUP"),
        "max_stdev": max((int(x) for x in info_d.get("STDEV", "").split(",")
                          if x not in (".", "")), default=None),
        "nclusters": int(info_d["NCLUSTERS"]) if "NCLUSTERS" in info_d else None,
        "n_outliers": len(outliers.split(",")) if outliers else 0,
        "quickref": "QUICKREF" in info_d,
        "clusterfailure": "CLUSTERFAILURE" in info_d,
    }


def run_replicate(binary, fasta, bam, bed, region, threads, debug_extra):
    """Run STRdust once with --debug and return its VCF stdout as a list of lines."""
    cmd = [binary, "--debug", "-t", str(threads)]
    if bed:
        cmd += ["-R", bed]
    if region:
        cmd += ["-r", region]
    cmd += list(debug_extra) + [fasta, bam]
    print(f"  $ {' '.join(cmd)}", file=sys.stderr)
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"STRdust failed (exit {proc.returncode}):\n{proc.stderr[-2000:]}")
    return proc.stdout.splitlines()


def collect(lines):
    """Map locus -> features for every timed data line in one VCF."""
    out = {}
    for line in lines:
        if line.startswith("#"):
            continue
        feat = locus_features(line)
        if feat is not None:
            out[feat["locus"]] = feat
    return out


def aggregate(replicates, warmup, out_prefix):
    """Aggregate per-locus times across replicates and write TSV reports."""
    # Drop warmup replicates (cold cache) from the timing aggregation.
    timed = replicates[warmup:]
    if not timed:
        sys.exit("No replicates left after dropping warmup; lower --warmup.")

    # Raw long-format dump: locus, replicate index, time.
    long_path = Path(f"{out_prefix}.replicates.tsv")
    with long_path.open("w") as fh:
        fh.write("locus\treplicate\ttime_s\n")
        for i, rep in enumerate(replicates):
            tag = "warmup" if i < warmup else str(i - warmup)
            for locus, feat in rep.items():
                fh.write(f"{locus}\t{tag}\t{feat['time_s']:.6f}\n")

    # Features come from the last non-warmup replicate (genotypes are deterministic).
    feature_src = timed[-1]
    loci = sorted(feature_src, key=lambda k: (feature_src[k]["chrom"], feature_src[k]["pos"]))

    rows = []
    for locus in loci:
        times = [rep[locus]["time_s"] for rep in timed if locus in rep]
        if not times:
            continue
        med = statistics.median(times)
        mad = statistics.median([abs(t - med) for t in times]) if len(times) > 1 else 0.0
        cv = (statistics.pstdev(times) / med) if len(times) > 1 and med > 0 else 0.0
        f = feature_src[locus]
        rows.append({
            "locus": locus, "chrom": f["chrom"], "pos": f["pos"],
            "n_rep": len(times),
            "median_s": med, "min_s": min(times), "max_s": max(times),
            "mad_s": mad, "cv": cv,
            "ref_len": f["ref_len"], "max_rb": f["max_rb"], "max_frb": f["max_frb"],
            "max_alt_len": f["max_alt_len"], "total_sup": f["total_sup"],
            "max_stdev": f["max_stdev"], "nclusters": f["nclusters"],
            "n_outliers": f["n_outliers"], "n_alt": f["n_alt"],
            "quickref": int(f["quickref"]), "clusterfailure": int(f["clusterfailure"]),
        })

    cols = list(rows[0].keys())
    per_locus_path = Path(f"{out_prefix}.per_locus.tsv")
    with per_locus_path.open("w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r[c] is None else
                               (f"{r[c]:.6f}" if isinstance(r[c], float) else str(r[c]))
                               for c in cols) + "\n")

    print(f"\nWrote {per_locus_path} ({len(rows)} loci) and {long_path}", file=sys.stderr)
    return rows, per_locus_path


def report(rows, top, per_locus_path, plot, out_prefix):
    """Print the consistently-slowest loci and feature correlations; optional scatter plot."""
    consistent = sorted(rows, key=lambda r: r["median_s"], reverse=True)
    print(f"\n=== Top {top} consistently slowest loci (by median thread-CPU time) ===")
    hdr = ("locus", "median_s", "cv", "total_sup", "max_frb", "nclusters", "n_outliers")
    print("\t".join(hdr))
    for r in consistent[:top]:
        print("\t".join(str(r[c] if r[c] is not None else "") for c in hdr))

    # Correlations need scipy/pandas; degrade gracefully if absent.
    try:
        import pandas as pd
        from scipy.stats import spearmanr
    except ImportError:
        print("\n(install pandas + scipy for the correlation report)", file=sys.stderr)
        return

    df = pd.read_csv(per_locus_path, sep="\t")
    numeric = ["ref_len", "max_rb", "max_frb", "max_alt_len", "total_sup",
               "max_stdev", "nclusters", "n_outliers", "n_alt", "quickref"]
    print("\n=== Spearman correlation of features with median CPU time ===")
    print("feature\trho\tp_value\tn")
    corrs = []
    for col in numeric:
        sub = df[[col, "median_s"]].dropna()
        if sub[col].nunique() < 2 or len(sub) < 3:
            continue
        rho, p = spearmanr(sub[col], sub["median_s"])
        corrs.append((col, rho, p, len(sub)))
    for col, rho, p, n in sorted(corrs, key=lambda x: abs(x[1]), reverse=True):
        print(f"{col}\t{rho:+.3f}\t{p:.2e}\t{n}")

    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        feats = [c for c, *_ in sorted(corrs, key=lambda x: abs(x[1]), reverse=True)][:6]
        n = len(feats)
        if n:
            ncol = 3
            nrow = (n + ncol - 1) // ncol
            fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.2 * nrow), squeeze=False)
            for ax, col in zip(axes.flat, feats):
                ax.scatter(df[col], df["median_s"], s=8, alpha=0.4)
                ax.set_xlabel(col)
                ax.set_ylabel("median CPU time (s)")
            for ax in axes.flat[n:]:
                ax.set_visible(False)
            fig.tight_layout()
            png = f"{out_prefix}.png"
            fig.savefig(png, dpi=130)
            print(f"\nWrote {png}", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--binary", default="./target/release/STRdust", help="STRdust binary")
    p.add_argument("--fasta", help="reference fasta")
    p.add_argument("--bam", help="BAM/CRAM file")
    p.add_argument("--bed", help="region bed file (-R)")
    p.add_argument("--region", help="single region chr:start-end (-r)")
    p.add_argument("--replicates", type=int, default=5, help="number of STRdust runs")
    p.add_argument("--warmup", type=int, default=1,
                   help="leading replicates to exclude from timing (cold cache)")
    p.add_argument("--threads", type=int, default=1,
                   help="threads per run; 1 gives the least cache contention")
    p.add_argument("--debug-extra", nargs=argparse.REMAINDER, default=[],
                   help="extra args passed verbatim before FASTA/BAM (e.g. --phasing-strategy dbscan)")
    p.add_argument("--from-vcfs", nargs="+",
                   help="skip running; aggregate these already-produced --debug VCFs instead")
    p.add_argument("--out-prefix", default="strdust_profile", help="output file prefix")
    p.add_argument("--top", type=int, default=25, help="how many slow loci to print")
    p.add_argument("--plot", action="store_true", help="write a feature-vs-time scatter PNG")
    args = p.parse_args()

    if args.from_vcfs:
        replicates = [collect(Path(v).read_text().splitlines()) for v in args.from_vcfs]
        print(f"Parsed {len(replicates)} VCF(s)", file=sys.stderr)
    else:
        if not (args.fasta and args.bam and (args.bed or args.region)):
            p.error("need --fasta, --bam and one of --bed/--region (or use --from-vcfs)")
        if args.threads != 1:
            print(f"warning: --threads {args.threads}: per-thread CPU time stays valid, "
                  "but cache contention between workers adds noise; --threads 1 is cleanest.",
                  file=sys.stderr)
        replicates = []
        for i in range(args.replicates):
            print(f"replicate {i + 1}/{args.replicates}", file=sys.stderr)
            lines = run_replicate(args.binary, args.fasta, args.bam, args.bed,
                                  args.region, args.threads, args.debug_extra)
            replicates.append(collect(lines))

    rows, per_locus_path = aggregate(replicates, args.warmup, args.out_prefix)
    report(rows, args.top, per_locus_path, args.plot, args.out_prefix)


if __name__ == "__main__":
    main()
