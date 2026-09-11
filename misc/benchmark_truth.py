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

If the truth VCF is a general call set that also contains SNVs and other non-repeat
variants, pass ``--repeat-bed`` with a BED of repeat loci. Each BED locus then becomes one
STRdust target, and its truth alleles are the summed per-haplotype length change of the
variants inside it - so an SNV in the repeat contributes 0, and a locus with no variant at
all (inside the confident regions) is truth-homozygous-reference, which is how the large
invariant class enters the benchmark. Loci whose phase cannot be resolved - two or more
unphased heterozygous length-changing records - are skipped rather than guessed.

Without ``--repeat-bed``, the truth VCF's own coordinates become the target BED handed to
STRdust. That makes
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

    # a release older than v0.21.0, which has no --mode flag (single alignment path)
    python misc/benchmark_truth.py ... --binary ./STRdust-v0.16.0b \
        --no-mode-flag --modes sensitive

Note that this is *not* a reimplementation of any published benchmark: it pairs alleles
positionally, scores exact/1bp/5bp length tolerances rather than motif-aware deviation
classes, and knows nothing about motifs or sequence-level (Levenshtein) similarity.

Outputs
-------
    <out-prefix>.per_locus.tsv   one row per locus per mode: truth vs called lengths
    <out-prefix>.summary.tsv     concordance overall and by stratum
    <out-prefix>.<mode>.vcf      the raw STRdust output, kept for re-analysis
    <out-prefix>.png             (with --plot) called vs truth, per mode
"""

import argparse
import bisect
import gzip
import hashlib
import json
import math
import random
import statistics
import subprocess
import sys
import time
from pathlib import Path

MODES = ("sensitive", "fast")
# how far back to scan for a long locus whose start precedes the variant
MAX_LOCUS_SPAN = 100_000


def open_maybe_gzip(path):
    """Open a plain or gzipped text file."""
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def load_confident_regions(path):
    """Read a BED into {chrom: (starts, max_end_so_far)}, both sorted by start.

    ``max_end_so_far[i]`` is the largest end among regions 0..i. That is what makes
    ``inside_confident`` a bisect instead of a scan: a locus is contained in some region iff
    the largest end seen up to the last region starting at or before it reaches past the
    locus. Storing it costs one extra list per chromosome and removes a quadratic factor
    over the repeat catalog -- see ``inside_confident``.
    """
    if path is None:
        return None
    raw = {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            raw.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    regions = {}
    for chrom, spans in raw.items():
        spans.sort()
        starts = [start for start, _end in spans]
        max_ends = []
        running = 0
        for _start, end in spans:
            running = end if end > running else running
            max_ends.append(running)
        regions[chrom] = (starts, max_ends)
    return regions


def inside_confident(regions, chrom, start, end):
    """Whether [start, end) is fully contained in one confident region.

    Only regions starting at or before ``start`` can contain the locus, and among those it
    is contained iff the largest end reaches ``end`` -- that region then satisfies both
    bounds. ``max_end_so_far`` makes this two bisected lookups rather than a walk from the
    start of the chromosome, which matters because ``read_truth_over_repeats`` asks this
    once per catalog locus: ~1.8M times against a comparable number of regions, where a
    scan whose length grows with position is quadratic and costs hours.
    """
    if regions is None:
        return True
    entry = regions.get(chrom)
    if not entry:
        return False
    starts, max_ends = entry
    index = bisect.bisect_right(starts, start) - 1
    if index < 0:
        return False
    return end <= max_ends[index]


def genotype_class(genotype):
    """Normalise a GT into a class: "0/0", "0/1", "1/1", "1/2", haploid "0"/"1", or None.

    Classes follow the usual convention: distinct non-reference allele indices give 1/2,
    the same index twice gives 1/1. Phasing and allele order are ignored.
    """
    if not genotype or genotype in (".", "./.", ".|."):
        return None
    parts = genotype.replace("|", "/").split("/")
    if any(part == "." for part in parts):
        return None
    try:
        indices = [int(part) for part in parts]
    except ValueError:
        return None
    if not indices:
        return None
    if len(indices) == 1:
        return "0" if indices[0] == 0 else "1"
    first, second = sorted(indices)[:2]
    if first == 0 and second == 0:
        return "0/0"
    if first == 0:
        return "0/1"
    return "1/1" if first == second else "1/2"


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


def read_truth(path, confident, regions_filter, max_loci, sample_loci=None, seed=0):
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
            loci[(chrom, pos)] = {
                "ref_len": len(ref),
                "lengths": lengths,
                "gt_class": genotype_class(sample.get("GT", ".")),
            }
            bed.append((chrom, start, end))
            if max_loci and not sample_loci and len(bed) >= max_loci:
                break
    loci, bed = subset_loci(loci, bed, max_loci, sample_loci, seed)
    return loci, bed, skipped


def load_repeat_loci(path):
    """Repeat catalog BED -> {chrom: [(start, end), ...]} sorted, 0-based half-open."""
    loci = {}
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")):
                continue
            fields = line.split()
            if len(fields) < 3:
                continue
            loci.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
    for spans in loci.values():
        spans.sort()
    return loci


def overlapping_loci(spans, starts, var_start, var_end):
    """Indices of loci in `spans` overlapping the 0-based half-open [var_start, var_end)."""
    hits = []
    index = bisect.bisect_right(starts, var_end) - 1
    while index >= 0:
        locus_start, locus_end = spans[index]
        if locus_end <= var_start:
            # loci are sorted by start; an earlier one may still be long enough to reach,
            # so only stop once we are clearly past any plausible overlap
            if var_start - locus_start > MAX_LOCUS_SPAN:
                break
        elif locus_start < var_end:
            hits.append(index)
        index -= 1
    return hits


def haplotype_deltas(ref, alts, genotype):
    """Per-haplotype (length delta, is_reference) for one record, or None if unusable."""
    if genotype in (".", "./.", ".|."):
        return None
    phased = "|" in genotype
    parts = genotype.replace("|", "/").split("/")
    if any(part == "." for part in parts):
        return None
    sequences = [ref] + alts
    out = []
    for part in parts:
        index = int(part)
        if index >= len(sequences):
            return None
        allele = sequences[index]
        if allele.startswith("<") or allele == "*":
            return None
        out.append((len(allele) - len(ref), index == 0))
    return out, phased


def clamp_to_locus(parsed, var_start, locus_start, locus_end):
    """Restrict a record's deleted bases to the ones actually inside the locus.

    A variant is attributed to a repeat when it *overlaps* it, but a deletion reaching well
    past the repeat only removes the repeat bases it covers. Counting its whole length
    produced truth alleles like a 61 bp locus contracting by 145,777 bp, which no genotyper
    can reproduce and which STRdust correctly no-called -- scored as failure. The artefacts
    were rare overall (0.15% of loci) but concentrated entirely in the long-expansion
    stratum, where they were ~40% of the loci and nearly all of its no-calls.

    Insertions are not clamped: the inserted bases enter at a position inside the locus, so
    their whole length belongs to it. Deleted bases start one base after the VCF position,
    the usual anchor-base convention.
    """
    deltas, phased = parsed
    clamped = []
    for delta, is_ref in deltas:
        if delta < 0:
            deleted_start = var_start + 1
            deleted_end = deleted_start - delta
            overlap = min(locus_end, deleted_end) - max(locus_start, deleted_start)
            delta = -max(0, overlap)
        clamped.append((delta, is_ref))
    return clamped, phased


def read_truth_over_repeats(
    vcf_path, repeat_bed, confident, regions_filter, max_loci, sample_loci=None, seed=0
):
    """Truth allele lengths per *repeat locus*, summing the variants inside each one.

    Use this when the truth VCF is a general variant call set (SNVs, indels, ...) rather
    than one record per tandem repeat: the repeat BED says which intervals are repeats,
    each locus becomes one STRdust target, and its truth alleles are the total length
    change of the variants on each haplotype. A locus inside the confident regions with no
    variant record is truth-homozygous-reference, which is how the ~84% invariant class
    enters the benchmark at all.

    Loci are skipped when phase cannot be resolved: two or more heterozygous
    length-changing records with unphased genotypes have no unique assignment to
    haplotypes, and guessing would invent truth.
    """
    catalog = load_repeat_loci(repeat_bed)
    per_locus = {}
    skipped = {
        "symbolic_or_missing_gt": 0,
        "outside_confident": 0,
        "other_chrom": 0,
        "unphased_multivariant": 0,
        "no_overlapping_locus": 0,
    }

    starts = {chrom: [span[0] for span in spans] for chrom, spans in catalog.items()}
    with open_maybe_gzip(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, ref, alt = fields[0], int(fields[1]), fields[3], fields[4]
            if regions_filter and chrom not in regions_filter:
                continue
            spans = catalog.get(chrom)
            if not spans:
                continue
            var_start, var_end = pos - 1, pos - 1 + len(ref)
            hits = overlapping_loci(spans, starts[chrom], var_start, var_end)
            if not hits:
                skipped["no_overlapping_locus"] += 1
                continue
            sample = dict(zip(fields[8].split(":"), fields[9].split(":")))
            parsed = haplotype_deltas(ref, alt.split(","), sample.get("GT", "."))
            if parsed is None:
                skipped["symbolic_or_missing_gt"] += 1
                for index in hits:
                    per_locus.setdefault((chrom, index), []).append(None)
                continue
            for index in hits:
                locus_start, locus_end = spans[index]
                per_locus.setdefault((chrom, index), []).append(
                    clamp_to_locus(parsed, var_start, locus_start, locus_end)
                )

    loci = {}
    bed = []
    for chrom, spans in catalog.items():
        if regions_filter and chrom not in regions_filter:
            skipped["other_chrom"] += len(spans)
            continue
        for index, (start, end) in enumerate(spans):
            if not inside_confident(confident, chrom, start, end):
                skipped["outside_confident"] += 1
                continue
            records = per_locus.get((chrom, index), [])
            if any(record is None for record in records):
                continue  # already counted as symbolic_or_missing_gt
            resolved = resolve_locus(records)
            if resolved is None:
                skipped["unphased_multivariant"] += 1
                continue
            lengths, gt_class = resolved
            # STRdust turns a BED start into VCF POS by adding 1 (repeats.rs::from_bed)
            loci[(chrom, start + 1)] = {
                "ref_len": end - start,
                "lengths": lengths,
                "gt_class": gt_class,
            }
            bed.append((chrom, start, end))
    bed.sort()
    loci, bed = subset_loci(loci, bed, max_loci, sample_loci, seed)
    return loci, bed, skipped


def resolve_locus(records):
    """Total per-haplotype length change for one locus, or None if phase is ambiguous."""
    if not records:
        return [0, 0], "0/0"
    ploidy = max(len(deltas) for deltas, _phased in records)
    # only *length-changing* heterozygous records need phase: a het SNV inside the repeat
    # contributes 0 to both haplotypes whichever way it is assigned
    ambiguous = [
        deltas
        for deltas, phased in records
        if not phased and len({delta for delta, _is_ref in deltas}) > 1
    ]
    if len(ambiguous) > 1:
        return None
    totals = [0] * ploidy
    for deltas, _phased in records:
        for hap in range(ploidy):
            delta, _allele_is_ref = deltas[hap % len(deltas)]
            totals[hap] += delta
    # the genotype class is derived from length, to match the axis this script scores:
    # a haplotype whose total length change is 0 counts as reference, so a tool is not
    # penalised for an SNV inside the repeat that it cannot express as a length
    if ploidy == 1:
        return totals, "0" if totals[0] == 0 else "1"
    if totals[0] == 0 and totals[1] == 0:
        gt_class = "0/0"
    elif totals[0] == 0 or totals[1] == 0:
        gt_class = "0/1"
    else:
        gt_class = "1/1" if totals[0] == totals[1] else "1/2"
    return sorted(totals), gt_class


def subset_loci(loci, bed, max_loci, sample_loci, seed):
    """Trim to `max_loci` in coordinate order, or to a random sample of `sample_loci`.

    Random sampling matters for a partial benchmark: taking the first N loci of a catalog
    samples one end of one chromosome, where repeat composition is not representative.
    """
    if sample_loci and len(bed) > sample_loci:
        rng = random.Random(seed)
        bed = sorted(rng.sample(bed, sample_loci))
    elif max_loci and len(bed) > max_loci:
        bed = bed[:max_loci]
    else:
        return loci, bed
    keep = {(chrom, start + 1) for chrom, start, _end in bed}
    return {key: value for key, value in loci.items() if key in keep}, bed


def write_bed(bed, path):
    with open(path, "w") as handle:
        for chrom, start, end in bed:
            handle.write(f"{chrom}\t{start}\t{end}\n")


def run_strdust(binary, fasta, bam, bed, mode, threads, extra_args, out_vcf, mode_flag=True):
    """Run one mode and return its CPU seconds (user + system of the child process).

    ``mode_flag=False`` omits ``--mode``, which only exists from v0.21.0 onwards. Older
    releases have a single (alignment) path, equivalent to ``sensitive``.
    """
    command = [
        str(binary),
        "-R", str(bed),
        *(("--mode", mode) if mode_flag else ()),
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


def binary_fingerprint(path):
    """md5 and size of the binary actually executed, or None when it cannot be read.

    Recorded with every run because a benchmark series compares *builds*, and a build is
    only identifiable by its content: a staged binary can be replaced, rebuilt from the
    wrong branch, or renamed, and the VCF header records `argv[0]`, which proves the
    filename and nothing else. Two runs of this harness that disagree are only interpretable
    if each says which bytes it ran.
    """
    try:
        resolved = Path(path).resolve()
        digest = hashlib.md5()
        with open(resolved, "rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return {"md5": digest.hexdigest(), "size": resolved.stat().st_size, "path": str(resolved)}
    except OSError as problem:
        print(f"[benchmark] could not fingerprint {path}: {problem}", file=sys.stderr)
        return None


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
                "gt_class": genotype_class(sample.get("GT", ".")),
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
                    "truth_gt": truth["gt_class"], "called_gt": None,
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
                    "truth_gt": truth["gt_class"], "called_gt": call["gt_class"],
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
                "truth_gt": truth["gt_class"],
                "called_gt": call["gt_class"],
                "max_abs_diff": max(abs(d) for d in diffs) if diffs else None,
                "status": status if diffs else "no_shared_allele",
                "support": call["support"],
            }
        )
    return rows


def summarise(rows, timings, fingerprint=None):
    """Concordance overall, per stratum and per truth genotype class, as a list of rows.

    Two axes are reported side by side, because a change can move one without the other:

    - **length** (``exact_pct`` and the tolerance bands) scores ``RB`` against the truth
      allele lengths. This is what a consensus-length change moves.
    - **genotype class** (``gt_match_pct``) scores 0/0 vs 0/1 vs 1/1 vs 1/2. This is what a
      change to ``is_similar_to_ref`` moves, and it is invisible on the length axis, since
      that function does not touch the consensus.

    Rows prefixed ``gt=`` stratify by the *truth* genotype class, which is where the
    published benchmark found STRdust weakest (0/0 loci, ~84% of any catalog).
    """
    import collections

    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[(row["mode"], "all")].append(row)
        grouped[(row["mode"], row["stratum"])].append(row)
        if row.get("truth_gt"):
            grouped[(row["mode"], f"gt={row['truth_gt']}")].append(row)

    summary = []
    for (mode, stratum_name), group in sorted(grouped.items()):
        scored = [r for r in group if r["max_abs_diff"] is not None]
        n = len(scored)
        gt_scored = [r for r in group if r.get("truth_gt") and r.get("called_gt")]
        gt_n = len(gt_scored)
        if n == 0 and gt_n == 0:
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
                "exact_pct": round(100 * sum(1 for d in diffs if d == 0) / n, 2) if n else None,
                "within_1bp_pct": round(100 * sum(1 for d in diffs if d <= 1) / n, 2) if n else None,
                "within_5bp_pct": round(100 * sum(1 for d in diffs if d <= 5) / n, 2) if n else None,
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
                )
                if n
                else None,
                # genotype-class axis: how often the called class equals the truth class,
                # over loci where both sides have a usable GT
                "gt_scored": gt_n,
                "gt_match_pct": round(
                    100 * sum(1 for r in gt_scored if r["called_gt"] == r["truth_gt"]) / gt_n, 2
                )
                if gt_n
                else None,
                "median_abs_diff": statistics.median(diffs) if n else None,
                # nearest-rank p90: int(0.9 * (n - 1)) could fall below the median at small n
                "p90_abs_diff": diffs[math.ceil(0.9 * n) - 1] if n else None,
                "max_abs_diff": diffs[-1] if n else None,
                # blank rather than a placeholder number when --from-vcfs skipped the run
                "cpu_seconds": (
                    round(timings[mode]["cpu_seconds"], 1) if mode in timings else None
                ),
                "binary_md5": fingerprint["md5"] if fingerprint else None,
            }
        )
    return summary


def confusion(rows):
    """Truth genotype class against called class, one row per (mode, truth, called)."""
    import collections

    counts = collections.Counter()
    for row in rows:
        if not row.get("truth_gt"):
            continue
        counts[(row["mode"], row["truth_gt"], row.get("called_gt") or "missing")] += 1
    totals = collections.Counter()
    for (mode, truth_gt, _called), count in counts.items():
        totals[(mode, truth_gt)] += count
    return [
        {
            "mode": mode,
            "truth_gt": truth_gt,
            "called_gt": called_gt,
            "loci": count,
            "pct_of_truth_class": round(100 * count / totals[(mode, truth_gt)], 2),
        }
        for (mode, truth_gt, called_gt), count in sorted(counts.items())
    ]


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


# bump when anything that changes the *content* of the truth set changes: the parsing, the
# confident-region semantics, the sampling. A pure speed-up does not need a bump.
TRUTH_CACHE_VERSION = 2


def truth_cache_key(args, regions_filter):
    """Everything that determines the truth set, so a stale cache is never reused.

    Input files are identified by path, size and mtime rather than by content: hashing a
    multi-gigabyte VCF on every run would cost more than the work being cached.
    """

    def identify(path):
        if path is None:
            return None
        resolved = Path(path).resolve()
        stat = resolved.stat()
        return {"path": str(resolved), "size": stat.st_size, "mtime_ns": stat.st_mtime_ns}

    return {
        "version": TRUTH_CACHE_VERSION,
        "truth_vcf": identify(args.truth_vcf),
        "repeat_bed": identify(args.repeat_bed),
        "confident_bed": identify(args.confident_bed),
        "regions": sorted(regions_filter) if regions_filter else None,
        "max_loci": args.max_loci,
        "sample_loci": args.sample_loci,
        "seed": args.seed,
    }


def load_truth_cache(path, key):
    """(loci, bed, skipped) from the cache, or None when it is absent, stale or unreadable."""
    if path is None or not Path(path).exists():
        return None
    try:
        with open(path) as handle:
            blob = json.load(handle)
    except (OSError, ValueError) as problem:
        print(f"[benchmark] ignoring unreadable truth cache {path}: {problem}", file=sys.stderr)
        return None
    if blob.get("key") != key:
        print(f"[benchmark] truth cache {path} is stale, rebuilding", file=sys.stderr)
        return None
    loci = {
        (chrom, pos): {"ref_len": ref_len, "lengths": lengths, "gt_class": gt_class}
        for chrom, pos, ref_len, lengths, gt_class in blob["loci"]
    }
    bed = [tuple(entry) for entry in blob["bed"]]
    return loci, bed, blob["skipped"]


def save_truth_cache(path, key, loci, bed, skipped):
    """Write the cache atomically, so an interrupted run cannot leave a truncated one."""
    payload = {
        "key": key,
        "loci": [
            [chrom, pos, value["ref_len"], value["lengths"], value["gt_class"]]
            for (chrom, pos), value in loci.items()
        ],
        "bed": [list(entry) for entry in bed],
        "skipped": skipped,
    }
    temporary = Path(f"{path}.tmp")
    try:
        with open(temporary, "w") as handle:
            json.dump(payload, handle)
        temporary.replace(path)
    except OSError as problem:
        print(f"[benchmark] could not write truth cache {path}: {problem}", file=sys.stderr)
        return
    print(f"[benchmark] wrote truth cache {path}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--binary", default="./target/release/STRdust", help="STRdust binary")
    parser.add_argument("--fasta", required=True, help="reference genome (indexed)")
    parser.add_argument("--bam", required=True, help="BAM/CRAM for the truth-set sample")
    parser.add_argument("--truth-vcf", required=True, help="truth VCF with sequence-resolved TR alleles")
    parser.add_argument("--confident-bed", help="benchmark regions; loci not fully inside are skipped")
    parser.add_argument(
        "--repeat-bed",
        help="BED of repeat loci. Use when the truth VCF is a general call set containing "
        "SNVs and non-repeat variants: each BED locus becomes one target, its truth alleles "
        "are the summed length change of the variants inside it, and a locus with no variant "
        "inside the confident regions counts as homozygous reference",
    )
    parser.add_argument("--out-prefix", required=True, help="prefix for the output files")
    parser.add_argument("--modes", default=",".join(MODES), help=f"comma-separated, from {MODES}")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--regions", help="comma-separated chromosomes to restrict to (e.g. chr1)")
    parser.add_argument("--max-loci", type=int, help="stop after this many truth loci (smoke runs)")
    parser.add_argument(
        "--sample-loci",
        type=int,
        help="score a random sample of this many loci instead of the first --max-loci. "
        "Use for any partial benchmark: the first N loci of a catalog are one end of one "
        "chromosome and are not representative",
    )
    parser.add_argument("--seed", type=int, default=0, help="seed for --sample-loci")
    parser.add_argument(
        "--expect-binary-md5",
        help="abort unless the binary's md5 is exactly this. A benchmark series compares "
        "builds, and a build is identified by its content, not by its filename: pass the "
        "md5 of the build you intend to run and a stale, renamed or wrongly rebuilt binary "
        "stops the run instead of quietly producing a mislabelled result",
    )
    parser.add_argument(
        "--truth-cache",
        help="reuse the parsed truth set across runs. Building it walks the whole truth VCF "
        "and repeat catalog, which is identical for every build being compared, so pass the "
        "same path to each run of a benchmark series. Rebuilt automatically when the inputs, "
        "the sampling or the parsing change",
    )
    parser.add_argument(
        "--strdust-arg",
        action="append",
        default=[],
        help="extra argument passed through to STRdust, repeatable (e.g. --strdust-arg=--unphased)",
    )
    parser.add_argument(
        "--no-mode-flag",
        action="store_true",
        help="omit --mode when calling STRdust, for versions before v0.21.0 that lack it "
        "(their single alignment path corresponds to --modes sensitive)",
    )
    parser.add_argument("--from-vcfs", action="store_true", help="re-analyse existing output, do not rerun")
    parser.add_argument("--plot", action="store_true", help="write a called-vs-truth scatter per mode")
    args = parser.parse_args()

    modes = [m.strip() for m in args.modes.split(",") if m.strip()]
    unknown = [m for m in modes if m not in MODES]
    if unknown:
        sys.exit(f"unknown mode(s) {unknown}; expected any of {MODES}")
    if args.no_mode_flag and len(modes) > 1:
        sys.exit(
            "--no-mode-flag means the binary has only one genotyping path, so requesting "
            f"{modes} would run it {len(modes)} times identically; pass --modes sensitive"
        )

    regions_filter = set(args.regions.split(",")) if args.regions else None

    fingerprint = None if args.from_vcfs else binary_fingerprint(args.binary)
    if fingerprint:
        print(
            f"[benchmark] binary {fingerprint['path']}\n"
            f"[benchmark]   md5 {fingerprint['md5']}  ({fingerprint['size']} bytes)",
            file=sys.stderr,
        )
    if args.expect_binary_md5:
        expected = args.expect_binary_md5.strip().lower()
        if args.from_vcfs:
            sys.exit("--expect-binary-md5 cannot be checked with --from-vcfs: no binary is run")
        if fingerprint is None:
            sys.exit(f"--expect-binary-md5 given but {args.binary} could not be read")
        if fingerprint["md5"] != expected:
            sys.exit(
                "binary md5 mismatch, refusing to run:\n"
                f"  expected {expected}\n"
                f"  found    {fingerprint['md5']}  ({fingerprint['path']})\n"
                "The binary is not the build this run is meant to measure. Re-copy it, or "
                "rebuild from the intended branch, and check the tree is on that branch."
            )
        print(f"[benchmark]   md5 matches --expect-binary-md5", file=sys.stderr)

    cache_key = truth_cache_key(args, regions_filter) if args.truth_cache else None
    cached = load_truth_cache(args.truth_cache, cache_key) if args.truth_cache else None
    if cached is not None:
        truth_loci, bed, skipped = cached
        print(f"[benchmark] reusing truth set from {args.truth_cache}", file=sys.stderr)
    else:
        print("[benchmark] reading truth set", file=sys.stderr)
        confident = load_confident_regions(args.confident_bed)
        if args.repeat_bed:
            truth_loci, bed, skipped = read_truth_over_repeats(
                args.truth_vcf,
                args.repeat_bed,
                confident,
                regions_filter,
                args.max_loci,
                args.sample_loci,
                args.seed,
            )
        else:
            truth_loci, bed, skipped = read_truth(
                args.truth_vcf, confident, regions_filter, args.max_loci, args.sample_loci, args.seed
            )
        if args.truth_cache:
            save_truth_cache(args.truth_cache, cache_key, truth_loci, bed, skipped)
    if not args.confident_bed:
        print(
            "[benchmark] WARNING: no --confident-bed. Every catalog locus without a truth "
            "record counts as homozygous reference, including loci the truth set does not "
            "cover, which inflates the 0/0 class with unverified loci.",
            file=sys.stderr,
        )
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
                args.binary,
                args.fasta,
                args.bam,
                bed_path,
                mode,
                args.threads,
                args.strdust_arg,
                out_vcf,
                mode_flag=not args.no_mode_flag,
            )
            print(
                f"[benchmark] --mode {mode}: {timings[mode]['cpu_seconds']:.1f}s CPU, "
                f"{timings[mode]['wall_seconds']:.1f}s wall",
                file=sys.stderr,
            )
        all_rows.extend(compare(truth_loci, read_strdust(out_vcf), mode))

    per_locus = Path(f"{args.out_prefix}.per_locus.tsv")
    write_tsv(all_rows, per_locus)
    summary = summarise(all_rows, timings, fingerprint)
    summary_path = Path(f"{args.out_prefix}.summary.tsv")
    write_tsv(summary, summary_path)
    confusion_rows = confusion(all_rows)
    confusion_path = Path(f"{args.out_prefix}.genotype_confusion.tsv")
    write_tsv(confusion_rows, confusion_path)

    print(
        f"\n[benchmark] wrote {per_locus}, {summary_path} and {confusion_path}\n",
        file=sys.stderr,
    )
    if summary:
        columns = list(summary[0].keys())
        def cell(value):
            return "" if value is None else str(value)

        widths = [max(len(c), max(len(cell(r[c])) for r in summary)) for c in columns]
        print("  ".join(c.ljust(w) for c, w in zip(columns, widths)))
        for row in summary:
            print("  ".join(cell(row[c]).ljust(w) for c, w in zip(columns, widths)))

    if confusion_rows:
        print("\ntruth genotype class vs called:", file=sys.stderr)
        for row in confusion_rows:
            print(
                f"  {row['mode']:<10} {row['truth_gt']:>4} -> {row['called_gt']:<8}"
                f" {row['loci']:>7} loci  {row['pct_of_truth_class']:>6}%",
                file=sys.stderr,
            )

    if args.plot:
        plot(all_rows, f"{args.out_prefix}.png")


if __name__ == "__main__":
    main()
