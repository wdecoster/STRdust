"""Tests for benchmark_truth.py, the harness every accuracy number of STRdust comes from.

Each case is a small synthetic fixture; no real data is needed. Several pin down bugs that
once changed what the benchmark measured without any error (issue #38):

- the truth set attributing a whole spanning deletion to a tiny repeat (clamp_to_locus)
- reference calls going unrecognised once they carried RB=0 (read_strdust / compare)
- the truth genotype class ignoring SNVs inside the repeat (resolve_locus)

Run with ``python3 -m unittest discover -s misc`` (standard library only).
"""

import argparse
import os
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import benchmark_truth as bt  # noqa: E402


def write(directory, name, text):
    path = Path(directory) / name
    path.write_text(text)
    return path


class ClampToLocus(unittest.TestCase):
    """A deletion only removes the locus bases it covers; deleted bases start after POS."""

    LOCUS = (200, 261)  # 61 bp, 0-based half-open

    def clamp(self, var_start, deltas):
        clamped, _phased = bt.clamp_to_locus((deltas, True), var_start, *self.LOCUS)
        return clamped

    def test_deletion_spanning_the_locus_is_cut_to_the_locus(self):
        # the original bug: a 61 bp repeat got a truth allele of -145,777
        self.assertEqual(self.clamp(100, [(-145_777, 1), (0, 0)]), [(-61, 1), (0, 0)])

    def test_deletion_half_covering_the_locus(self):
        # deletes [231, 291): 30 of the locus bases
        self.assertEqual(self.clamp(230, [(-60, 1)]), [(-30, 1)])

    def test_deletion_inside_the_locus_is_untouched(self):
        self.assertEqual(self.clamp(210, [(-5, 1)]), [(-5, 1)])

    def test_deletion_abutting_the_locus_leaves_it_reference(self):
        # deletes [190, 200), ending exactly where the locus starts: only the anchor-free
        # overlap counts, so no length change and the allele is reference for the class
        self.assertEqual(self.clamp(189, [(-10, 1)]), [(0, 0)])

    def test_deletion_whose_anchor_is_the_last_locus_base(self):
        # VCF anchor at 260 (inside), deleted bases start at 261 (outside)
        self.assertEqual(self.clamp(260, [(-10, 1)]), [(0, 0)])

    def test_insertions_and_snvs_are_not_clamped(self):
        self.assertEqual(self.clamp(150, [(500, 1), (0, 2)]), [(500, 1), (0, 2)])


class ResolveLocus(unittest.TestCase):
    """Lengths are summed per haplotype; the genotype class is sequence-aware."""

    def test_no_records_is_homozygous_reference(self):
        self.assertEqual(bt.resolve_locus([]), ([0, 0], "0/0"))

    def test_heterozygous_snv_is_0_1_with_unchanged_lengths(self):
        # the length-only class called this 0/0, scoring STRdust wrong for reporting it
        self.assertEqual(bt.resolve_locus([([(0, 0), (0, 1)], True)]), ([0, 0], "0/1"))

    def test_homozygous_snv_is_1_1(self):
        self.assertEqual(bt.resolve_locus([([(0, 1), (0, 1)], True)]), ([0, 0], "1/1"))

    def test_different_alts_on_each_haplotype_are_1_2(self):
        records = [([(0, 1), (0, 0)], True), ([(0, 0), (0, 1)], True)]
        self.assertEqual(bt.resolve_locus(records), ([0, 0], "1/2"))

    def test_same_length_change_with_different_sequence_is_1_2(self):
        records = [([(3, 1), (0, 0)], True), ([(0, 0), (3, 1)], True)]
        self.assertEqual(bt.resolve_locus(records), ([3, 3], "1/2"))

    def test_phased_lengths_sum_per_haplotype(self):
        records = [([(4, 1), (0, 0)], True), ([(2, 1), (-1, 1)], True)]
        self.assertEqual(bt.resolve_locus(records), ([-1, 6], "1/2"))

    def test_single_unphased_length_change_is_resolved(self):
        self.assertEqual(bt.resolve_locus([([(0, 0), (5, 1)], False)]), ([0, 5], "0/1"))

    def test_unphased_length_change_beside_another_varying_record_is_skipped(self):
        # phased 0|10 plus unphased 5/0 fits both [5, 10] and [0, 15]
        records = [([(0, 0), (10, 1)], True), ([(5, 1), (0, 0)], False)]
        self.assertIsNone(bt.resolve_locus(records))

    def test_two_unphased_length_changes_are_skipped(self):
        records = [([(0, 0), (10, 1)], False), ([(5, 1), (0, 0)], False)]
        self.assertIsNone(bt.resolve_locus(records))

    def test_unphased_snv_beside_a_phased_length_change_leaves_the_class_unscored(self):
        # the lengths do not depend on the SNV's phase, but 0/1 vs 1/2 does: the locus is
        # kept (so --sample-loci draws the same set) with its class left out
        records = [([(0, 0), (4, 1)], True), ([(0, 1), (0, 0)], False)]
        self.assertEqual(bt.resolve_locus(records), ([0, 4], None))

    def test_homozygous_record_needs_no_phase(self):
        records = [([(2, 1), (2, 1)], False), ([(0, 0), (0, 1)], True)]
        self.assertEqual(bt.resolve_locus(records), ([2, 2], "1/2"))

    def test_haploid(self):
        self.assertEqual(bt.resolve_locus([([(0, 1)], False)]), ([0], "1"))
        self.assertEqual(bt.resolve_locus([([(0, 0)], False)]), ([0], "0"))


class ReadTruthOverRepeats(unittest.TestCase):
    """End to end over a tiny VCF and catalog, through clamping and resolution."""

    def test_snv_spanning_deletion_and_invariant_locus(self):
        with tempfile.TemporaryDirectory() as tmp:
            vcf = write(
                tmp,
                "truth.vcf",
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
                # SNV inside locus 1 ([100, 120))
                "chr1\t111\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\n"
                # 1000 bp deletion swallowing locus 2 ([300, 320)) whole
                "chr1\t201\t.\t" + "A" * 1001 + "\tA\t.\tPASS\t.\tGT\t1|0\n",
            )
            bed = write(tmp, "repeats.bed", "chr1\t100\t120\nchr1\t300\t320\nchr1\t5000\t5010\n")
            loci, bed_out, skipped = bt.read_truth_over_repeats(
                str(vcf), str(bed), None, None, None
            )
        self.assertEqual(loci[("chr1", 101)], {"ref_len": 20, "lengths": [0, 0], "gt_class": "0/1"})
        self.assertEqual(loci[("chr1", 301)]["lengths"], [-20, 0])
        self.assertEqual(loci[("chr1", 5001)]["gt_class"], "0/0")
        self.assertEqual(len(bed_out), 3)


class ReadStrdust(unittest.TestCase):
    """Reference calls are recognised by their flag, whether or not RB is present."""

    HEADER = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
    )

    def parse(self, records):
        with tempfile.TemporaryDirectory() as tmp:
            return bt.read_strdust(str(write(tmp, "calls.vcf", self.HEADER + records)))

    def test_quickref_with_rb_zero_is_a_reference_call(self):
        # the second bug: once reference records carried RB=0 they were labelled "called"
        calls = self.parse(
            "chr1\t101\t.\tACACAC\t.\t.\t.\tEND=106;QUICKREF\tGT:RB:SUP\t0|0:0,0:.\n"
        )
        rows = bt.compare({("chr1", 101): {"ref_len": 6, "lengths": [0, 0], "gt_class": "0/0"}}, calls, "fast")
        self.assertTrue(calls[("chr1", 101)]["quickref"])
        self.assertEqual(rows[0]["status"], "quickref_reference")

    def test_quickref_without_rb_is_a_reference_call(self):
        calls = self.parse("chr1\t101\t.\tACACAC\t.\t.\t.\tEND=106;QUICKREF\tGT\t0|0\n")
        rows = bt.compare({("chr1", 101): {"ref_len": 6, "lengths": [0, 0], "gt_class": "0/0"}}, calls, "fast")
        self.assertEqual(rows[0]["status"], "quickref_reference")
        self.assertEqual(rows[0]["max_abs_diff"], 0)

    def test_genotyped_reference_and_variant_calls(self):
        calls = self.parse(
            "chr1\t101\t.\tACACAC\t.\t.\t.\tEND=106\tGT:RB:SUP\t0|0:0,0:7,6\n"
            "chr1\t201\t.\tACACAC\tACACACAC\t.\t.\tEND=206\tGT:RB:SUP\t0|1:0,2:7,6\n"
            "chr1\t301\t.\tACACAC\t.\t.\t.\tEND=306\tGT:RB:SUP\t.|.:.,.:1,1\n"
        )
        self.assertEqual(calls[("chr1", 101)]["lengths"], [0, 0])
        self.assertFalse(calls[("chr1", 101)]["quickref"])
        self.assertEqual(calls[("chr1", 201)]["gt_class"], "0/1")
        self.assertIsNone(calls[("chr1", 301)]["lengths"])
        truth = {
            ("chr1", 101): {"ref_len": 6, "lengths": [0, 0], "gt_class": "0/0"},
            ("chr1", 201): {"ref_len": 6, "lengths": [0, 2], "gt_class": "0/1"},
            ("chr1", 301): {"ref_len": 6, "lengths": [0, 0], "gt_class": "0/0"},
            ("chr1", 401): {"ref_len": 6, "lengths": [0, 0], "gt_class": "0/0"},
        }
        status = {row["pos"]: row["status"] for row in bt.compare(truth, calls, "fast")}
        self.assertEqual(status, {101: "called", 201: "called", 301: "no_call", 401: "not_reported"})


class InsideConfident(unittest.TestCase):
    """The bisect over running max ends agrees with a plain scan, nested regions included."""

    def test_matches_a_linear_scan(self):
        rng = random.Random(38)
        spans = []
        for _ in range(200):
            start = rng.randrange(0, 10_000)
            spans.append((start, start + rng.choice([5, 50, 500, 3000])))
        with tempfile.TemporaryDirectory() as tmp:
            bed = write(tmp, "conf.bed", "".join(f"chr1\t{s}\t{e}\n" for s, e in spans))
            regions = bt.load_confident_regions(str(bed))
        for _ in range(5000):
            start = rng.randrange(0, 12_000)
            end = start + rng.randrange(1, 400)
            expected = any(s <= start and end <= e for s, e in spans)
            self.assertEqual(bt.inside_confident(regions, "chr1", start, end), expected, (start, end))

    def test_edges_and_other_chromosomes(self):
        with tempfile.TemporaryDirectory() as tmp:
            regions = bt.load_confident_regions(str(write(tmp, "c.bed", "chr1\t100\t200\n")))
        self.assertTrue(bt.inside_confident(regions, "chr1", 100, 200))
        self.assertFalse(bt.inside_confident(regions, "chr1", 99, 150))
        self.assertFalse(bt.inside_confident(regions, "chr1", 150, 201))
        self.assertFalse(bt.inside_confident(regions, "chr2", 120, 130))
        self.assertTrue(bt.inside_confident(None, "chr2", 120, 130))


class TruthCache(unittest.TestCase):
    """The cache round-trips, and is rebuilt when any input that defines the truth changes."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = Path(self.tmp.name)
        self.vcf = write(self.dir, "truth.vcf", "#header\n")
        self.args = argparse.Namespace(
            truth_vcf=str(self.vcf), repeat_bed=None, confident_bed=None,
            max_loci=None, sample_loci=1000, seed=5,
        )
        self.cache = self.dir / "truth.json"
        self.loci = {
            ("chr1", 101): {"ref_len": 20, "lengths": [0, 3], "gt_class": "0/1"},
            ("chr1", 501): {"ref_len": 8, "lengths": [0, 0], "gt_class": None},
        }
        self.bed = [("chr1", 100, 120), ("chr1", 500, 508)]

    def tearDown(self):
        self.tmp.cleanup()

    def test_round_trip(self):
        key = bt.truth_cache_key(self.args, None)
        bt.save_truth_cache(self.cache, key, self.loci, self.bed, {"outside_confident": 2})
        self.assertEqual(
            bt.load_truth_cache(self.cache, bt.truth_cache_key(self.args, None)),
            (self.loci, self.bed, {"outside_confident": 2}),
        )

    def test_stale_when_the_truth_vcf_changes(self):
        bt.save_truth_cache(self.cache, bt.truth_cache_key(self.args, None), self.loci, self.bed, {})
        self.vcf.write_text("#header\nchr1\t1\t.\tA\tG\t.\t.\t.\tGT\t0|1\n")
        self.assertIsNone(bt.load_truth_cache(self.cache, bt.truth_cache_key(self.args, None)))

    def test_stale_when_sampling_or_version_changes(self):
        bt.save_truth_cache(self.cache, bt.truth_cache_key(self.args, None), self.loci, self.bed, {})
        self.args.seed = 6
        self.assertIsNone(bt.load_truth_cache(self.cache, bt.truth_cache_key(self.args, None)))
        self.args.seed = 5
        key = bt.truth_cache_key(self.args, None)
        key["version"] -= 1
        self.assertIsNone(bt.load_truth_cache(self.cache, key))

    def test_absent_or_unreadable_cache(self):
        self.assertIsNone(bt.load_truth_cache(self.cache, {}))
        self.cache.write_text("{not json")
        self.assertIsNone(bt.load_truth_cache(self.cache, {}))


if __name__ == "__main__":
    unittest.main()
