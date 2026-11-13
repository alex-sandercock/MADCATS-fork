import argparse
import csv
import dataclasses
import gzip
import json
import random
import tempfile
import unittest
from pathlib import Path

import madcats as mp


def rule(min_global=1, min_quality=20, clean_reads=0, rep_reads=10, rep_samples=0, depth_reads=0, freq=0.0):
    return mp.EvidenceRule(
        min_global_count=min_global,
        discovery_min_quality=min_quality,
        quality_min_clean_reads=clean_reads,
        reproducible_min_reads=rep_reads,
        reproducible_min_samples=rep_samples,
        depth_min_reads=depth_reads,
        min_supporting_freq=freq,
    )


class MadcPipelineTests(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(mp.reverse_complement("ACGTRY"), "RYACGT")

    def test_iupac_compatibility(self):
        self.assertEqual(mp.compatible_mismatches("ACGT", "AYGN"), 0)
        self.assertEqual(mp.compatible_mismatches("ACGT", "ARGN"), 1)
        self.assertTrue(mp.bases_compatible("R", "A"))
        self.assertFalse(mp.bases_compatible("R", "C"))

    def test_iupac_aware_classification(self):
        locus = mp.Locus(0, "m1", "AAAA", "", "R", "C", 5, "AAAAA" + "R" + "AAAA", "AAAAA" + "C" + "AAAA")
        self.assertEqual(mp.classify_sequence(locus, "AAAAA" + "A" + "AAAA"), "RefMatch")
        self.assertEqual(mp.classify_sequence(locus, "AAAAA" + "G" + "AAAA"), "RefMatch")
        self.assertEqual(mp.classify_sequence(locus, "AAAAA" + "C" + "AAAA"), "Alt")
        self.assertEqual(mp.classify_sequence(locus, "AAAAA" + "T" + "AAAA"), "Other")
        # An undetermined observed base is missing data, not a Ref/Alt match.
        self.assertEqual(mp.classify_sequence(locus, "AAAAA" + "N" + "AAAA"), "Other")
        concrete = mp.Locus(0, "m2", "AAAA", "", "A", "C", 5, "AAAAA" + "A" + "AAAA", "AAAAA" + "C" + "AAAA")
        self.assertEqual(mp.classify_sequence(concrete, "AAAAA" + "N" + "AAAA"), "Other")

    def test_evidence_rule_routes(self):
        combined = rule(clean_reads=5, rep_samples=2, depth_reads=30)
        self.assertEqual(combined.evaluate(0, 0, 0, 0), (False, "below_min_global_count"))
        self.assertEqual(combined.evaluate(20, 5, 1, 20), (True, "quality"))
        self.assertEqual(combined.evaluate(20, 0, 2, 10), (True, "reproducibility"))
        self.assertEqual(combined.evaluate(40, 0, 1, 30), (True, "depth"))
        self.assertEqual(combined.evaluate(40, 5, 2, 30), (True, "reproducibility+depth+quality"))
        self.assertEqual(combined.evaluate(20, 4, 1, 29), (False, "failed_evidence"))
        exhaustive = rule(clean_reads=0, rep_samples=0, depth_reads=0)
        self.assertEqual(exhaustive.evaluate(1, 0, 0, 1), (True, "no_route_required"))

    def test_denoise_absorbs_low_abundance_neighbor_and_protects_alt(self):
        ref = "A" * 10
        error = "A" * 9 + "C"
        alt = "C" + "A" * 9
        centroids, mapping = mp.denoise_counts(
            {ref: 100, error: 3, alt: 1}, max_distance=1, alpha=2.0, protected=(ref, alt)
        )
        self.assertEqual(centroids[ref], 103)
        self.assertEqual(mapping[error], ref)
        self.assertIn(alt, centroids)

    def test_quality_supported_sequence_is_protected_before_denoising(self):
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        novel = "A" * 10 + "C" + "A" * 9
        locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
        sample = mp.Sample("s1", "s1", "unused.fastq", {})
        result = mp.SampleResult(
            sample=sample,
            counts={(0, ref): 100, (0, novel): 5},
            quality_support={(0, ref): (100, 100, 100, 100), (0, novel): (5, 5, 5, 5)},
            partial_counts={},
            partial_quality_support={},
            qc={},
        )
        finalized = mp.finalize_counts(
            [result],
            [locus],
            mp.DenoiseSettings(sample_distance=1, global_distance=1, alpha=2.0),
            rule=rule(clean_reads=5),
        )
        self.assertEqual(finalized[0].counts[(0, novel)], 5)

    def test_parallel_finalize_matches_the_sequential_passes(self):
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        ref2 = "C" * 20
        alt2 = "C" * 9 + "T" + "C" * 10
        # A neighbour of Alt private to one sample: 27 reads against 80 Alt reads
        # is below the per-sample skew, so it survives that pass, while the
        # project-wide 400 Alt reads absorb it in the pooled pass (it is not
        # reproducible, so it is not protected). The random Ref neighbours
        # exercise the per-sample pass; this exercises the pooled one.
        stray = alt[:15] + "C" + alt[16:]
        loci = [mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt), mp.Locus(1, "m2", "CCCCCC", "", "C", "T", 9, ref2, alt2)]
        rng = random.Random(7)
        results = []
        for sample_id in ("s1", "s2", "s3", "s4", "s5"):
            sample = mp.Sample(sample_id, sample_id, "unused.fastq", {})
            counts = {(0, ref): 200, (0, alt): 80, (1, ref2): 150, (1, alt2): 60}
            for _ in range(12):
                pos = rng.randrange(20)
                noisy = ref[:pos] + rng.choice("CGT") + ref[pos + 1:]
                counts[(0, noisy)] = counts.get((0, noisy), 0) + rng.randrange(1, 6)
                noisy2 = ref2[:pos] + rng.choice("AGT") + ref2[pos + 1:]
                counts[(1, noisy2)] = counts.get((1, noisy2), 0) + rng.randrange(1, 6)
            if sample_id == "s5":
                counts[(0, stray)] = 27
            results.append(
                mp.SampleResult(
                    sample=sample,
                    counts=counts,
                    quality_support={key: (c, c, c, c) for key, c in counts.items()},
                    partial_counts={},
                    partial_quality_support={},
                    qc={},
                )
            )
        settings = mp.DenoiseSettings(sample_distance=1, global_distance=1, alpha=2.0)
        sequential = mp.finalize_counts(results, loci, settings, rule(rep_reads=5, rep_samples=2), jobs=1)
        parallel = mp.finalize_counts(results, loci, settings, rule(rep_reads=5, rep_samples=2), jobs=2)
        # Both passes moved reads in the sequential run: the per-sample pass in
        # every sample, the pooled pass only for the private neighbour in s5.
        self.assertTrue(all(r.qc["reads_absorbed_sample_denoise"] > 0 for r in sequential))
        self.assertEqual([r.qc["reads_absorbed_global_denoise"] for r in sequential], [0, 0, 0, 0, 27])
        self.assertNotIn((0, stray), sequential[4].counts)
        self.assertEqual([r.counts for r in parallel], [r.counts for r in sequential])
        self.assertEqual([r.quality_support for r in parallel], [r.quality_support for r in sequential])
        self.assertEqual([r.qc for r in parallel], [r.qc for r in sequential])

    def test_reproducible_low_quality_sequence_is_protected_before_denoising(self):
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        novel = "A" * 10 + "C" + "A" * 9
        locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
        results = []
        for sample_id in ("s1", "s2"):
            sample = mp.Sample(sample_id, sample_id, "unused.fastq", {})
            results.append(
                mp.SampleResult(
                    sample=sample,
                    counts={(0, ref): 100, (0, novel): 5},
                    # One tag base is systematically low quality: zero clean reads.
                    quality_support={(0, ref): (100, 100, 100, 100), (0, novel): (0, 0, 0, 0)},
                    partial_counts={},
                    partial_quality_support={},
                    qc={},
                )
            )
        evidence = rule(rep_reads=5, rep_samples=2)
        finalized = mp.finalize_counts(
            results,
            [locus],
            mp.DenoiseSettings(sample_distance=1, global_distance=1, alpha=2.0),
            rule=evidence,
        )
        for result in finalized:
            self.assertEqual(result.counts[(0, novel)], 5)
        # Without the reproducibility route the same sequence is absorbed.
        absorbed = mp.finalize_counts(
            results,
            [locus],
            mp.DenoiseSettings(sample_distance=1, global_distance=1, alpha=2.0),
            rule=rule(clean_reads=5),
        )
        for result in absorbed:
            self.assertNotIn((0, novel), result.counts)
            self.assertEqual(result.qc["reads_absorbed_sample_denoise"], 5)

    def test_partial_assign_never_invents_an_allele(self):
        # A prefix of nothing learned is unmapped; a prefix shared by Ref and
        # Alt is ambiguous; neither becomes a row.
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        novel = "A" * 10 + "C" + "A" * 9
        novel_prefix = novel[:15]
        ambiguous_prefix = ref[:6]
        locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
        sample = mp.Sample("s1", "s1", "unused.fastq", {})
        result = mp.SampleResult(
            sample=sample,
            counts={(0, ref): 10, (0, alt): 10},
            quality_support={(0, ref): (10, 10, 10, 10), (0, alt): (10, 10, 10, 10)},
            partial_counts={(0, novel_prefix): 5, (0, ambiguous_prefix): 2},
            partial_quality_support={
                (0, novel_prefix): (5, 5, 5, 5),
                (0, ambiguous_prefix): (2, 2, 2, 2),
            },
            qc={},
        )
        updated, summary = mp.apply_partial_reads(
            [result],
            [locus],
            mode="on",
            rule=rule(clean_reads=5),
        )
        self.assertEqual(summary["reads_unmapped"], 5)
        self.assertEqual(summary["reads_ambiguous"], 2)
        self.assertNotIn("extended_alleles", summary)
        self.assertEqual(set(updated[0].counts), {(0, ref), (0, alt)})

    @staticmethod
    def _partial_locus(with_refmatch=True, refmatch_positions=(106, 107, 108)):
        """A 109-bp locus: Ref, Alt (target at index 60) and optionally a RefMatch."""
        random.seed(7)
        ref = "".join(random.choice("ACGT") for _ in range(109))
        alt = ref[:60] + ("A" if ref[60] != "A" else "C") + ref[61:]
        locus = mp.Locus(0, "m1", ref[:24], "", ref[60], alt[60], 60, ref, alt)
        counts = {(0, ref): 10, (0, alt): 10}
        refmatch = None
        if with_refmatch:
            bases = list(ref)
            for position in refmatch_positions:
                bases[position] = "A" if bases[position] != "A" else "C"
            refmatch = "".join(bases)
            counts[(0, refmatch)] = 10
        quality = {key: (10, 10, 10, 10) for key in counts}
        return locus, ref, alt, refmatch, counts, quality

    def _apply_partials(self, locus, counts, quality, partials):
        sample = mp.Sample("s1", "s1", "unused.fastq", {})
        result = mp.SampleResult(
            sample=sample,
            counts=dict(counts),
            quality_support=dict(quality),
            partial_counts={(0, partial): n for partial, n in partials.items()},
            partial_quality_support={(0, partial): (n, n, n, n) for partial, n in partials.items()},
            qc={},
        )
        return mp.apply_partial_reads([result], [locus], mode="on", rule=rule(clean_reads=5))

    def test_partial_assign_locus_guard_drops_both_alleles_when_a_tail_variant_exists(self):
        # Ref, Alt and a RefMatch differing from Ref only at positions 107-109:
        # 106-bp partials of Ref and of Alt are both dropped, so the Ref:Alt
        # ratio cannot move even though the Alt partial alone has a unique prefix.
        locus, ref, alt, refmatch, counts, quality = self._partial_locus(with_refmatch=True)
        updated, summary = self._apply_partials(locus, counts, quality, {ref[:106]: 7, alt[:106]: 5})
        self.assertEqual(summary["reads_ambiguous"], 12)
        self.assertNotIn("reads_assigned", summary)
        self.assertEqual(updated[0].counts[(0, ref)], 10)
        self.assertEqual(updated[0].counts[(0, alt)], 10)
        self.assertEqual(updated[0].qc["partial_reads_ambiguous"], 12)

    def test_partial_assign_credits_both_alleles_without_a_tail_variant(self):
        locus, ref, alt, _, counts, quality = self._partial_locus(with_refmatch=False)
        updated, summary = self._apply_partials(locus, counts, quality, {ref[:106]: 7, alt[:106]: 5})
        self.assertEqual(summary["reads_assigned"], 12)
        self.assertEqual(updated[0].counts[(0, ref)], 17)
        self.assertEqual(updated[0].counts[(0, alt)], 15)
        self.assertEqual(set(updated[0].counts), {(0, ref), (0, alt)})

    def test_partial_not_reaching_the_target_is_always_ambiguous(self):
        locus, ref, alt, _, counts, quality = self._partial_locus(with_refmatch=False)
        # 60 bases stop just before the target (index 60); 61 bases include it.
        updated, summary = self._apply_partials(locus, counts, quality, {ref[:60]: 4, alt[:61]: 3})
        self.assertEqual(summary["reads_ambiguous"], 4)
        self.assertEqual(summary["reads_assigned"], 3)
        self.assertEqual(updated[0].counts[(0, alt)], 13)

    def test_partial_shorter_than_locus_shared_prefix_is_never_credited(self):
        # The RefMatch differs from Ref at index 80 only, so the locus's longest
        # shared prefix is 80.  A 70-bp Alt partial reaches the target and is a
        # prefix of Alt alone, yet it is still dropped by the locus-level guard.
        locus, ref, alt, refmatch, counts, quality = self._partial_locus(with_refmatch=True, refmatch_positions=(80,))
        self.assertEqual(mp._longest_shared_prefix([ref, alt, refmatch]), 80)
        updated, summary = self._apply_partials(locus, counts, quality, {alt[:70]: 6, alt[:81]: 2})
        self.assertEqual(summary["reads_ambiguous"], 6)
        self.assertEqual(summary["reads_assigned"], 2)
        self.assertEqual(updated[0].counts[(0, alt)], 12)

    def test_partial_defaults_resolve_from_tag_length(self):
        parser = mp.build_parser()
        args = parser.parse_args(["generate", "--fastq-dir", "x", "--reference-tags", "y", "--tag-length", "109", "--output", "z"])
        self.assertEqual(args.partial_read_mode, "on")
        legacy = parser.parse_args(["generate", "--fastq-dir", "x", "--reference-tags", "y", "--output", "z", "--partial-read-mode", "assign"])
        self.assertEqual(legacy.partial_read_mode, "on")
        self.assertIsNone(args.min_partial_length)
        self.assertEqual(max(1, args.tag_length - mp.PARTIAL_TAIL_DEFAULT), 100)

    def test_one_reading_for_category_and_variants(self):
        # A deletion inside the last 12 bases is not a verified reading: both the
        # category (fixed offset) and the Variants (substitutions) ignore it.
        random.seed(11)
        ref = "".join(random.choice("ACGT") for _ in range(109))
        alt = ref[:60] + ("A" if ref[60] != "A" else "C") + ref[61:]
        locus = mp.Locus(0, "m1", ref[:24], "", ref[60], alt[60], 60, ref, alt)
        tail_deletion = ref[:100] + ref[102:] + "GT"
        self.assertEqual(mp.classify_sequence(locus, tail_deletion, 20), "RefMatch")
        events = mp.variant_events(locus, tail_deletion, 20)
        self.assertTrue(events and all(kind == "snp" for kind, *_ in events))
        self.assertEqual(mp.marker_type_label(locus, tail_deletion, 20, events), "snp")
        # A verified upstream deletion shifts the target for both consumers alike.
        shifted = ref[:30] + ref[31:60] + alt[60] + ref[61:] + "G"
        self.assertEqual(mp.classify_sequence(locus, shifted, 20), "AltMatch")
        events = mp.variant_events(locus, shifted, 20)
        kinds = [kind for kind, *_ in events]
        self.assertIn("del", kinds)
        self.assertIn(("snp", 61, ref[60], alt[60]), events)
        self.assertEqual(mp.marker_type_label(locus, shifted, 20, events), "indel")

    def test_category_reads_the_base_where_the_variants_alignment_puts_it(self):
        # Invariant over random indel rows: classify_sequence reads the target at
        # the read index that the (shared) accepted reading maps it to.
        random.seed(5)
        checked = 0
        for _ in range(300):
            ref = "".join(random.choice("ACGT") for _ in range(109))
            alt = ref[:60] + ("A" if ref[60] != "A" else "C") + ref[61:]
            locus = mp.Locus(0, "m1", ref[:24], "", ref[60], alt[60], 60, ref, alt)
            position = random.randrange(24, 90)
            size = random.randrange(1, 6)
            if random.random() < 0.5:
                row = ref[:position] + ref[position + size :] + "".join(random.choice("ACGT") for _ in range(size))
            else:
                row = (ref[:position] + "".join(random.choice("ACGT") for _ in range(size)) + ref[position:])[:109]
            if row in (ref, alt):
                continue
            reading = mp._reading(row, ref, 20)
            if reading is not None:
                mapped = reading.ref_to_read.get(60)
                expected = "Other" if mapped is None else (
                    "RefMatch" if row[mapped] == ref[60] else "AltMatch" if row[mapped] == alt[60] else "Other")
            else:
                expected = "RefMatch" if row[60] == ref[60] else "AltMatch" if row[60] == alt[60] else "Other"
            self.assertEqual(mp.classify_sequence(locus, row, 20), expected, (row, reading is not None))
            checked += 1
        self.assertGreater(checked, 250)

    def test_designed_substitution_does_not_count_against_the_indel_guard(self):
        # AltMatch carrying the designed Alt base, a verified 1-bp deletion and
        # six further substitutions: against Ref that is seven substitutions,
        # against Alt six, so the reading is accepted and both consumers see
        # the deletion.
        random.seed(13)
        ref = "".join(random.choice("ACGT") for _ in range(109))
        alt = ref[:60] + ("A" if ref[60] != "A" else "C") + ref[61:]
        locus = mp.Locus(0, "m1", ref[:24], "", ref[60], alt[60], 60, ref, alt)
        bases = list(alt)
        for position in (26, 33, 40, 47, 70, 80):
            bases[position] = "A" if bases[position] != "A" else "C"
        row = "".join(bases[:30] + bases[31:]) + "G"
        self.assertEqual(mp.classify_sequence(locus, row, 20), "AltMatch")
        events = mp.variant_events(locus, row, 20)
        self.assertIn("del", [kind for kind, *_ in events])
        self.assertEqual(mp.marker_type_label(locus, row, 20, events), "both")

    def test_guard_judges_the_fewest_substitution_explanation(self):
        # Fewest-openings alignment: one gap plus eight substitutions (fails the
        # guard); fewest-substitution alignment of the same edit count uses
        # more gaps and fewer substitutions.  The guard must ask the latter.
        ref = "ACGTACGGTCCAGTTGACCTGATCCGTAGCATGCATGCCTAGGATCCATGACCAGTCTGACTGGACCATCAGGTCCAGATCGGACTACCGTAGCCTTAGGACTGCACGT"
        self.assertEqual(len(ref), 109)
        alt = ref[:60] + ("A" if ref[60] != "A" else "C") + ref[61:]
        locus = mp.Locus(0, "m1", ref[:24], "", ref[60], alt[60], 60, ref, alt)
        # delete one base at 30, then flip 6 bases spread downstream, then delete 2 at 90
        bases = list(ref[:30] + ref[31:])
        for position in (35, 42, 50, 57, 66, 75):
            bases[position] = "A" if bases[position] != "A" else "C"
        row = ("".join(bases[:89] + bases[91:]) + "GTC")[:109]
        fewest_openings = mp.gapped_alignment(row, ref, 20)
        fewest_subs = mp.fewest_substitution_alignment(row, ref, 20)
        self.assertEqual(fewest_openings[3], fewest_subs[3])  # same edit count
        subs = lambda aln: sum(1 for k, i, j in aln[0] if k == "M" and row[i] != ref[j])
        self.assertLessEqual(subs(fewest_subs), subs(fewest_openings))
        self.assertLessEqual(subs(fewest_subs), mp.GAP_MAX_SUBSTITUTIONS)
        self.assertIsNotNone(mp.row_reading(locus, row, 20))

    def test_write_madc_evidence_routes_and_accounting(self):
        with tempfile.TemporaryDirectory() as temporary:
            ref = "A" * 20
            alt = "A" * 6 + "G" + "A" * 13
            reproducible_novel = ref[:10] + "C" + ref[11:]
            quality_novel = ref[:11] + "C" + ref[12:]
            weak_novel = ref[:12] + "C" + ref[13:]
            locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
            results = []
            for sample_id, plan in (
                # (sequence -> (count, clean_count))
                ("s1", {reproducible_novel: (12, 0), quality_novel: (30, 30), weak_novel: (3, 0)}),
                ("s2", {reproducible_novel: (10, 0), weak_novel: (2, 0)}),
            ):
                sample = mp.Sample(sample_id, sample_id, "unused.fastq", {})
                counts = {(0, sequence): count for sequence, (count, _) in plan.items()}
                quality = {
                    (0, sequence): (clean, clean, clean, clean)
                    for sequence, (_, clean) in plan.items()
                }
                results.append(mp.SampleResult(sample, counts, quality, {}, {}, {}))
            output = Path(temporary) / "evidence.csv"
            summary = mp.write_madc(
                output,
                [locus],
                results,
                rule=rule(clean_reads=5, rep_reads=10, rep_samples=2),
            )
            sequences = {row["sequence"] for row in mp.read_madc(output)["records"]}
            self.assertIn(ref, sequences)
            self.assertIn(alt, sequences)
            self.assertIn(reproducible_novel, sequences)
            self.assertIn(quality_novel, sequences)
            self.assertNotIn(weak_novel, sequences)
            self.assertEqual(
                summary["novel_rows_by_evidence_route"],
                {"reproducibility": 1, "quality": 1},
            )
            self.assertEqual(summary["excluded_novel"]["failed_evidence"], {"rows": 1, "reads": 5})
            table_reads = sum(sum(result.counts.values()) for result in results)
            self.assertEqual(
                table_reads,
                summary["reads_in_emitted_rows"]
                + summary["excluded_novel"]["failed_evidence"]["reads"]
                + summary["excluded_novel"]["below_min_global_count"]["reads"],
            )
            self.assertEqual(summary["per_sample_reads_in_emitted_rows"], [42, 10])

    def test_depth_route_discovers_single_sample_private_allele(self):
        with tempfile.TemporaryDirectory() as temporary:
            ref = "A" * 20
            alt = "A" * 6 + "G" + "A" * 13
            private_novel = ref[:10] + "C" + ref[11:]
            locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
            sample = mp.Sample("s1", "s1", "unused.fastq", {})
            results = [
                mp.SampleResult(
                    sample,
                    {(0, private_novel): 35},
                    {(0, private_novel): (0, 0, 0, 0)},
                    {},
                    {},
                    {},
                ),
                mp.SampleResult(
                    mp.Sample("s2", "s2", "unused.fastq", {}), {}, {}, {}, {}, {}
                ),
            ]
            with_depth = mp.write_madc(
                Path(temporary) / "with_depth.csv",
                [locus],
                results,
                rule=rule(rep_samples=2, depth_reads=30),
            )
            self.assertEqual(with_depth["novel_rows_by_evidence_route"], {"depth": 1})
            without_depth = mp.write_madc(
                Path(temporary) / "without_depth.csv",
                [locus],
                results,
                rule=rule(rep_samples=2, depth_reads=0),
            )
            self.assertEqual(without_depth["excluded_novel"]["failed_evidence"]["rows"], 1)

    def test_prefix_nested_primer_cannot_shadow_locus(self):
        ref_a = "AAAACCGGTTGG"
        alt_a = ref_a[:7] + "T" + ref_a[8:]
        ref_b = "AAAACCGGAAAA"
        alt_b = ref_b[:8] + "C" + ref_b[9:]
        locus_a = mp.Locus(0, "a", "AAAACC", "", "G", "T", 7, ref_a, alt_a)
        locus_b = mp.Locus(1, "b", "AAAACCGG", "", "A", "C", 8, ref_b, alt_b)
        index = mp.PrimerIndex([locus_a, locus_b])
        candidates = index.match(ref_a)
        self.assertEqual(sorted(candidates), [0, 1])
        marker, status = mp.select_marker(ref_a, candidates, [locus_a, locus_b], max_mismatches=2)
        self.assertEqual((marker, status), (0, "assigned"))
        marker, status = mp.select_marker(ref_b, index.match(ref_b), [locus_a, locus_b], max_mismatches=2)
        self.assertEqual((marker, status), (1, "assigned"))

    def test_build_panel_rejects_variant_inside_primer(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RightPrimerSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                writer.writerow(["m1", "AACCGG", "TTTT", "C", "T", 3])
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", "AA[C/T]CGG" + "A" * 40])
            with self.assertRaisesRegex(mp.MadcError, "inside its left primer"):
                mp.build_panel(primers, designs, 20)

    def test_build_panel_rejects_alt_absent_from_design(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RightPrimerSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                # Design bracket is [A/G] but the primer table claims Alt=C.
                writer.writerow(["m1", "AACCGG", "TTTT", "A", "C", 7])
            target = "TTTTAACCGGA" + "C" * 80 + "TTTT"
            target = target[:10] + "[A/G]" + target[11:]
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            with self.assertRaisesRegex(mp.MadcError, "probe-design bracket expansion"):
                mp.build_panel(primers, designs, 20)

    def test_reference_tags_validation_errors(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "tags.csv"

            def write_tags(position):
                with path.open("w", newline="") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(
                        ["MarkerName", "LeftPrimerSequence", "RefSequence", "AltSequence", "VariantPositionInTagSequence"]
                    )
                    writer.writerow(["m1", "AAAA", "AAAATTTT", "AAAATTTG", position])

            write_tags("abc")
            with self.assertRaisesRegex(mp.MadcError, "Invalid VariantPositionInTagSequence"):
                mp.load_reference_tags(path)
            write_tags("99")
            with self.assertRaisesRegex(mp.MadcError, "outside its 8-bp tag"):
                mp.load_reference_tags(path)
            # Position 2 is inside the primer, but the position-must-differ
            # invariant catches it first (both tags carry the same primer).
            write_tags("2")
            with self.assertRaisesRegex(mp.MadcError, "do not differ"):
                mp.load_reference_tags(path)
            write_tags("8")
            loci = mp.load_reference_tags(path)
            self.assertEqual(loci[0].variant_pos, 7)
            # A declared position where Ref and Alt agree must be rejected.
            write_tags("7")
            with self.assertRaisesRegex(mp.MadcError, "do not differ"):
                mp.load_reference_tags(path)
            # Identical Ref and Alt sequences must be rejected even with an
            # explicit position.
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["MarkerName", "LeftPrimerSequence", "RefSequence", "AltSequence", "VariantPositionInTagSequence"]
                )
                writer.writerow(["m1", "AAAA", "AAAATTTT", "AAAATTTT", "8"])
            with self.assertRaisesRegex(mp.MadcError, "identical"):
                mp.load_reference_tags(path)
            # Explicit identical VariantRef/VariantAlt must be rejected.
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RefSequence",
                        "AltSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                writer.writerow(["m1", "AAAA", "AAAATTTT", "AAAATTTG", "T", "T", "8"])
            with self.assertRaisesRegex(mp.MadcError, "identical"):
                mp.load_reference_tags(path)

    def test_build_panel_rejects_identical_alleles(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RightPrimerSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                writer.writerow(["m1", "AACCGG", "TTTT", "A", "A", 7])
            target = "TTTTAACCGGA" + "C" * 80 + "TTTT"
            target = target[:10] + "[A/G]" + target[11:]
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            with self.assertRaisesRegex(mp.MadcError, "identical"):
                mp.build_panel(primers, designs, 20)

    def test_fuzzy_primer_match_returns_all_candidates_within_budget(self):
        ref_a = "AAATCCCCCCCC"
        alt_a = ref_a[:8] + "G" + ref_a[9:]
        locus_a = mp.Locus(0, "a", "AAAT", "", "C", "G", 8, ref_a, alt_a)
        ref_b = "AAAAGGGGTTTT"
        alt_b = ref_b[:8] + "A" + ref_b[9:]
        locus_b = mp.Locus(1, "b", "AAAAGGGG", "", "T", "A", 8, ref_b, alt_b)
        index = mp.PrimerIndex([locus_a, locus_b])
        # A read from locus b with one error inside b's primer that makes the
        # prefix exactly match a's short primer.
        read = "AAAT" + ref_b[4:]
        # An exact short primer must not suppress the longer one-mismatch
        # primer; full-tag identity resolves the locus.
        self.assertEqual(index.match(read, max_mismatches=1), [0, 1])
        marker, status = mp.select_marker(read, [0, 1], [locus_a, locus_b], max_mismatches=2)
        self.assertEqual((marker, status), (1, "assigned"))

    def test_min_supporting_freq_gates_dilute_sequences(self):
        with tempfile.TemporaryDirectory() as temporary:
            ref = "A" * 20
            alt = "A" * 6 + "G" + "A" * 13
            dilute = ref[:10] + "C" + ref[11:]
            private_deep = ref[:11] + "C" + ref[12:]
            locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
            results = []
            for sample_id, dilute_count, private_count in (("s1", 40, 100), ("s2", 40, 0)):
                counts = {(0, ref): 10000, (0, dilute): dilute_count}
                if private_count:
                    counts[(0, private_deep)] = private_count
                quality = {key: (0, 0, 0, 0) for key in counts}
                results.append(
                    mp.SampleResult(
                        mp.Sample(sample_id, sample_id, "unused.fastq", {}),
                        counts,
                        quality,
                        {},
                        {},
                        {},
                    )
                )
            without_floor = mp.write_madc(
                Path(temporary) / "no_floor.csv",
                [locus],
                results,
                rule=rule(rep_samples=2, depth_reads=30),
            )
            self.assertEqual(
                without_floor["novel_rows_by_evidence_route"],
                {"reproducibility+depth": 1, "depth": 1},
            )
            with_floor = mp.write_madc(
                Path(temporary) / "floor.csv",
                [locus],
                results,
                rule=rule(rep_samples=2, depth_reads=30, freq=0.01),
            )
            self.assertEqual(with_floor["novel_rows_by_evidence_route"], {})
            self.assertEqual(with_floor["excluded_novel"]["failed_evidence"]["rows"], 2)

    def test_read_madc_prefers_unique_numeric_preamble_row(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "vendor.csv"
            with path.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(["*", "*", "*", "*", "4403907", "4403908"])
                writer.writerow(["*", "*", "*", "*", "1", "1"])  # plate row, closest to header
                writer.writerow(
                    ["AlleleID", "CloneID", "AlleleSequence", "ClusterConsensusSequence", "X", "Y"]
                )
                writer.writerow(["m1|Ref", "m1", "ACGT", "ACGT", "3", "4"])
            parsed = mp.read_madc(path)
            self.assertEqual(parsed["samples"], ["4403907", "4403908"])
            self.assertEqual(parsed["records"][0]["counts"], [3, 4])

    def test_enumerate_variant_windows(self):
        target = "AAAACCCC" + "[A/G]" + "TTTTGGGG"
        windows = mp.enumerate_variant_windows(target, 6)
        # Six forward and six reverse windows fit around the SNP.
        self.assertEqual(len(windows), 12)
        ref, alt, offset = windows[("F", 4)]
        self.assertEqual((ref, alt, offset), ("CCCCAT", "CCCCGT", 4))
        ref, alt, offset = windows[("R", 4)]
        self.assertEqual(ref, mp.reverse_complement("CATTTT"))
        self.assertEqual(offset, 4)
        # Longer tags than the flanks allow yield nothing, as do identical alleles.
        self.assertEqual(mp.enumerate_variant_windows(target, 60), {})
        self.assertEqual(mp.enumerate_variant_windows("AAAA[C/C]TTTT", 4), {})

    def test_enumerate_variant_windows_indels(self):
        # Anchored insertion: [A/ACT] in AAAG[A/ACT]TCCC.
        windows = mp.enumerate_variant_windows("AAAG[A/ACT]TCCC", 6)
        ref, alt, offset = windows[("F", 3)]
        self.assertEqual((ref, alt, offset), ("GATCCC", "GACTTC", 2))
        # The shifted-suffix pair reproduces the vendor's fixed-length model.
        ref, alt, offset = windows[("F", 0)]
        self.assertEqual((ref, alt, offset), ("AAAGAT", "AAAGAC", 5))
        # Dash notation for a deletion behaves like the anchored form.
        dash = mp.enumerate_variant_windows("AAAGA[CT/-]CCCT", 6)
        ref, alt, offset = dash[("F", 3)]
        self.assertEqual((ref, alt, offset), ("GACTCC", "GACCCT", 3))
        # Homopolymer insertion: divergence shifts past the repeated bases.
        homo = mp.enumerate_variant_windows("AAG[T/TT]TTCA", 5)
        for _key, (r, a, off) in homo.items():
            self.assertNotEqual(r, a)
            self.assertEqual(r[off] == a[off], False)

    def test_build_panel_indel_marker(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            left = "AACCGGTT"
            # Anchored bracket [A/ACTT]: 4-bp insertion after the anchor base.
            target = "GGGG" + left + "A" * 6 + "[A/ACTT]" + "TGCATGCATGCA"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RightPrimerSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                # DArT's dash notation: insertion content in VariantAlt.
                writer.writerow(["m1", left, "", "-", "CTT", 16])
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            loci = mp.build_panel(primers, designs, 20)
            l = loci[0]
            ref_genomic = "GGGG" + left + "A" * 6 + "A" + "TGCATGCATGCA"
            alt_genomic = "GGGG" + left + "A" * 6 + "ACTT" + "TGCATGCATGCA"
            self.assertEqual(l.ref_tag, ref_genomic[4:24])
            self.assertEqual(l.alt_tag, alt_genomic[4:24])
            self.assertTrue(l.ref_tag.startswith(left) and l.alt_tag.startswith(left))
            self.assertEqual(l.ref_tag[: l.variant_pos], l.alt_tag[: l.variant_pos])
            self.assertNotEqual(l.ref_tag[l.variant_pos], l.alt_tag[l.variant_pos])
            self.assertGreaterEqual(l.variant_pos + 1, 16)
            # Round-trips through the reference-tags contract.
            out = root / "tags.csv"
            mp.write_reference_tags(out, loci)
            reloaded = mp.load_reference_tags(out, tag_length=20)
            self.assertEqual(reloaded[0].ref_tag, l.ref_tag)
            self.assertEqual(reloaded[0].alt_tag, l.alt_tag)

    def test_build_panel_degenerate_primer_overlay(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            # Oligo carries R over a flanking SNP; the genome shows A there.
            left_oligo = "AACRGG"
            target = "TTTTAACAGGA" + "C" * 40 + "[A/G]" + "C" * 40 + "TTTT"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["MarkerName", "LeftPrimerSequence", "RightPrimerSequence",
                     "VariantRef", "VariantAlt", "VariantPositionInTagSequence"]
                )
                writer.writerow(["m1", left_oligo, "", "A", "G", 48])
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            loci = mp.build_panel(primers, designs, 60)
            # The tag carries the synthesized oligo verbatim, R included,
            # matching the vendor's canonical representation.
            self.assertTrue(loci[0].ref_tag.startswith("AACRGG"))
            self.assertEqual(loci[0].ref_tag[6:11], "A" + "C" * 4)

    def test_build_panel_duplicated_flank_keeps_the_window_that_carries_the_variant(self):
        # The primer occurs twice in the target (a duplicated flank).  Only the
        # second copy carries the bracketed variant; the first copy shows the
        # Ref base at the same offset but differs one base later, so both
        # placements used to survive the Ref-base filter and the marker was
        # reported as ambiguous.  A tag must carry the designed variant.
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            primer = "ACGTACGTAC"
            copy_without = primer + "GGGGG" + "C" + "TTTTTTTATTTTTT"
            copy_with = primer + "GGGGG" + "[C/A]" + "TTTTTTTTTTTTTT"
            target = "AAAA" + copy_without + "AAAAAAAA" + copy_with + "AAAA"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    ["MarkerName", "LeftPrimerSequence", "RightPrimerSequence",
                     "VariantRef", "VariantAlt", "VariantPositionInTagSequence"]
                )
                writer.writerow(["m1", primer, "", "C", "A", 16])
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            loci = mp.build_panel(primers, designs, 30)
            self.assertEqual(loci[0].ref_tag, primer + "GGGGG" + "C" + "TTTTTTTTTTTTTT")
            self.assertEqual(loci[0].alt_tag, primer + "GGGGG" + "A" + "TTTTTTTTTTTTTT")

    def test_infer_design_resolves_indel_locus(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTACGTACGTACGTACGTACGTACGTAC"
            right = "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"
            target = left + "[A/ACGGT]" + right
            true_ref = (left + "A" + right)[23:43]
            true_alt = (left + "ACGGT" + right)[23:43]
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index, sequence in enumerate([true_ref] * 3 + [true_alt] * 3):
                    handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            output = root / "inferred.csv"
            exit_code = mp.main(
                [
                    "infer-design",
                    "--probe-design", str(design),
                    "--fastq-dir", str(fastq_dir),
                    "--tag-length", "20",
                    "--jobs", "1",
                    "--fuzzy-min-reads", "0",
                    "--output", str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            loci = mp.load_reference_tags(output, tag_length=20)
            self.assertEqual(len(loci), 1)
            self.assertEqual(loci[0].ref_tag, true_ref)
            self.assertEqual(loci[0].alt_tag, true_alt)

    def test_choose_anchor_lengths(self):
        resolved = {
            "a": ("AAAACCCCGGGGTTTT", 10),
            "b": ("AAAACCCCGGGTTTTT", 12),  # collides with a up to length 11
            "c": ("TTTTAAAACCCCGGGG", 9),
        }
        lengths = mp.choose_anchor_lengths(resolved, minimum=4, maximum=30)
        self.assertEqual(lengths["c"], 4)
        # "a" stops at its variant-position cap; "b" grows one base past it,
        # leaving prefix-nested anchors, which the trie resolves by identity.
        self.assertEqual(lengths["a"], 10)
        self.assertEqual(lengths["b"], 11)
        # Caps stop growth even when a collision cannot be resolved.
        twins = {"x": ("AAAAAAAAAA", 3), "y": ("AAAAAAAAAT", 3)}
        lengths = mp.choose_anchor_lengths(twins, minimum=2, maximum=30)
        self.assertEqual(lengths, {"x": 3, "y": 3})

    def _write_design(self, path, markers):
        with path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["MarkerName", "TargetSequence"])
            for marker, target in markers:
                writer.writerow([marker, target])

    def test_infer_design_end_to_end_exact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # SNP with 30-bp flanks; the true tag starts 7 bases before it.
            left = "ACGTACGTACGTACGTACGTACGTACGTAC"
            right = "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"
            target = left + "[A/G]" + right
            true_ref = (left + "A" + right)[23:43]
            true_alt = (left + "G" + right)[23:43]
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index, sequence in enumerate([true_ref] * 3 + [true_alt] * 2):
                    handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            output = root / "inferred.csv"
            exit_code = mp.main(
                [
                    "infer-design",
                    "--probe-design", str(design),
                    "--fastq-dir", str(fastq_dir),
                    "--tag-length", "20",
                    "--jobs", "1",
                    "--fuzzy-min-reads", "0",
                    "--output", str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            loci = mp.load_reference_tags(output, tag_length=20)
            self.assertEqual(len(loci), 1)
            self.assertEqual(loci[0].ref_tag, true_ref)
            self.assertEqual(loci[0].alt_tag, true_alt)
            self.assertEqual(loci[0].variant_pos, 7)
            self.assertTrue(true_ref.startswith(loci[0].left_primer))
            with (root / "inferred.csv.qc.json").open() as handle:
                report = json.load(handle)
            self.assertEqual(report["summary"]["resolved_exact"], 1)
            # The inferred panel must drive a full generate run.
            madc = root / "out.csv"
            exit_code = mp.main(
                [
                    "generate",
                    "--fastq-dir", str(fastq_dir),
                    "--reference-tags", str(output),
                    "--tag-length", "20",
                    "--jobs", "1",
                    "--reproducible-min-samples", "0",
                    "--quality-min-clean-reads", "1",
                    "--output", str(madc),
                ]
            )
            self.assertEqual(exit_code, 0)
            counts = {
                mp.allele_type(row["allele_id"]): row["counts"][0]
                for row in mp.read_madc(madc)["records"]
            }
            self.assertEqual((counts["Ref"], counts["Alt"]), (3, 2))

    def test_infer_design_fuzzy_recovers_divergent_locus(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            right = "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCA"
            target = left + "[C/T]" + right
            true_ref = (left + "C" + right)[30:60]  # 30-bp tag, variant offset 10
            # Every read carries the Ref target base plus two substitutions
            # placed to corrupt the later seeds but never the tag start.
            divergent = list(true_ref)
            divergent[13] = "A" if divergent[13] != "A" else "C"
            divergent[20] = "A" if divergent[20] != "A" else "C"
            divergent = "".join(divergent)
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            for sample in ("s1", "s2"):
                with gzip.open(fastq_dir / f"{sample}.fastq.gz", "wt") as handle:
                    for index in range(4):
                        handle.write(f"@r{index}\n{divergent}\n+\n{'I' * len(divergent)}\n")
            output = root / "inferred.csv"
            exit_code = mp.main(
                [
                    "infer-design",
                    "--probe-design", str(design),
                    "--fastq-dir", str(fastq_dir),
                    "--tag-length", "30",
                    "--jobs", "1",
                    "--fuzzy-min-reads", "3",
                    "--fuzzy-min-samples", "2",
                    "--output", str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            loci = mp.load_reference_tags(output, tag_length=30)
            self.assertEqual(len(loci), 1)
            self.assertEqual(loci[0].ref_tag, true_ref)
            self.assertEqual(loci[0].variant_pos, 10)
            with (root / "inferred.csv.qc.json").open() as handle:
                report = json.load(handle)
            self.assertEqual(report["per_locus"][0]["route"], "fuzzy")
            self.assertEqual(report["per_locus"][0]["supporting_samples"], 2)

    def test_infer_design_ambiguous_windows_excluded(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTACGTACGTACGTACGTACGTACGTAC"
            right = "GTCAGTCAGTCAGTCAGTCAGTCAGTCAGT"
            target = left + "[A/G]" + right
            window_a = (left + "A" + right)[23:43]
            window_b = (left + "A" + right)[15:35]
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index, sequence in enumerate([window_a] * 10 + [window_b] * 5):
                    handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            output = root / "inferred.csv"
            exit_code = mp.main(
                [
                    "infer-design",
                    "--probe-design", str(design),
                    "--fastq-dir", str(fastq_dir),
                    "--tag-length", "20",
                    "--jobs", "1",
                    "--fuzzy-min-reads", "0",
                    "--output", str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            with (root / "inferred.csv.qc.json").open() as handle:
                report = json.load(handle)
            self.assertEqual(report["summary"]["resolved"], 0)
            self.assertEqual(report["summary"]["excluded_by_reason"], {"ambiguous_windows": 1})

    def test_choose_tag_windows_invariants(self):
        ref_g = "ACGT" * 10 + "A" + "TGCA" * 10  # SNP at index 40
        expansions = {"m1": (ref_g, ref_g[:40] + "G" + ref_g[41:], 40)}
        # center placement honors the anchor-safety invariant
        w = mp.choose_tag_windows(expansions, {}, 21, 8, "center")["m1"]
        self.assertLessEqual(w["start"], 40 - 8)
        self.assertGreaterEqual(40 - w["start"], 8)  # variant outside anchor
        self.assertLess(40 - w["start"], 21)         # variant inside window
        # auto placement follows read coverage
        cover = {("m1", 30, 40): 100}  # reads span target coords [30, 70)
        w = mp.choose_tag_windows(expansions, cover, 21, 8, "auto")["m1"]
        self.assertGreaterEqual(w["start"], 30)
        self.assertLessEqual(w["start"] + 21, 70)
        self.assertEqual(w["coverage"], 1.0)

    def test_choose_tag_windows_str_penalty(self):
        self.assertEqual(mp.max_tandem_run("ACACACACACAC"), 12)   # (AC)6
        self.assertEqual(mp.max_tandem_run("AAAAA"), 5)           # homopolymer
        self.assertEqual(mp.max_tandem_run("ACGTACGTACGT"), 1)    # period 4: not counted
        # SNP at 40 with a dinucleotide repeat immediately downstream: a window
        # holding the full repeat has slightly better coverage than one that
        # truncates it; the soft penalty must prefer the truncating window.
        ref_g = "GATTCTGACCGGTTAACCGGATCCGGAATTCGGCCATTGA" + "A" + "ACACACACACAC" + "GGATTCCGGATCAGGATTAACCGGTT"
        expansions = {"m1": (ref_g, ref_g[:40] + "G" + ref_g[41:], 40)}
        # near-tied coverage: the repeat-spanning start (30+) barely beats the
        # repeat-truncating start (28); only the soft penalty can flip it
        spanning = {("m1", 30, 41): 1, ("m1", 22, 41): 1000}
        with_pen = mp.choose_tag_windows(expansions, spanning, 21, 8, "auto")["m1"]
        without = mp.choose_tag_windows(expansions, spanning, 21, 8, "auto", str_penalty=0.0)["m1"]
        self.assertGreater(without["max_run"], with_pen["max_run"])
        self.assertLess(with_pen["start"], without["start"])

    def test_choose_tag_windows_group_dedup(self):
        # two tiles of one locus share the frame and the coverage profile;
        # the profile pushes both onto one start -> identical Ref tags ->
        # every ref read would tie between the markers and be dropped.
        # The dedup must hand the second member its best DISTINCT start.
        ref_g = ("GATTCTGACCGGTTAACCGGATCCGGAATT" + "A"
                 + "CGT" + "G" + "GCCATTGACCGGATTCCGGATCAGGATT")
        alt1 = ref_g[:30] + "C" + ref_g[31:]      # snp at 30
        alt2 = ref_g[:34] + "T" + ref_g[35:]      # snp at 34
        expansions = {"tA": (ref_g, alt1, 30), "tB": (ref_g, alt2, 34)}
        groups = {"tA": "locus1", "tB": "locus1"}
        cover = {("tA", 10, 31): 100, ("tB", 10, 31): 100}  # full cover only for starts 10..20
        w = mp.choose_tag_windows(expansions, cover, 21, 8, "auto", groups=groups)
        self.assertNotEqual(w["tA"]["start"], w["tB"]["start"])
        self.assertIn("dedup_nudged_from", w["tB"])
        self.assertEqual(w["tB"]["dedup_nudged_from"], w["tA"]["start"])
        # without groups the same inputs are allowed to collide (pre-0.9 behavior)
        w2 = mp.choose_tag_windows(expansions, cover, 21, 8, "auto")
        self.assertEqual(w2["tA"]["start"], w2["tB"]["start"])
        # a homopolymer frame has ONE window content everywhere: no distinct
        # start exists and the guard must refuse rather than emit twins
        poly = "A" * 60
        exp_poly = {"pA": (poly, poly[:30] + "C" + poly[31:], 30),
                    "pB": (poly, poly[:34] + "C" + poly[35:], 34)}
        with self.assertRaises(mp.MadcError):
            mp.choose_tag_windows(exp_poly, {}, 21, 8, "auto",
                                  groups={"pA": "locus1", "pB": "locus1"})

    def test_normalize_reads_end_to_end(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTTGCAACGGATCC"
            right = "GGATCCAAGGTTCCAA"
            filler_a = "TTGACCTGAGGACTGTAGCA"
            filler_b = "CATGGCTAAGCTTGGCATCA"
            target = left + filler_a + "[A/G]" + filler_b + right
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            ref_g = left + filler_a + "A" + filler_b + right
            alt_g = left + filler_a + "G" + filler_b + right
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            reads = [
                ref_g[2:66],                        # forward, variable start
                mp.reverse_complement(alt_g[5:69]), # reverse orientation
                ref_g[0:30],                        # locates but cannot cover any SNP window
                "T" * 50,                           # no panel k-mer
            ]
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index, sequence in enumerate(reads):
                    handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            out_dir = root / "normalized"
            exit_code = mp.main(
                [
                    "normalize-reads",
                    "--probe-design", str(design),
                    "--fastq-dir", str(fastq_dir),
                    "--tag-length", "21",
                    "--anchor-length", "8",
                    "--jobs", "1",
                    "--out-dir", str(out_dir),
                ]
            )
            self.assertEqual(exit_code, 0)
            report = json.loads((out_dir / "normalize_reads.qc.json").read_text())
            self.assertEqual(report["totals"]["reads_extracted"], 2)
            self.assertEqual(report["totals"]["reads_no_panel_kmer"], 1)
            self.assertEqual(report["totals"]["reads_window_not_covered"], 1)
            loci = mp.load_reference_tags(out_dir / "window_reference_tags.csv", tag_length=21)
            self.assertEqual(len(loci), 1)
            self.assertGreaterEqual(loci[0].variant_pos, 8)
            madc = root / "out.csv"
            exit_code = mp.main(
                [
                    "generate",
                    "--fastq-dir", str(out_dir),
                    "--fastq-glob", "*.fastq",
                    "--reference-tags", str(out_dir / "window_reference_tags.csv"),
                    "--tag-length", "21",
                    "--jobs", "1",
                    "--reproducible-min-samples", "0",
                    "--depth-min-reads", "1",
                    "--output", str(madc),
                ]
            )
            self.assertEqual(exit_code, 0)
            counts = {
                mp.allele_type(row["allele_id"]): row["counts"][0]
                for row in mp.read_madc(madc)["records"]
            }
            self.assertEqual((counts["Ref"], counts["Alt"]), (1, 1))

    def test_normalize_reads_reuse_tags(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTTGCAACGGATCC"
            right = "GGATCCAAGGTTCCAA"
            filler_a = "TTGACCTGAGGACTGTAGCA"
            filler_b = "CATGGCTAAGCTTGGCATCA"
            target = left + filler_a + "[A/G]" + filler_b + right
            design = root / "design.tsv"
            self._write_design(design, [("m1", target)])
            ref_g = left + filler_a + "A" + filler_b + right
            alt_g = left + filler_a + "G" + filler_b + right
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index, sequence in enumerate([ref_g[2:66], mp.reverse_complement(alt_g[5:69])]):
                    handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            base = ["--probe-design", str(design), "--fastq-dir", str(fastq_dir),
                    "--tag-length", "21", "--jobs", "1"]
            out_fresh = root / "fresh"
            self.assertEqual(
                mp.main(["normalize-reads", *base, "--anchor-length", "8", "--out-dir", str(out_fresh)]), 0
            )
            panel = out_fresh / "window_reference_tags.csv"

            # identity reuse: byte-identical panel and extracted FASTQ; the coverage
            # pass now scores the frozen windows (freeze health) without touching them
            out_reuse = root / "reuse"
            self.assertEqual(
                mp.main(["normalize-reads", *base, "--reuse-tags", str(panel), "--out-dir", str(out_reuse)]), 0
            )
            self.assertEqual(panel.read_bytes(), (out_reuse / "window_reference_tags.csv").read_bytes())
            self.assertEqual((out_fresh / "s1.fastq").read_bytes(), (out_reuse / "s1.fastq").read_bytes())
            report = json.loads((out_reuse / "normalize_reads.qc.json").read_text())
            self.assertEqual(report["parameters"]["window_placement"], "reuse")
            health = report["freeze_health"]
            self.assertEqual(health["markers_scored"], 1)
            self.assertEqual(health["markers_moved_by_fresh_placement"], 0)
            self.assertEqual(health["markers_losing_over_0.10"], 0)
            self.assertAlmostEqual(health["read_weighted_full_coverage_frozen"], health["read_weighted_full_coverage_fresh"])
            marker = report["markers"]["m1"]
            self.assertEqual(marker["fresh_window_start"], marker["window_start"])
            self.assertIsInstance(marker["predicted_full_coverage"], float)
            # center placement with reuse skips the scoring pass
            out_center = root / "reuse_center"
            self.assertEqual(
                mp.main(["normalize-reads", *base, "--reuse-tags", str(panel), "--window-placement", "center",
                         "--out-dir", str(out_center)]), 0
            )
            center_report = json.loads((out_center / "normalize_reads.qc.json").read_text())
            self.assertNotIn("freeze_health", center_report)
            self.assertIsNone(center_report["markers"]["m1"]["predicted_full_coverage"])
            self.assertEqual((out_fresh / "s1.fastq").read_bytes(), (out_center / "s1.fastq").read_bytes())

            # shifted reuse: a hand-shifted window is honored verbatim
            loci = mp.load_reference_tags(panel, tag_length=21)
            start = ref_g.find(loci[0].ref_tag) - 2
            shifted = mp.Locus(
                index=0, marker="m1",
                left_primer=ref_g[start : start + 8], right_primer="",
                ref_base="A", alt_base="G",
                variant_pos=loci[0].variant_pos + 2,
                ref_tag=ref_g[start : start + 21], alt_tag=alt_g[start : start + 21],
            )
            shifted_panel = root / "shifted.csv"
            mp.write_reference_tags(shifted_panel, [shifted])
            out_shift = root / "reuse_shifted"
            self.assertEqual(
                mp.main(["normalize-reads", *base, "--reuse-tags", str(shifted_panel), "--out-dir", str(out_shift)]), 0
            )
            extracted = (out_shift / "s1.fastq").read_text().splitlines()[1]
            self.assertEqual(extracted, shifted.ref_tag)
            shift_report = json.loads((out_shift / "normalize_reads.qc.json").read_text())
            self.assertEqual(shift_report["freeze_health"]["markers_moved_by_fresh_placement"], 1)
            self.assertNotEqual(shift_report["markers"]["m1"]["fresh_window_start"], shift_report["markers"]["m1"]["window_start"])

            # a panel that disagrees with the design target is refused
            bad = mp.Locus(
                index=0, marker="m1", left_primer=shifted.left_primer, right_primer="",
                ref_base="A", alt_base="G", variant_pos=shifted.variant_pos,
                ref_tag=shifted.ref_tag[:-1] + ("C" if shifted.ref_tag[-1] != "C" else "T"),
                alt_tag=shifted.alt_tag[:-1] + ("C" if shifted.alt_tag[-1] != "C" else "T"),
            )
            bad_panel = root / "bad.csv"
            mp.write_reference_tags(bad_panel, [bad])
            with self.assertRaises(SystemExit) as ctx:
                mp.main(["normalize-reads", *base, "--reuse-tags", str(bad_panel), "--out-dir", str(root / "x")])
            self.assertEqual(ctx.exception.code, 2)

    def test_normalize_reads_variable_anchors(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            common = "ACGTTGCAAC"
            # m1/m2 windows share their first 8 bases then diverge before the
            # variant -> anchors must grow beyond the 8-bp minimum, uniquely.
            # m3/m4 windows are identical all the way to the variant -> they
            # legitimately share a capped anchor (resolved by tag identity).
            t1 = "T" * 10 + common + "GT" + "[A/G]" + "CATGGCTAAGCTTGG"
            t2 = "T" * 10 + common + "CT" + "[A/G]" + "GGCCAATTGGCCAAT"
            t3 = "C" * 10 + common + "AA" + "[A/G]" + "CCGGTTAACCGGTTA"
            t4 = "C" * 10 + common + "AA" + "[A/G]" + "TTAACCGGTTAACCG"
            design = root / "design.tsv"
            self._write_design(design, [("m1", t1), ("m2", t2), ("m3", t3), ("m4", t4)])
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                handle.write("@r0\n" + "T" * 30 + "\n+\n" + "I" * 30 + "\n")
            out_dir = root / "normalized"
            exit_code = mp.main(
                ["normalize-reads", "--probe-design", str(design),
                 "--fastq-dir", str(fastq_dir), "--tag-length", "21",
                 "--anchor-length", "8", "--window-placement", "center",
                 "--jobs", "1", "--out-dir", str(out_dir)]
            )
            self.assertEqual(exit_code, 0)
            loci = {l.marker: l for l in mp.load_reference_tags(out_dir / "window_reference_tags.csv")}
            self.assertEqual(len(loci), 4)
            self.assertGreater(len(loci["m1"].left_primer), 8)      # grew past minimum
            self.assertNotEqual(loci["m1"].left_primer, loci["m2"].left_primer)
            self.assertEqual(loci["m3"].left_primer, loci["m4"].left_primer)  # capped share
            report = json.loads((out_dir / "normalize_reads.qc.json").read_text())
            self.assertEqual(report["anchor_summary"]["shared_anchor_groups"], [["m3", "m4"]])

    def test_locus_group_tile_routing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # one 120bp locus, SNPs at 30 (A/G) and 90 (C/T): two tiles
            left = "ACGTTGCAACGGATCCTTGACCTGAGGACT"          # 0-29
            mid = "CATGGCTAAGCTTGGCATCAGGCCAATTGGCCAATACGTACGTTGGAACCTTGGAACCT"  # 31-90
            right = "TTGACGTAGCATCAGCCATGACTTAGCAT"
            frame = left + "A" + mid[:-1] + "C" + right
            t1 = left + "[A/G]" + mid[:-1] + "C" + right       # bracket at 30
            t2 = left + "A" + mid[:-1] + "[C/T]" + right       # bracket at 90
            design = root / "design.tsv"
            with design.open("w") as fh:
                fh.write("MarkerName\tTargetSequence\tLocusGroup\n")
                fh.write(f"tileA\t{t1}\tlocus1\n")
                fh.write(f"tileB\t{t2}\tlocus1\n")
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            # left-end fragments cover only tileA's window; right-end only tileB's
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                left_read = frame[:55]
                right_read = mp.reverse_complement(frame[-55:])
                for index in range(4):
                    handle.write(f"@l{index}\n{left_read}\n+\n{'I'*len(left_read)}\n")
                for index in range(6):
                    handle.write(f"@r{index}\n{right_read}\n+\n{'I'*len(right_read)}\n")
            out_dir = root / "normalized"
            exit_code = mp.main(
                ["normalize-reads", "--probe-design", str(design),
                 "--fastq-dir", str(fastq_dir), "--tag-length", "31",
                 "--anchor-length", "8", "--window-placement", "center",
                 "--jobs", "1", "--out-dir", str(out_dir)]
            )
            self.assertEqual(exit_code, 0)
            report = json.loads((out_dir / "normalize_reads.qc.json").read_text())
            self.assertEqual(report["totals"]["reads_extracted"], 10)   # every read, once
            self.assertEqual(
                report["locus_groups"],
                {"groups": 1, "tiled_groups": 1, "window_dedup_nudged": []},
            )
            loci = {l.marker: l for l in mp.load_reference_tags(out_dir / "window_reference_tags.csv")}
            self.assertEqual({loci["tileA"].locus_group, loci["tileB"].locus_group}, {"locus1"})
            madc = root / "out.csv"
            exit_code = mp.main(
                ["generate", "--fastq-dir", str(out_dir), "--fastq-glob", "*.fastq",
                 "--reference-tags", str(out_dir / "window_reference_tags.csv"),
                 "--tag-length", "31", "--jobs", "1", "--output", str(madc)]
            )
            self.assertEqual(exit_code, 0)
            counts = {}
            for row in mp.read_madc(madc)["records"]:
                counts[(row["clone_id"], mp.allele_type(row["allele_id"]))] = row["counts"][0]
            self.assertEqual(counts[("tileA", "Ref")], 4)   # left fragments -> tileA only
            self.assertEqual(counts[("tileB", "Ref")], 6)   # right fragments -> tileB only

    def test_locus_group_covered_routing(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # same two-tile locus as above, but FULL-SPAN fragments that
            # cover both windows (the merged-pairs geometry)
            left = "ACGTTGCAACGGATCCTTGACCTGAGGACT"
            mid = "CATGGCTAAGCTTGGCATCAGGCCAATTGGCCAATACGTACGTTGGAACCTTGGAACCT"
            right = "TTGACGTAGCATCAGCCATGACTTAGCAT"
            frame = left + "A" + mid[:-1] + "C" + right
            t1 = left + "[A/G]" + mid[:-1] + "C" + right
            t2 = left + "A" + mid[:-1] + "[C/T]" + right
            design = root / "design.tsv"
            with design.open("w") as fh:
                fh.write("MarkerName\tTargetSequence\tLocusGroup\n")
                fh.write(f"tileA\t{t1}\tlocus1\n")
                fh.write(f"tileB\t{t2}\tlocus1\n")
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            with gzip.open(fastq_dir / "s1.fastq.gz", "wt") as handle:
                for index in range(10):
                    handle.write(f"@f{index}\n{frame}\n+\n{'I'*len(frame)}\n")

            def run(routing, out_name):
                out_dir = root / out_name
                exit_code = mp.main(
                    ["normalize-reads", "--probe-design", str(design),
                     "--fastq-dir", str(fastq_dir), "--tag-length", "31",
                     "--anchor-length", "8", "--window-placement", "center",
                     "--group-routing", routing,
                     "--jobs", "1", "--out-dir", str(out_dir)]
                )
                self.assertEqual(exit_code, 0)
                madc = root / f"{out_name}.csv"
                exit_code = mp.main(
                    ["generate", "--fastq-dir", str(out_dir), "--fastq-glob", "*.fastq",
                     "--reference-tags", str(out_dir / "window_reference_tags.csv"),
                     "--tag-length", "31", "--jobs", "1", "--output", str(madc)]
                )
                self.assertEqual(exit_code, 0)
                counts = {}
                for row in mp.read_madc(madc)["records"]:
                    counts[(row["clone_id"], mp.allele_type(row["allele_id"]))] = row["counts"][0]
                report = json.loads((out_dir / "normalize_reads.qc.json").read_text())
                return counts, report["totals"]

            # nearest: full-span fragments all share one center -> winner-take-
            # all at a single tile (the measured production pathology)
            counts, totals = run("nearest", "norm_nearest")
            self.assertEqual(totals["reads_extracted"], 10)
            self.assertNotIn("tile_extra_emissions", totals)
            self.assertEqual(counts.get(("tileA", "Ref"), 0)
                             + counts.get(("tileB", "Ref"), 0), 10)
            self.assertIn(0, (counts.get(("tileA", "Ref"), 0),
                              counts.get(("tileB", "Ref"), 0)))
            # covered: each fragment genotypes BOTH tiles, once per window
            counts, totals = run("covered", "norm_covered")
            self.assertEqual(totals["reads_extracted"], 10)
            self.assertEqual(totals["tile_extra_emissions"], 10)
            self.assertEqual(counts[("tileA", "Ref")], 10)
            self.assertEqual(counts[("tileB", "Ref")], 10)

    def test_window_scan_prefers_informative_length(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            # target SNP at 40, a linked SNP at 60: a 21bp window sees only
            # the target, a 41bp window captures the linked haplotype too ->
            # equal retention, higher information at 41
            left = "ACGTTGCAACGGATCCTTGACCTGAGGACTGTAGCACATG"      # 0-39
            mid = "GCTAAGCTTGGCATCAGGC"                             # 41-59
            right = "TATTGACCGGATTCCGGATCAGGATTAACC"                # 61-90
            frame = left + "A" + mid + "C" + right
            t1 = left + "[A/G]" + mid + "C" + right
            design = root / "design.tsv"
            self._write_design(design, [("m1", t1)])
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            linked = frame[:60] + "T" + frame[61:]
            for sample in ("s1", "s2"):
                with gzip.open(fastq_dir / f"{sample}.fastq.gz", "wt") as handle:
                    for index in range(15):
                        handle.write(f"@a{index}\n{frame}\n+\n{'I'*len(frame)}\n")
                    for index in range(15):
                        handle.write(f"@b{index}\n{linked}\n+\n{'I'*len(linked)}\n")
            out_dir = root / "scan"
            exit_code = mp.main(
                ["window-scan", "--probe-design", str(design),
                 "--fastq-dir", str(fastq_dir), "--lengths", "21,41",
                 "--scan-samples", "2", "--anchor-length", "8",
                 "--jobs", "1", "--out-dir", str(out_dir)]
            )
            self.assertEqual(exit_code, 0)
            report = json.loads((out_dir / "window_scan.qc.json").read_text())
            by_len = {r["tag_length"]: r for r in report["candidates"]}
            self.assertEqual(set(by_len), {21, 41})
            self.assertGreater(by_len[41]["novel_rows_per_100_markers"],
                               by_len[21]["novel_rows_per_100_markers"])
            self.assertEqual(report["recommended_tag_length"], 41)
            self.assertTrue((out_dir / "window_scan.tsv").exists())
            self.assertTrue((out_dir / "len_41" / "scan_MADC.csv").exists())

    def test_window_scan_draws_report_agreement(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            left = "ACGTTGCAACGGATCCTTGACCTGAGGACTGTAGCACATG"
            mid = "GCTAAGCTTGGCATCAGGC"
            right = "TATTGACCGGATTCCGGATCAGGATTAACC"
            frame = left + "A" + mid + "C" + right
            design = root / "design.tsv"
            self._write_design(design, [("m1", left + "[A/G]" + mid + "C" + right)])
            fastq_dir = root / "raw"
            fastq_dir.mkdir()
            linked = frame[:60] + "T" + frame[61:]
            for sample in ("s1", "s2", "s3"):
                with gzip.open(fastq_dir / f"{sample}.fastq.gz", "wt") as handle:
                    for index in range(15):
                        handle.write(f"@a{index}\n{frame}\n+\n{'I'*len(frame)}\n")
                    for index in range(15):
                        handle.write(f"@b{index}\n{linked}\n+\n{'I'*len(linked)}\n")
            out_dir = root / "scan"
            exit_code = mp.main(
                ["window-scan", "--probe-design", str(design), "--fastq-dir", str(fastq_dir),
                 "--lengths", "21,41", "--scan-samples", "2", "--scan-draws", "2", "--scan-seed", "7",
                 "--anchor-length", "8", "--jobs", "1", "--out-dir", str(out_dir)]
            )
            self.assertEqual(exit_code, 0)
            report = json.loads((out_dir / "window_scan.qc.json").read_text())
            self.assertEqual(len(report["draws"]), 2)
            self.assertEqual([d["seed"] for d in report["draws"]], [7, 8])
            self.assertEqual(report["recommended_tag_length"], 41)
            self.assertEqual(report["recommendation_votes"], {"41": 2})
            self.assertEqual(report["recommendation_stability"], 1.0)
            self.assertEqual(report["by_length"]["41"]["draws_recommended"], 2)
            for draw in (1, 2):
                self.assertTrue((out_dir / f"draw_{draw}" / "window_scan.qc.json").exists())
                self.assertTrue((out_dir / f"draw_{draw}" / "len_41" / "scan_MADC.csv").exists())
            self.assertTrue((out_dir / "window_scan.tsv").exists())

    def test_merge_read_pair_prefers_longest_overlap_in_repeats(self):
        # 120bp fragment whose true 40bp overlap is pure (AC)n: every
        # in-phase shorter overlap also matches perfectly, so shortest-first
        # acceptance would insert phantom repeat units (120 -> 144bp).
        left = "GATTACAGGCTTACGATCGTAGCCTAGGTTAACGGATCCAA"[:40]
        right = "TTGGATCCGTTAACCTAGGCTACGATCGTAAGCCTGTAATC"[:40]
        fragment = left + "AC" * 20 + right
        r1, r2 = fragment[:80], mp.reverse_complement(fragment[40:])
        merged = mp.merge_read_pair(r1, "I" * 80, r2, "I" * 80, 15, 0.1)
        self.assertEqual(merged[0], fragment)

    def test_primer_index_iupac_degenerate_primers(self):
        ref = "AARACCGGTTGGCCAA"          # R = A/G at position 2
        alt = ref[:8] + "A" + ref[9:]
        locus = mp.Locus(0, "deg", "AARACC", "", "G", "A", 8, ref, alt)
        index = mp.PrimerIndex([locus])
        # either expansion matches with a zero budget: the degenerate base
        # is the oligo pool's own variability, not a mismatch
        self.assertEqual(index.match("AAAACCGGTTGGCCAA", 0), [0])
        self.assertEqual(index.match("AAGACCGGTTGGCCAA", 0), [0])
        # a real sequencing error elsewhere still fits the default budget
        self.assertEqual(index.match("AATACCGGTTGGCCAA", 1), [0])
        # a read N satisfies nothing: budget spent, second error rejected
        self.assertEqual(index.match("AANACCGGTTGGCCAA", 0), [])
        self.assertEqual(index.match("AANTCCGGTTGGCCAA", 1), [])

    def test_probe_design_rejects_multiallelic_bracket(self):
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "design.csv"
            path.write_text(
                "MarkerName,TargetSequence\ntri,AACCGG[A/G/T]TTGGCC\n"
            )
            with self.assertRaisesRegex(mp.MadcError, "more than two alleles"):
                mp.read_probe_design(path)

    def test_merge_read_pair_geometries(self):
        rng_frag = "ACGTTGCAACGGATCCTTGACCTGAGGACTGTAGCACATGGCTAAGCTTGGCATCAGGCCAATTGGCCAATACGTACGTA"
        # extended: 80bp fragment, 50bp reads, 20bp true overlap
        r1, r2 = rng_frag[:50], mp.reverse_complement(rng_frag[30:])
        merged = mp.merge_read_pair(r1, "I"*50, r2, "I"*50, 15, 0.1)
        self.assertEqual(merged[0], rng_frag)
        # engulfed / read-through: 40bp fragment inside 50bp reads with junk tails
        frag = rng_frag[:40]
        r1 = frag + "GGGGGGGGGG"
        r2 = mp.reverse_complement("TTTTTTTTTT" + frag)
        merged = mp.merge_read_pair(r1, "I"*50, r2, "I"*50, 15, 0.1)
        self.assertEqual(merged[0], frag)
        # quality-resolved disagreement inside the overlap
        r1 = rng_frag[:50]
        bad = list(rng_frag[30:])
        bad[5] = "A" if bad[5] != "A" else "C"      # error in R2 at overlap pos
        r2 = mp.reverse_complement("".join(bad))
        merged = mp.merge_read_pair(r1, "I"*50, r2, "#"*50, 15, 0.1)
        self.assertEqual(merged[0], rng_frag)       # high-qual R1 base wins
        # no overlap -> None
        self.assertIsNone(mp.merge_read_pair("A"*50, "I"*50, "C"*50, "I"*50, 15, 0.1))

    def test_dart_preamble_and_sample_lookup(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            ref_tag = "ACGTTGCA" + "TT" + "A" + "GGCCAATTGG"
            alt_tag = ref_tag[:10] + "G" + ref_tag[11:]
            locus = mp.Locus(
                index=0, marker="m1", left_primer="ACGTTGCA", right_primer="",
                ref_base="A", alt_base="G", variant_pos=10, ref_tag=ref_tag, alt_tag=alt_tag,
            )
            tags = root / "tags.csv"
            mp.write_reference_tags(tags, [locus])
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            for name, sequence in (("s1", ref_tag), ("s2", alt_tag)):
                with gzip.open(fastq_dir / f"{name}.fastq.gz", "wt") as handle:
                    for index in range(3):
                        handle.write(f"@r{index}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            lookup = root / "lookup.csv"
            lookup.write_text("fastq_id,sample_name\ns1,SampleA\ns2,SampleB\n")
            output = root / "out.csv"
            base = ["generate", "--fastq-dir", str(fastq_dir), "--reference-tags", str(tags),
                    "--tag-length", "21", "--jobs", "1", "--output", str(output)]
            exit_code = mp.main(
                base + ["--dart-preamble", "--preamble-project", "TESTPROJ",
                        "--sample-lookup", str(lookup)]
            )
            self.assertEqual(exit_code, 0)
            with output.open() as handle:
                rows = list(csv.reader(handle))
            width = len(mp.MADC_METADATA_COLUMNS)
            # provenance lives in the first cell only; every other valueless
            # preamble cell is '*', vendor-style
            self.assertTrue(rows[0][0].startswith(f"* MADCATS v{mp.__version__} "))
            self.assertIn("tag_length=21", rows[0][0])
            self.assertIn("min_supporting_freq=0.0", rows[0][0])
            self.assertEqual(rows[0][1:width], ["*"] * (width - 1))
            for row in rows[1:6]:
                self.assertEqual(row, ["*"] * len(row))      # plate rows all starred
            self.assertEqual(rows[6][:width], ["*"] * width)
            self.assertEqual(rows[0][width:], ["TESTPROJ", "TESTPROJ"])
            self.assertEqual(rows[6][width:], ["s1", "s2"])  # targetid row: FASTQ ids
            self.assertEqual(rows[7][0], "AlleleID")
            self.assertEqual(rows[7][width:], ["SampleA", "SampleB"])
            # body: uncomputed statistics columns are '.', never blank
            one_ratio_ref = mp.MADC_METADATA_COLUMNS.index("OneRatioRef")
            for row in rows[8:]:
                self.assertEqual(row[one_ratio_ref], ".")
                self.assertNotIn("", row[:width])
            parsed = mp.read_madc(output)
            self.assertEqual(parsed["samples"], ["SampleA", "SampleB"])
            counts = {
                mp.allele_type(row["allele_id"]): row["counts"] for row in parsed["records"]
            }
            self.assertEqual(counts["Ref"], [3, 0])
            self.assertEqual(counts["Alt"], [0, 3])

            # incomplete lookup: unmapped run samples must be a hard error
            (root / "partial.csv").write_text("fastq_id,sample_name\ns1,SampleA\n")
            with self.assertRaises(SystemExit) as ctx:
                mp.main(base[:-2] + ["--output", str(root / "out2.csv"),
                                     "--sample-lookup", str(root / "partial.csv")])
            self.assertEqual(ctx.exception.code, 2)

    def test_design_reconstruction_and_small_generation(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            primers = root / "primers.csv"
            designs = root / "design.tsv"
            with primers.open("w", newline="") as handle:
                writer = csv.writer(handle)
                writer.writerow(
                    [
                        "MarkerName",
                        "LeftPrimerSequence",
                        "RightPrimerSequence",
                        "VariantRef",
                        "VariantAlt",
                        "VariantPositionInTagSequence",
                    ]
                )
                writer.writerow(["m1", "AACCGG", "TTTT", "A", "G", 7])
            target = "TTTTAACCGGA" + "C" * 80 + "TTTT"
            target = target[:10] + "[A/G]" + target[11:]
            with designs.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(["MarkerName", "TargetSequence"])
                writer.writerow(["m1", target])
            loci = mp.build_panel(primers, designs, 20)
            self.assertEqual(loci[0].ref_tag[:8], "AACCGGAC")
            self.assertEqual(loci[0].alt_tag[6], "G")

            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            ref = loci[0].ref_tag
            alt = loci[0].alt_tag
            mutate = lambda text: "".join("A" if base != "A" else "C" for base in text)
            high_quality_novel = ref[:10] + mutate(ref[10:14]) + ref[14:]
            low_quality_novel = ref[:14] + mutate(ref[14:18]) + ref[18:]
            with gzip.open(fastq_dir / "sample1.FASTQ.gz", "wt") as handle:
                reads = [
                    (ref, "I" * len(ref)),
                    (ref, "I" * len(ref)),
                    (alt, "I" * len(alt)),
                    (high_quality_novel, "I" * len(high_quality_novel)),
                    (low_quality_novel, "#" * len(low_quality_novel)),
                ]
                # Repeat the read list to exercise the memoized decision path.
                for index, (sequence, quality) in enumerate(reads * 3):
                    handle.write(f"@r{index}\n{sequence}\n+\n{quality}\n")
            output = root / "out.csv"
            exit_code = mp.main(
                [
                    "generate",
                    "--fastq-dir",
                    str(fastq_dir),
                    "--primers",
                    str(primers),
                    "--probe-design",
                    str(designs),
                    "--tag-length",
                    "20",
                    "--jobs",
                    "1",
                    "--quality-min-clean-reads",
                    "1",
                    "--reproducible-min-samples",
                    "0",
                    "--output",
                    str(output),
                ]
            )
            self.assertEqual(exit_code, 0)
            parsed = mp.read_madc(output)
            counts = {mp.allele_type(row["allele_id"]): row["counts"][0] for row in parsed["records"]}
            self.assertEqual(counts["Ref"], 6)
            self.assertEqual(counts["Alt"], 3)
            sequences = {row["sequence"] for row in parsed["records"]}
            self.assertIn(high_quality_novel, sequences)
            self.assertNotIn(low_quality_novel, sequences)
            with (root / "out.csv.qc.json").open() as handle:
                qc = json.load(handle)
            accounting = qc["read_accounting"]
            self.assertEqual(accounting["reads_total"], 15)
            self.assertEqual(
                accounting["reads_in_final_count_table"],
                accounting["reads_in_emitted_rows"]
                + accounting["reads_dropped_below_min_global_count"]
                + accounting["reads_dropped_failed_evidence"],
            )


if __name__ == "__main__":
    unittest.main()


class ReviewFixTests(unittest.TestCase):
    """0.10: gapped rescue, cheaper-explanation classification, and the review fixes."""

    @staticmethod
    def _prefix_distances_dp(pattern, text):
        m = len(pattern)
        previous = list(range(m + 1))  # D[i][0] = i
        out = [m]
        for j, tb in enumerate(text, start=1):
            current = [j] + [0] * m
            for i in range(1, m + 1):
                current[i] = min(previous[i - 1] + (pattern[i - 1] != tb), previous[i] + 1, current[i - 1] + 1)
            out.append(current[m])
            previous = current
        return out

    def test_myers_prefix_distances_match_dynamic_programming(self):
        import random as _random
        rng = _random.Random(11)
        for _ in range(300):
            a = "".join(rng.choice("ACGT") for _ in range(rng.randint(1, 45)))
            b = list(a)
            for _ in range(rng.randint(0, 5)):
                op = rng.choice("sid")
                p = rng.randrange(len(b) + 1)
                if op == "s" and p < len(b):
                    b[p] = rng.choice("ACGT")
                elif op == "i":
                    b.insert(p, rng.choice("ACGT"))
                elif op == "d" and b:
                    del b[min(p, len(b) - 1)]
            b = "".join(b)
            expected = self._prefix_distances_dp(a, b)
            self.assertEqual(mp._myers_prefix_distances(a, b, design_is_pattern=False), expected, (a, b))
            self.assertEqual(mp._myers_prefix_distances(a, b, design_is_pattern=True), expected, (a, b))

    def test_gapped_prescreen_never_exceeds_banded_distance(self):
        import random as _random
        rng = _random.Random(5)
        for _ in range(200):
            ref = "".join(rng.choice("ACGT") for _ in range(41))
            read = list(ref)
            p = rng.randrange(5, 35)
            if rng.random() < 0.5:
                del read[p]
                read.append(rng.choice("ACGT"))
            else:
                read.insert(p, rng.choice("ACGT"))
                read = read[:41]
            read = "".join(read)
            banded = mp.gapped_stats(read, ref, 10)[0]
            self.assertLessEqual(mp.gapped_prescreen(read, ref, 10), banded)
            self.assertLessEqual(banded, 1)

    def _indel_locus(self):
        ref = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"  # 41 bp
        variant_pos = 20
        alt = ref[:variant_pos] + ("C" if ref[variant_pos] != "C" else "T") + ref[variant_pos + 1:]
        return mp.Locus(0, "m1", ref[:8], "", ref[variant_pos], alt[variant_pos], variant_pos, ref, alt)

    def test_gapped_rescue_assigns_and_classifies_through_the_gap(self):
        locus = self._indel_locus()
        ref = locus.ref_tag
        budget = 8  # floor(0.2 * 41)
        # 1-bp deletion at position 10: the suffix shifts, Hamming explodes.
        deleted = ref[:10] + ref[11:] + "A"
        self.assertGreater(mp.tag_mismatches(deleted, ref), budget)
        self.assertEqual(mp.select_marker(deleted, [0], [locus], budget), (None, "offtarget"))
        self.assertEqual(mp.select_marker(deleted, [0], [locus], budget, gap_band=10), (0, "assigned_gapped"))
        # The target (ref position 20) sits at read index 19 and carries the ref base.
        self.assertEqual(mp.classify_sequence(locus, deleted), "RefMatch")
        # A read carrying the alt base beyond the same deletion is AltMatch, not a shifted misread.
        alt_deleted = locus.alt_tag[:10] + locus.alt_tag[11:] + "A"
        self.assertEqual(mp.classify_sequence(locus, alt_deleted), "AltMatch")
        # A deletion that removes the target itself is missing data.
        target_deleted = ref[:20] + ref[21:] + "A"
        self.assertEqual(mp.select_marker(target_deleted, [0], [locus], budget, gap_band=10)[1], "assigned_gapped")
        self.assertEqual(mp.classify_sequence(locus, target_deleted), "Other")
        # Two events (deletion + insertion) are still one gapped reading.
        two_events = ref[:10] + ref[11:30] + "G" + ref[30:]
        self.assertEqual(len(two_events), 41)
        self.assertEqual(mp.select_marker(two_events, [0], [locus], budget, gap_band=10), (0, "assigned_gapped"))

    def test_gapped_rescue_rejects_substitution_heavy_reads(self):
        locus = self._indel_locus()
        ref = list(locus.ref_tag)
        # Seven substitutions plus one deletion reach the edit budget (8) only
        # through substitutions: the composition guard keeps this paralog-like
        # read off-target.
        for p in (12, 15, 18, 24, 27, 30, 33):
            ref[p] = "A" if ref[p] != "A" else "C"
        read = "".join(ref[:10] + ref[11:]) + "A"
        self.assertEqual(mp.select_marker(read, [0], [locus], 8, gap_band=10), (None, "offtarget"))

    def test_tail_gap_keeps_the_substitution_reading(self):
        # A shift inside the last 12 bases is not believed for classification:
        # clustered end-of-read errors would otherwise masquerade as a gap.
        ref = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"
        variant_pos = 39
        alt = ref[:variant_pos] + "T" + ref[variant_pos + 1:]  # ref[39] == 'A'
        locus = mp.Locus(0, "m1", ref[:8], "", "A", "T", variant_pos, ref, alt)
        # Deletion at 37 pushes ref[40]='C' into index 39 and appends a T.
        read = ref[:37] + ref[38:] + "T"
        self.assertEqual(read[variant_pos], "C")
        self.assertEqual(mp.classify_sequence(locus, read), "Other")  # fixed-offset reading, not the gapped one

    def test_denoise_prefers_closer_parent_before_deeper(self):
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        child = "A" * 6 + "G" + "A" * 3 + "C" + "A" * 9  # distance 1 from alt, 2 from ref
        _, mapping = mp.denoise_counts({ref: 1000, alt: 100, child: 3}, max_distance=2, alpha=2.0, protected=(ref, alt))
        self.assertEqual(mapping[child], alt)

    def _pooled_case(self, child, others_alt):
        ref = "A" * 20
        alt = "A" * 6 + "G" + "A" * 13
        locus = mp.Locus(0, "m1", "AAAAAA", "", "A", "G", 6, ref, alt)
        results = []
        for index in range(21):
            sample = mp.Sample(f"s{index}", f"s{index}", "unused.fastq", {})
            if index == 0:
                counts = {(0, ref): 25, (0, child): 25}
            else:
                counts = {(0, ref): 25, (0, alt): 25} if others_alt else {(0, ref): 25}
            results.append(mp.SampleResult(sample=sample, counts=counts,
                                           quality_support={k: (v, v, v, v) for k, v in counts.items()},
                                           partial_counts={}, partial_quality_support={}, qc={}))
        return locus, results

    def test_pooled_pass_never_moves_reads_across_alleles(self):
        # A private alt-bearing haplotype (one substitution beyond the alt
        # base) fails every evidence route; 0.9.x moved it onto the far deeper
        # Ref two edits away, turning a 25/25 heterozygote into Ref=50.
        alt_child = "A" * 6 + "G" + "A" * 3 + "C" + "A" * 9
        locus, results = self._pooled_case(alt_child, others_alt=True)
        finalized = mp.finalize_counts(results, [locus], mp.DenoiseSettings(1, 2, 2.0), rule(rep_samples=2, depth_reads=30))
        s0 = finalized[0].counts
        self.assertEqual(s0[(0, locus.ref_tag)], 25)
        self.assertEqual(s0.get((0, locus.alt_tag), 0), 25)  # folded into the same-allele parent, depth kept
        self.assertNotIn((0, alt_child), s0)

    def test_pooled_pass_folds_same_allele_child_and_keeps_depth(self):
        ref_child = "A" * 10 + "C" + "A" * 9  # ref-bearing, one substitution away from Ref
        locus, results = self._pooled_case(ref_child, others_alt=False)
        finalized = mp.finalize_counts(results, [locus], mp.DenoiseSettings(1, 1, 2.0), rule(rep_samples=2, depth_reads=30))
        self.assertEqual(finalized[0].counts[(0, locus.ref_tag)], 50)
        self.assertEqual(finalized[0].qc["reads_absorbed_global_denoise"], 25)

    def test_n_bearing_tags_are_credited_or_dropped_never_counted_as_haplotypes(self):
        ref = "ACGTACGTACGTACGTACGT"
        alt = ref[:6] + "T" + ref[7:]  # ref[6] is G
        locus = mp.Locus(0, "m1", ref[:6], "", "G", "T", 6, ref, alt)
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "s1.fastq"
            with fastq.open("w") as handle:
                records = [ref] * 30 + [ref[:10] + "N" + ref[11:]] * 5 + [ref[:6] + "N" + ref[7:]] * 4
                for number, sequence in enumerate(records):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            sample = mp.Sample("s1", "s1", str(fastq), {})
            result = mp.process_sample(
                sample, [locus], mp.PrimerIndex([locus]),
                mp.AssignmentSettings(tag_length=20, max_offtarget_mismatches=4, primer_mismatches=1, partial_min_length=30, gap_band=10),
            )
        self.assertEqual(result.counts[(0, ref)], 35)          # N away from the target: credited to Ref
        self.assertFalse(any("N" in key[1] for key in result.counts))
        self.assertEqual(result.qc["reads_n_mapped"], 5)
        self.assertEqual(result.qc["reads_n_dropped"], 4)      # N at the target fits Ref and Alt: ambiguous

    def test_leading_anchor_cuts_raw_reads_at_the_common_sequence(self):
        # Raw DArT deliveries keep an inline barcode and the common sequence in
        # front of the amplicon; with --leading-anchor the read starts where
        # that sequence ends, and reads without it are counted, not guessed.
        ref = "ACGTACGTACGTACGTACGT"
        alt = ref[:6] + "T" + ref[7:]
        locus = mp.Locus(0, "m1", ref[:6], "", "G", "T", 6, ref, alt)
        anchor = "GATGGATCCAGTGC"
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "s1.fastq"
            records = (["TTCATAAGTTGC" + anchor + ref + "CCGATCTGGG"] * 20   # barcode + anchor + tag + tail
                       + ["GGTAAGCAT" + anchor + alt + "CCGATC"] * 10          # a different barcode length
                       + [ref + "CCGATCTGGG"] * 5)                            # no anchor at all
            with fastq.open("w") as handle:
                for number, sequence in enumerate(records):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            sample = mp.Sample("s1", "s1", str(fastq), {})
            settings = mp.AssignmentSettings(tag_length=20, max_offtarget_mismatches=4, primer_mismatches=1,
                                             partial_min_length=30, gap_band=10, leading_anchor=anchor)
            result = mp.process_sample(sample, [locus], mp.PrimerIndex([locus]), settings)
        self.assertEqual(result.counts[(0, ref)], 20)
        self.assertEqual(result.counts[(0, alt)], 10)
        self.assertEqual(result.qc["reads_no_leading_anchor"], 5)
        self.assertEqual(result.qc["reads_assigned"], 30)
        self.assertIn("leading_anchor", settings.cache_arguments())

    def test_fastq_records_leading_anchor_cuts_and_skips(self):
        anchor = "GATGGATCCAGTGC"
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "s.fastq"
            with fastq.open("w") as handle:
                for number, sequence in enumerate(["TTCA" + anchor + "ACGTACGT", "ACGTACGT", "GG" + anchor + "TTTT"]):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            records = list(mp.fastq_records(fastq, anchor))
        self.assertEqual([seq for seq, _ in records], ["ACGTACGT", "TTTT"])
        self.assertTrue(all(len(seq) == len(qual) for seq, qual in records))

    def test_raw_delivery_is_detected_and_named(self):
        anchor = "GATGGATCCAGTGC"
        tag = "ACGTACGTACGTACGTACGT"
        with tempfile.TemporaryDirectory() as tmp:
            # One inline barcode per sample, constant within a file: the
            # anchor is what the files share, so the check samples several.
            raws = []
            for barcode in ("TTCATAAGTTGC", "GGTAAGCATGCA"):  # barcodes that end differently, so only the anchor is shared
                raw = Path(tmp) / f"raw_{barcode}.fastq"
                with raw.open("w") as handle:
                    for number in range(200):
                        handle.write(f"@r{number}\n{barcode}{anchor}{tag}CCGATC\n+\n{'I' * (len(barcode) + len(anchor) + len(tag) + 6)}\n")
                raws.append(raw)
            raw = raws
            clean = Path(tmp) / "clean.fastq"
            with clean.open("w") as handle:
                for number in range(200):
                    handle.write(f"@r{number}\n{tag}CCGATC\n+\n{'I' * (len(tag) + 6)}\n")
            prefixes = {tag[:mp.LEADING_ANCHOR_PREFIX]}
            report = mp.detect_leading_anchor(raw, prefixes)
            self.assertEqual(report["anchor"], anchor)
            self.assertEqual(report["anchored_rate"], 0.0)
            self.assertEqual(report["shifted_rate"], 1.0)
            self.assertEqual(mp.detect_leading_anchor([clean], prefixes)["anchor"], "")
            # no option: a raw delivery stops the run and names the anchor
            args = argparse.Namespace(leading_anchor=None)
            with self.assertRaises(mp.MadcError) as context:
                mp.resolve_leading_anchor(args, raw, prefixes, "generate")
            self.assertIn(f"--leading-anchor {anchor}", str(context.exception))
            # a preprocessed delivery passes silently
            self.assertEqual(mp.resolve_leading_anchor(argparse.Namespace(leading_anchor=None), [clean], prefixes, "generate"), "")
            # auto applies the detected anchor and records it
            args = argparse.Namespace(leading_anchor="auto")
            self.assertEqual(mp.resolve_leading_anchor(args, raw, prefixes, "generate"), anchor)
            self.assertEqual(args.leading_anchor, anchor)
            self.assertTrue(args.leading_anchor_detected)
            # an explicit value is used as given, without reading the FASTQ
            self.assertEqual(mp.resolve_leading_anchor(argparse.Namespace(leading_anchor="acgt"), raw, prefixes, "generate"), "ACGT")

    def test_denoise_neighbour_lookup_matches_the_scan(self):
        # The distance-1 path enumerates neighbours instead of scanning
        # every centroid; it must pick exactly the parents the scan picks.
        def scan(counts, alpha, protected):
            protected_set = set(protected); centroids = {}; mapping = {}
            for sequence, abundance in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
                if sequence in protected_set:
                    centroids[sequence] = centroids.get(sequence, 0) + abundance; mapping[sequence] = sequence; continue
                best = None
                for parent, pa in centroids.items():
                    if pa < abundance * 2 ** (alpha + 1): continue
                    d = mp.hamming_distance(sequence, parent, stop_after=1)
                    if d != 1: continue
                    cand = (1, -pa, parent)
                    if best is None or cand < best: best = cand
                if best is None: centroids[sequence] = abundance; mapping[sequence] = sequence
                else: centroids[best[2]] += abundance; mapping[sequence] = best[2]
            return centroids, mapping
        rng = random.Random(7)
        for trial in range(40):
            base = "".join(rng.choice("ACGT") for _ in range(12))
            counts = {}
            for _ in range(rng.randint(5, 60)):
                seq = list(base)
                for _ in range(rng.randint(0, 3)):
                    seq[rng.randrange(12)] = rng.choice("ACGTN")
                seq = "".join(seq)
                counts[seq] = counts.get(seq, 0) + rng.choice([1, 1, 1, 2, 3, 5, 9, 40, 300])
            protected = [base] if rng.random() < 0.5 else []
            alpha = rng.choice([1.0, 2.0])
            self.assertEqual(mp.denoise_counts(counts, 1, alpha, protected), scan(counts, alpha, protected))

    def test_intra_target_repeated_kmers_are_not_anchors(self):
        unit = "ACGTTGCAAGCTTGACCGATA"  # 21 bp, repeated twice inside the target
        target = unit + "GGCTTACGGATCAGTTCGAC[A/G]TTGCCATGGAATCCGATTAC" + unit + "TTGACCA"
        with tempfile.TemporaryDirectory() as tmp:
            design = Path(tmp) / "design.tsv"
            design.write_text("MarkerName\tTargetSequence\n" + f"m1\t{target}\n")
            mp._normalize_init(str(design), 21)
            index = mp._NORM_STATE["index"]
        self.assertNotIn(unit, index)
        self.assertIn("GGCTTACGGATCAGTTCGACA", index)

    def test_merge_keeps_r2_extension_with_unequal_read_lengths(self):
        import random as _random
        rng = _random.Random(3)
        fragment = "".join(rng.choice("ACGT") for _ in range(70))
        r1 = fragment[:60]
        r2 = mp.reverse_complement("GGGGG" + fragment)  # R2 read through the fragment start into 5 adapter bases
        merged = mp.merge_read_pair(r1, "I" * len(r1), r2, "I" * len(r2), 15, 0.1)
        self.assertIsNotNone(merged)
        self.assertEqual(merged[0], fragment)

    def test_merge_pairs_rejects_mates_out_of_step(self):
        with tempfile.TemporaryDirectory() as tmp:
            r1 = Path(tmp) / "a_R1.fastq"
            r2 = Path(tmp) / "a_R2.fastq"
            seq = "ACGT" * 10
            r1.write_text(f"@read1/1\n{seq}\n+\n{'I' * 40}\n")
            r2.write_text(f"@read2/2\n{mp.reverse_complement(seq)}\n+\n{'I' * 40}\n")
            task = (str(r1), str(r2), str(Path(tmp) / "a.fastq.gz"), "a", 15, 0.1, "drop")
            with self.assertRaises(mp.MadcError):
                mp._merge_pairs_sample(task)
            r2.write_text(f"@read1/2\n{mp.reverse_complement(seq)}\n+\n{'I' * 40}\n")
            sample, stats = mp._merge_pairs_sample(task)
            self.assertEqual(stats["pairs_merged"], 1)

    def test_scan_subset_is_seeded_or_sorted_prefix(self):
        paths = [Path(f"s{index:02d}.fastq") for index in range(10)]
        self.assertEqual(mp.select_scan_subset(paths, 3, None), paths[:3])
        drawn = mp.select_scan_subset(paths, 3, 7)
        self.assertEqual(len(drawn), 3)
        self.assertEqual(drawn, mp.select_scan_subset(paths, 3, 7))
        self.assertEqual(drawn, sorted(drawn))


class ReviewRoundTwoTests(unittest.TestCase):
    def test_n_resolution_never_credits_an_ambiguity_coded_canonical(self):
        ref = "ACGTACGTACRTACGTACGT"  # R at index 10
        alt = ref[:6] + "T" + ref[7:]
        locus = mp.Locus(0, "m1", ref[:6], "", "G", "T", 6, ref, alt)
        n_tag = ref[:3] + "N" + ref[4:]  # compatible with the R-coded canonical, no concrete sequence observed
        counts = {(0, n_tag): 5}
        support = {(0, n_tag): [5, 5, 5, 5]}
        qc = mp.Counter()
        mp.resolve_n_tags(counts, support, [locus], qc)
        self.assertNotIn((0, ref), counts)              # the template acquired no invented counts
        self.assertEqual(qc["reads_n_dropped"], 5)
        # With the concrete realisation observed, the reads go there — no false ambiguity.
        concrete = ref.replace("R", "A")
        counts = {(0, concrete): 30, (0, n_tag): 5}
        support = {(0, concrete): [30] * 4, (0, n_tag): [5] * 4}
        qc = mp.Counter()
        mp.resolve_n_tags(counts, support, [locus], qc)
        self.assertEqual(counts[(0, concrete)], 35)
        self.assertEqual(qc["reads_n_mapped"], 5)

    def test_gapped_scoring_honours_design_iupac(self):
        ref = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"
        design = ref[:25] + "R" + ref[26:]      # ref[25] is 'T' in the concrete read below? make it A/G-compatible
        read = ref[:10] + ref[11:25] + "A" + ref[26:] + "C"  # one deletion at 10, A at the R position
        edits, subs, gaps, events, _, _ = mp.gapped_stats(read, design, 10)
        self.assertEqual((edits, subs, gaps, events), (1, 0, 1, 1))
        self.assertLessEqual(mp.gapped_prescreen(read, design, 10), 1)
        # an observed N still counts as a mismatch in the gapped path
        n_read = read[:30] + "N" + read[31:]
        self.assertEqual(mp.gapped_stats(n_read, design, 10)[1], 1)

    def test_scan_churn_counts_dropped_n_reads(self):
        ra = {"reads_moved_by_sample_denoise": 3, "reads_moved_by_global_denoise": 2,
              "reads_dropped_failed_evidence": 5, "reads_n_dropped": 200}
        self.assertEqual(mp.scan_churn(ra), 210)
        del ra["reads_n_dropped"]
        self.assertEqual(mp.scan_churn(ra), 10)

    def test_fastq_records_named_validates_lengths(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = Path(tmp) / "bad.fastq"
            bad.write_text("@r1/1\nACGTACGTAC\n+\nIIIIIIIII\n")  # 10 bases, 9 qualities
            with self.assertRaises(mp.MadcError):
                list(mp.fastq_records_named(bad))


class MarkerTypeTests(unittest.TestCase):
    def test_marker_type_labels(self):
        ref = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"
        alt = ref[:20] + "C" + ref[21:]
        locus = mp.Locus(0, "m1", ref[:8], "", ref[20], "C", 20, ref, alt)
        self.assertEqual(mp.marker_type_label(locus, ref), ".")
        self.assertEqual(mp.marker_type_label(locus, alt), "snp")
        self.assertEqual(mp.marker_type_label(locus, ref[:30] + ("G" if ref[30] != "G" else "T") + ref[31:]), "snp")
        self.assertEqual(mp.marker_type_label(locus, ref[:10] + ref[11:] + "A"), "indel")
        self.assertEqual(mp.marker_type_label(locus, alt[:10] + alt[11:] + "A"), "indel")  # nearest canonical is Alt
        both = ref[:10] + ref[11:30] + "T" + ref[31:] + "A"
        self.assertEqual(mp.marker_type_label(locus, both[:41]), "both")
        designed_indel = mp.Locus(0, "m2", ref[:8], "", "T", "TG", 20, ref, ref[:20] + "G" + ref[20:40])
        self.assertEqual(mp.marker_type_label(designed_indel, designed_indel.alt_tag), "indel")

    def test_marker_type_column_written_and_switchable(self):
        ref = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"
        alt = ref[:20] + "C" + ref[21:]
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "s1.fastq"
            with fastq.open("w") as handle:
                for number, sequence in enumerate([ref] * 30 + [alt] * 20 + [ref[:10] + ref[11:] + "A"] * 40):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * 41}\n")
            tags = Path(tmp) / "tags.csv"
            tags.write_text("MarkerName,LeftPrimerSequence,RightPrimerSequence,RefSequence,AltSequence,VariantRef,VariantAlt,VariantPositionInTagSequence\n"
                            f"m1,{ref[:8]},,{ref},{alt},{ref[20]},C,21\n")
            for mode, header in (("on", "MarkerType"), ("off", "AvgPIC")):
                out = Path(tmp) / f"{mode}.csv"
                self.assertEqual(mp.main(["generate", "--fastq-dir", tmp, "--fastq-glob", "*.fastq", "--reference-tags", str(tags),
                                          "--tag-length", "41", "--jobs", "1", "--depth-min-reads", "10", "--annotations", mode,
                                          "--output", str(out)]), 0)
                with out.open() as handle:
                    rows = list(csv.reader(handle))
                self.assertEqual(rows[0][12], header)
                self.assertEqual(rows[0][11], "Variants" if mode == "on" else "PICSnp")
                labels = {row[0].split("|")[-1].rstrip("0123456789_"): row[12] for row in rows[1:]}
                if mode == "on":
                    self.assertEqual(labels["Ref"], "."); self.assertEqual(labels["Alt"], "snp"); self.assertEqual(labels["RefMatch"], "indel")
                else:
                    self.assertEqual(set(labels.values()), {"."})


class VariantsTests(unittest.TestCase):
    REF = "ACGTTGCAAGCTTGACCGATAGGCTTACGGATCAGTTCGAC"  # 41 bp; target index 20 ('A')

    def _locus(self):
        ref = self.REF
        alt = ref[:20] + "C" + ref[21:]
        return mp.Locus(0, "m1", ref[:8], "", "A", "C", 20, ref, alt)

    def _cell(self, locus, sequence):
        return mp.format_variants(mp.variant_events(locus, sequence))

    def test_tokens_for_snp_insertion_and_deletion(self):
        locus = self._locus()
        ref = locus.ref_tag
        self.assertEqual(self._cell(locus, ref), ".")
        self.assertEqual(self._cell(locus, locus.alt_tag), "21:A>C")
        two = ref[:5] + "A" + ref[6:20] + "C" + ref[21:]  # ref[5] is G
        self.assertEqual(self._cell(locus, two), "6:G>A;21:A>C")
        # a 2-bp deletion of ref[10:12] pulls two downstream bases into the window
        deletion = ref[:10] + ref[12:] + "AA"
        self.assertEqual(self._cell(locus, deletion), "10:GCT>G;overhang=2")
        self.assertEqual(mp.variant_events(locus, deletion)[-1], ("overhang", 42, ".", "AA"))
        # a 2-bp insertion after ref[9] pushes the last two Ref bases out; inside the GG run it
        # left-normalises onto the base before the run
        insertion = ref[:10] + "GG" + ref[10:39]
        self.assertEqual(self._cell(locus, insertion), "9:A>AGG;clip=2")
        self.assertEqual(mp.variant_events(locus, insertion)[-1], ("clip", 40, "AC", "."))
        # an inserted block whose bases match the reference (a tandem duplication) is one event,
        # not fragments: fewest gap openings among equal-edit explanations
        duplication = ref[:20] + ref[12:20] + ref[20:33]
        events = mp.variant_events(locus, duplication)
        self.assertEqual([event[0] for event in events], ["ins", "clip"])
        self.assertEqual(len(events[0][3]) - len(events[0][2]), 8)
        self.assertEqual(events[1], ("clip", 34, ref[33:], "."))
        # a row nearest to Alt is read the way its label is (gapped), then described against Ref,
        # so the designed substitution appears alongside the gap
        alt_deletion = locus.alt_tag[:10] + locus.alt_tag[11:] + "A"
        self.assertEqual(mp.marker_type_label(locus, alt_deletion), "indel")
        self.assertEqual(self._cell(locus, alt_deletion), "10:GC>G;21:A>C;overhang=1")
        # a marker whose Alt is a designed insertion: every row is described in the Ref frame, so an
        # AltMatch with one extra SNP lists the inherited insertion plus the SNP and is labelled 'snp'
        ins_alt = ref[:20] + "G" + ref[20:40]
        designed = mp.Locus(0, "m2", ref[:8], "", "T", "TG", 20, ref, ins_alt)
        self.assertEqual(mp.format_variants(mp.variant_events(designed, ins_alt)), "20:T>TG;clip=1")
        self.assertEqual(mp.marker_type_label(designed, ins_alt), "indel")
        altmatch = ins_alt[:30] + "A" + ins_alt[31:]  # ins_alt[30] is ref[29] == "G"
        self.assertEqual(mp.format_variants(mp.variant_events(designed, altmatch)), "20:T>TG;30:G>A;clip=1")
        self.assertEqual(mp.marker_type_label(designed, altmatch), "snp")
        self.assertEqual(mp.classify_sequence(designed, altmatch), "AltMatch")
        # canonical tags sharing an ambiguity code: the shared column is not a substitution anywhere
        shared_ref = ref[:12] + "R" + ref[13:]
        shared_alt = shared_ref[:20] + "G" + shared_ref[20:40]
        shared = mp.Locus(0, "m3", shared_ref[:8], "", "T", "TG", 20, shared_ref, shared_alt)
        self.assertEqual(mp.format_variants(mp.variant_events(shared, shared_alt)), "20:T>TG;clip=1")
        self.assertEqual(mp.marker_type_label(shared, shared_alt), "indel")
        self.assertEqual(mp.gapped_stats(shared_alt, shared_ref, 10)[1], 0)
        self.assertLessEqual(mp.gapped_prescreen(shared_alt, shared_ref, 10), 1)
        # gap plus substitution: both are listed, in tag order
        both = ref[:10] + ref[12:30] + "T" + ref[31:] + "AA"
        self.assertEqual(self._cell(locus, both), "10:GCT>G;31:A>T;overhang=2")
        self.assertEqual(mp.marker_type_label(locus, both), "both")

    def test_left_normalisation_inside_repeats(self):
        locus = self._locus()
        ref = locus.ref_tag  # ref[11:13] == "TT"
        one_t_deleted = ref[:12] + ref[13:] + "A"
        self.assertEqual(one_t_deleted, ref[:11] + ref[12:] + "A")  # the two placements are one string
        self.assertEqual(self._cell(locus, one_t_deleted), "11:CT>C;overhang=1")
        one_t_inserted = ref[:12] + "T" + ref[12:40]
        self.assertEqual(self._cell(locus, one_t_inserted), "11:C>CT;clip=1")
        # a gap beside a substitution takes the leftmost equal-cost placement and never merges with
        # it: ref[21:23] is "GG"; reading "T" there is one deletion (anchored on ref[20]) plus one
        # substitution at ref[22], and the substituted base blocks any further shift
        beside = ref[:21] + "T" + ref[23:] + "A"
        self.assertEqual(self._cell(locus, beside), "21:AG>A;23:G>T;overhang=1")

    def test_alignment_path_matches_stats(self):
        rng = random.Random(7)
        for _ in range(300):
            ref = "".join(rng.choice("ACGT") for _ in range(60))
            read = list(ref)
            for _ in range(rng.randint(0, 3)):
                op = rng.random()
                index = rng.randrange(len(read))
                if op < 0.4:
                    read[index] = rng.choice("ACGT")
                elif op < 0.7:
                    del read[index]
                else:
                    read.insert(index, rng.choice("ACGT"))
            read = "".join(read)[:60].ljust(60, "A")
            edits, subs, gap_bases, gap_events, _, _ = mp.gapped_stats(read, ref, 6)
            ops, end_read, end_ref, path_edits = mp.gapped_alignment(read, ref, 6)
            ops_subs = sum(1 for kind, i, j in ops if kind == "M" and read[i] != ref[j])
            ops_gaps = sum(1 for kind, _, _ in ops if kind != "M")
            ops_opens = sum(1 for k, op in enumerate(ops) if op[0] != "M" and (k == 0 or ops[k - 1][0] != op[0]))
            self.assertEqual(path_edits, edits)  # same banded edit distance as the rescue decision
            self.assertEqual(ops_subs + ops_gaps, edits)  # unaligned tails are free
            self.assertLessEqual(ops_opens, gap_events + 1)  # never more openings than the unit-cost path (+ a leading gap it does not count)
            normalised = mp._left_align(ops, read, ref)
            self.assertEqual(sum(1 for kind, _, _ in normalised if kind != "M"), ops_gaps)
            self.assertEqual(sum(1 for kind, i, j in normalised if kind == "M" and read[i] != ref[j]), ops_subs)
            # the normalised path still spells the read: every read base appears exactly once
            self.assertEqual(sorted(i for kind, i, _ in normalised if kind != "D"), list(range(end_read)))
            self.assertEqual(sorted(j for kind, _, j in normalised if kind != "I"), list(range(end_ref)))

    def test_variants_column_and_sidecar(self):
        ref = self.REF
        alt = ref[:20] + "C" + ref[21:]
        deletion = ref[:10] + ref[11:] + "A"  # deletes ref[10] (C) after G: no shift
        with tempfile.TemporaryDirectory() as tmp:
            fastq = Path(tmp) / "s1.fastq"
            with fastq.open("w") as handle:
                for number, sequence in enumerate([ref] * 30 + [alt] * 20 + [deletion] * 40):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * 41}\n")
            tags = Path(tmp) / "tags.csv"
            tags.write_text("MarkerName,LeftPrimerSequence,RightPrimerSequence,RefSequence,AltSequence,VariantRef,VariantAlt,VariantPositionInTagSequence\n"
                            f"m1,{ref[:8]},,{ref},{alt},{ref[20]},C,21\n")
            out = Path(tmp) / "run.csv"
            for mode in ("on", "off"):  # the same output path: the second run must retire the first run's sidecar
                self.assertEqual(mp.main(["generate", "--fastq-dir", tmp, "--fastq-glob", "*.fastq", "--reference-tags", str(tags),
                                          "--tag-length", "41", "--jobs", "1", "--depth-min-reads", "10", "--annotations", mode,
                                          "--output", str(out)]), 0)
                with out.open() as handle:
                    rows = list(csv.reader(handle))
                sidecar = Path(str(out) + ".variants.csv")
                with (Path(str(out) + ".qc.json")).open() as handle:
                    qc = json.load(handle)
                cells = {row[0].split("|")[-1].rstrip("0123456789_"): row[11] for row in rows[1:]}
                if mode == "on":
                    self.assertEqual(rows[0][11], "Variants")
                    self.assertEqual(rows[0][12], "MarkerType")
                    self.assertEqual(cells["Ref"], ".")
                    self.assertEqual(cells["Alt"], "21:A>C")
                    self.assertEqual(cells["RefMatch"], "10:GC>G;overhang=1")
                    with sidecar.open() as handle:
                        lines = list(csv.reader(handle))
                    self.assertEqual(lines[0], mp.VARIANTS_SIDECAR_COLUMNS)
                    self.assertIn(["m1", "m1|Alt", alt, "21", "A", "C", "snp"], lines)
                    self.assertIn(["m1", "m1|RefMatch", deletion, "10", "GC", "G", "del"], lines)
                    self.assertIn(["m1", "m1|RefMatch", deletion, "42", ".", "A", "overhang"], lines)
                    summary = qc["output_summary"]["variants"]
                    self.assertEqual(summary["rows_with_events"], 2)
                    self.assertEqual(summary["events"], {"snp": 1, "del": 1, "overhang": 1})
                    self.assertEqual(summary["sidecar"], str(sidecar))
                    self.assertEqual(qc["parameters"]["annotations"], "on")
                else:
                    self.assertEqual(rows[0][11], "PICSnp")
                    self.assertEqual(rows[0][12], "AvgPIC")
                    self.assertEqual(set(cells.values()), {"."})
                    self.assertFalse(sidecar.exists())
                    self.assertIsNone(qc["output_summary"]["variants"]["sidecar"])


class TileFamilyTests(unittest.TestCase):
    """Tiles: two target SNPs inside one amplicon, read through one shared tag window."""

    # SNPs at 10 (C/T) and 20 (G/A); both tiles start at the same primer.
    AMPLICON = "ACGTACGGTA" + "C" + "TTGACCAGT" + "G" + "ATCCGTAGC"

    def tiles(self, groups=("", ""), iupac=True):
        """Tile A (SNP at 10) and tile B (SNP at 20); with ``iupac`` each tag carries
        the other tile's SNP as an ambiguity code, as DArTag panels write them."""
        amplicon = self.AMPLICON
        ref_a = amplicon[:20] + ("R" if iupac else "G") + amplicon[21:]
        ref_b = amplicon[:10] + ("Y" if iupac else "C") + amplicon[11:]
        primer = amplicon[:8]
        return [
            mp.Locus(0, "tileA", primer, "", "C", "T", 10, ref_a, ref_a[:10] + "T" + ref_a[11:], locus_group=groups[0]),
            mp.Locus(1, "tileB", primer, "", "G", "A", 20, ref_b, ref_b[:20] + "A" + ref_b[21:], locus_group=groups[1]),
        ]

    def test_tied_read_is_credited_to_every_tile(self):
        loci = self.tiles()
        self.assertEqual(mp.tile_families(loci), ((0, 1),))
        family_of = {0: (0, 1), 1: (0, 1)}
        for read in (self.AMPLICON, self.AMPLICON[:10] + "T" + self.AMPLICON[11:]):
            self.assertEqual(mp.select_markers(read, [0, 1], loci, 3, family_of=family_of), ((0, 1), "assigned"))
            # Without the family the read ties and is dropped, as before.
            self.assertEqual(mp.select_marker(read, [0, 1], loci, 3), (None, "ambiguous"))

    def test_untied_read_reaches_every_tile_it_fits(self):
        # Tiles with identical tags and no ambiguity codes: a read carrying
        # tile A's Alt fits A exactly and B with one substitution.  It no longer
        # goes to A alone; B sees it too, as it would in a panel of its own.
        loci = self.tiles(iupac=False)
        read = self.AMPLICON[:10] + "T" + self.AMPLICON[11:]
        self.assertEqual(mp.select_marker(read, [0, 1], loci, 3), (0, "assigned"))
        families = mp.tile_families(loci)
        family_of = {member: family for family in families for member in family}
        self.assertEqual(mp.select_markers(read, [0, 1], loci, 3, family_of=family_of), ((0, 1), "assigned"))

    def test_tile_group_modes(self):
        self.assertEqual(mp.tile_families(self.tiles(), "declared"), ())
        self.assertEqual(mp.tile_families(self.tiles(), "off"), ())
        self.assertEqual(mp.tile_families(self.tiles(("amp1", "amp1")), "declared"), ((0, 1),))
        self.assertEqual(mp.tile_families(self.tiles(("amp1", "")), "auto"), ((0, 1),))
        # Markers declared in different groups are never linked.
        self.assertEqual(mp.tile_families(self.tiles(("amp1", "amp2")), "auto"), ())
        # Different windows of one LocusGroup (normalize-reads tiling) are not one tag window.
        window_a, window_b = self.tiles(("amp1", "amp1"))
        shifted = mp.Locus(1, "tileB", self.AMPLICON[5:13], "", "G", "A", 15, self.AMPLICON[5:] + "AAAAA",
                           self.AMPLICON[5:20] + "A" + self.AMPLICON[21:] + "AAAAA", locus_group="amp1")
        self.assertEqual(mp.tile_families([window_a, shifted]), ())

    def test_tie_across_families_stays_ambiguous(self):
        tile_a, tile_b = self.tiles(("amp1", "amp1"))
        elsewhere = dataclasses.replace(tile_a, index=2, marker="elsewhere", locus_group="amp2")
        loci = [tile_a, tile_b, elsewhere]
        families = mp.tile_families(loci)
        self.assertEqual(families, ((0, 1),))
        family_of = {member: family for family in families for member in family}
        self.assertEqual(mp.select_markers(self.AMPLICON, [0, 1, 2], loci, 3, family_of=family_of), ((), "ambiguous"))

    def test_generate_counts_tiles_in_one_panel(self):
        amplicon = self.AMPLICON
        haplotypes = [amplicon] * 12 + [amplicon[:10] + "T" + amplicon[11:]] * 6 + [amplicon[:20] + "A" + amplicon[21:]] * 4
        partials = [amplicon[:26]] * 2
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fastq_dir = root / "fastq"
            fastq_dir.mkdir()
            with (fastq_dir / "s1.fastq").open("w") as handle:
                for number, sequence in enumerate(haplotypes + partials):
                    handle.write(f"@r{number}\n{sequence}\n+\n{'I' * len(sequence)}\n")
            tags = root / "tags.csv"
            mp.write_reference_tags(tags, self.tiles())  # no LocusGroup: the tiles are detected
            outcomes = {}
            for mode in ("auto", "off"):
                output = root / f"{mode}.csv"
                self.assertEqual(mp.main(["generate", "--fastq-dir", str(fastq_dir), "--reference-tags", str(tags),
                                          "--tag-length", "30", "--jobs", "1", "--reproducible-min-samples", "0",
                                          "--depth-min-reads", "1", "--tile-groups", mode,
                                          "--cache-dir", str(root / "cache"), "--output", str(output)]), 0)
                totals = mp.Counter()
                for record in mp.read_madc(output)["records"]:
                    totals[record["clone_id"]] += record["counts"][0]
                with Path(str(output) + ".qc.json").open() as handle:
                    outcomes[mode] = (totals, json.load(handle))
            totals, qc = outcomes["auto"]
            self.assertEqual(totals, {"tileA": 24, "tileB": 24})
            self.assertEqual(qc["tile_families"], [["tileA", "tileB"]])
            accounting = qc["read_accounting"]
            self.assertEqual(accounting["reads_assigned_full_length"], 22)   # each read once
            self.assertEqual(accounting["tile_extra_credits_full_length"], 22)
            self.assertEqual(accounting["tile_extra_credits_partial"], 2)
            self.assertEqual(
                accounting["reads_in_final_count_table"],
                accounting["reads_assigned_full_length"] + accounting["tile_extra_credits_full_length"]
                - accounting["reads_n_dropped"] + accounting["partial_reads_added_to_counts"],
            )
            totals, qc = outcomes["off"]
            self.assertEqual(totals, {"tileA": 0, "tileB": 0})
            self.assertEqual(qc["tile_families"], [])
            self.assertEqual(qc["sample_qc"][0]["reads_ambiguous"], 22)
            # The families are part of the cache key: two runs, two cache entries.
            self.assertEqual(len(list((root / "cache").glob("s1.*.json.gz"))), 2)
