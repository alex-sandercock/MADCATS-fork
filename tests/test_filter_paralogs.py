"""filter-paralogs (madcats_paralogs.py): genome labels, placement, MADC
rewriting, and the optional mappy aligner.

Most tests use a fake aligner with preset hits, so they run without mappy.
GenomeAlignerTests builds a tiny synthetic genome and runs the real command;
they are skipped when mappy is not installed.
"""
import contextlib
import csv
import importlib.util
import io
import json
import random
import subprocess
import sys
import tempfile
import unittest
from collections import Counter
from pathlib import Path
from unittest import mock

import madcats as mp
import madcats_paralogs as fp

HAS_MAPPY = importlib.util.find_spec("mappy") is not None
REPOSITORY = Path(__file__).resolve().parents[1]
SAMPLES = ["S1", "S2", "S3"]


def hit(contig, start, end, score=162, strand=1):
    return fp.GenomeHit(contig, start, end, strand, score)


class FakeAligner:
    """Stands in for GenomeAligner: genome-wide hits are given up front, and the
    window realignment scans the window ungapped on both strands."""

    def __init__(self, genome, hits):
        self.genome = genome
        self.contigs = frozenset(genome)
        self._hits = hits
        self.windows = []  # every window realigned against, in order

    def new_buffer(self):
        return None

    def hits(self, sequence, buffer):
        return list(self._hits.get(sequence, ()))

    def fetch(self, contig, start, end):
        return self.genome[contig][start:end]

    def window(self, contig, start, end):
        self.windows.append((contig, start, end))
        return FakeWindow(contig, start, self.fetch(contig, start, end))


class FakeWindow:
    def __init__(self, contig, start, sequence):
        self.contig, self.start, self.sequence = contig, start, sequence

    def hits(self, sequence, buffer=None):
        found = []
        for strand, query in ((1, sequence), (-1, mp.reverse_complement(sequence))):
            for offset in range(len(self.sequence) - len(query) + 1):
                mismatches = sum(a != b for a, b in zip(query, self.sequence[offset:]))
                if mismatches <= 6:
                    score = 2 * (len(query) - mismatches) - 8 * mismatches
                    found.append(fp.GenomeHit(self.contig, self.start + offset, self.start + offset + len(query), strand, score))
        return found


def madc_row(clone, category, sequence, counts, vendor=False):
    stats = ["0.9", "0.5", "0.5", "0.2", "0.2", "0.6", "0.37", "0.4"] if vendor else ["."] * 8
    return [f"{clone}|{category}", clone, sequence, sequence, "1"] + stats + ["10", "10", "."] + [str(c) for c in counts]


def write_madc(path, rows, vendor=False, preamble=()):
    header = list(mp.MADC_METADATA_COLUMNS)
    if not vendor:
        header = [mp.ANNOTATION_COLUMNS.get(column, column) for column in header]
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerows(preamble)
        writer.writerow(header + SAMPLES)
        writer.writerows(rows)


def read_rows(path):
    with path.open(newline="") as handle:
        return list(csv.reader(handle))


def read_dicts(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def mutate(sequence, position, base):
    return sequence[:position] + base + sequence[position + 1 :]


def changed(sequence, *positions):
    for position in positions:
        sequence = mutate(sequence, position, "A" if sequence[position] != "A" else "C")
    return sequence


class HitScoreTests(unittest.TestCase):
    def test_score_matches_minimap2_short_read_scoring(self):
        # (CIGAR, NM, matching bases, block length) as mappy reports them.
        self.assertEqual(fp.hit_score([[81, 0]], 0, 81, 81), 162)
        self.assertEqual(fp.hit_score([[81, 0]], 1, 80, 81), 80 * 2 - 8)
        # An ambiguous base is outside the block, inside NM, and costs 1.
        self.assertEqual(fp.hit_score([[81, 0]], 1, 80, 80), 80 * 2 - 1)
        # A 1-base insertion costs 12 + 2; the inserted base is not a mismatch.
        self.assertEqual(fp.hit_score([[40, 0], [1, 1], [40, 0]], 1, 80, 81), 160 - 14)
        # A long deletion takes the cheaper long-gap piece: min(12 + 2*20, 24 + 20) = 44.
        self.assertEqual(fp.hit_score([[40, 0], [20, 2], [41, 0]], 20, 81, 101), 162 - 44)

    def test_summarize_keeps_distinct_sites_best_first(self):
        placement = fp.summarize_hits([
            hit("chrB", 500, 581), hit("chrA", 100, 181), hit("chrA", 103, 184), hit("chrC", 0, 81, score=150),
        ])
        self.assertEqual([site.contig for site in placement.best], ["chrA", "chrB"])
        self.assertEqual([site.score for site in placement.sites], [162, 162, 150])
        self.assertEqual(fp.summarize_hits([]).sites, ())

    def test_adjacent_tandem_copies_stay_distinct(self):
        # 81-base tags on adjacent copies of an 80-base repeat share one base: two loci, not one.
        ref = fp.summarize_hits([hit("chr1", 500, 581), hit("chr1", 580, 661)], keep=None)
        self.assertEqual(len(ref.sites), 2)
        marker = fp.place_marker("chr1_540", ref, {"chr1"}, 81)
        self.assertEqual((marker.source, marker.ref_tag_sites, marker.sites), ("ref_tag+clone_id", 2, ref.sites[:1]))
        variant = fp.summarize_hits([hit("chr1", 500, 581, score=152), hit("chr1", 580, 661, score=152)])
        self.assertEqual(fp.genome_label(variant, marker).label, "ambiguous")


class PlacementTests(unittest.TestCase):
    contigs = {"Chr1", "Chr2", "scaffold_12"}

    def test_clone_position_parsing(self):
        self.assertEqual(fp.clone_position("Chr1_000333995", self.contigs), (("Chr1", 333995), "agree"))
        self.assertEqual(fp.clone_position("scaffold_12_000345", self.contigs), (("scaffold_12", 345), "agree"))
        self.assertEqual(fp.clone_position("marker7", self.contigs), (None, "no_position"))
        self.assertEqual(fp.clone_position("Chr9_100", self.contigs), (None, "contig_absent"))

    def test_unique_ref_tag_is_cross_checked_against_the_clone_id(self):
        ref = fp.GenomePlacement((hit("Chr1", 1000, 1081),))
        placed = fp.place_marker("Chr1_1060", ref, self.contigs, 81)
        self.assertEqual((placed.source, placed.clone_check, placed.sites), ("ref_tag", "agree", ref.best))

    def test_disagreeing_clone_id_is_checked_in_its_own_window(self):
        # The Ref tag's best genome-wide site is on Chr2 (the search missed the design site, or this assembly
        # carries the Alt allele there and a Ref-like copy on Chr2); the CloneID names Chr1:1043.
        ref = fp.GenomePlacement((hit("Chr2", 7000, 7081, score=96),))
        windows = []

        def align_ref_tag(*window):
            windows.append(window)
            return [hit("Chr1", 900, 981), hit("Chr1", 1000, 1081, score=156)]

        placed = fp.place_marker("Chr1_1043", ref, self.contigs, 81, align_ref_tag=align_ref_tag)
        self.assertEqual(windows, [("Chr1", 1043 - 1 - 81 - fp.SITE_PAD, 1043 + 81 + fp.SITE_PAD)])
        self.assertEqual((placed.source, placed.ref_tag_sites, placed.clone_check), ("clone_id_site", 1, "disagree"))
        self.assertEqual(placed.sites, (hit("Chr1", 1000, 1081, score=156),))  # the alignment covering the position
        # A multi-copy Ref tag none of whose sites covers the position is checked the same way.
        copies = fp.GenomePlacement((hit("Chr2", 7000, 7081), hit("Chr2", 9000, 9081)))
        tied = fp.place_marker("Chr1_1043", copies, self.contigs, 81, align_ref_tag=align_ref_tag)
        self.assertEqual((tied.source, tied.ref_tag_sites, tied.sites[0].contig), ("clone_id_site", 2, "Chr1"))
        # Where the Ref tag does not align over the CloneID position, the marker is unplaced, so it keeps its rows.
        lost = fp.place_marker("Chr1_1043", ref, self.contigs, 81, align_ref_tag=lambda *window: [hit("Chr1", 900, 981)])
        self.assertEqual((lost.source, lost.sites, lost.ref_tag_sites, lost.clone_check), ("unplaced", (), 1, "disagree"))

    def test_multicopy_ref_tag_resolved_by_clone_id_or_left_unresolved(self):
        ref = fp.GenomePlacement((hit("Chr1", 1000, 1081), hit("Chr2", 5000, 5081)))
        resolved = fp.place_marker("Chr2_5040", ref, self.contigs, 81)
        self.assertEqual((resolved.source, resolved.ref_tag_sites), ("ref_tag+clone_id", 2))
        self.assertEqual([site.contig for site in resolved.sites], ["Chr2"])
        unresolved = fp.place_marker("marker7", ref, self.contigs, 81)
        self.assertEqual((unresolved.source, len(unresolved.sites)), ("ref_tag_multicopy", 2))

    def test_unaligned_ref_tag_falls_back_to_the_clone_id_window(self):
        fallback = fp.place_marker("Chr1_2000", fp.GenomePlacement(()), self.contigs, 81)
        self.assertEqual((fallback.source, fallback.clone_check), ("clone_id", "fallback"))
        self.assertEqual(fallback.sites[0][:3], ("Chr1", 2000 - 1 - 81, 2000 + 81))
        self.assertEqual(fp.place_marker("marker7", None, self.contigs, 81).source, "unplaced")

    def test_clone_positions_can_be_ignored(self):
        def place(clone, *hits):
            return fp.place_marker(clone, fp.GenomePlacement(hits), self.contigs, 81, use_clone_position=False)

        unique = place("Chr2_1060", hit("Chr1", 1000, 1081))
        self.assertEqual((unique.source, unique.clone_check), ("ref_tag", "ignored"))
        tied = place("Chr2_5040", hit("Chr1", 1000, 1081), hit("Chr2", 5000, 5081))
        self.assertEqual((tied.source, len(tied.sites)), ("ref_tag_multicopy", 2))  # no CloneID tie-break
        unaligned = place("Chr1_2000")
        self.assertEqual((unaligned.source, unaligned.sites, unaligned.clone_check), ("unplaced", (), "ignored"))

    def test_labels(self):
        marker = fp.place_marker("Chr1_1060", fp.GenomePlacement((hit("Chr1", 1000, 1081),)), self.contigs, 81)

        def label(*hits, site_score=None):
            return fp.genome_label(fp.summarize_hits(hits), marker, site_score).label

        self.assertEqual(label(hit("Chr1", 1000, 1081)), "on_target")
        self.assertEqual(label(hit("Chr1", 1004, 1081, score=150)), "on_target")  # clipped, same locus
        self.assertEqual(label(hit("Chr2", 7000, 7081), hit("Chr1", 1000, 1081, score=150)), "paralog")
        self.assertEqual(label(hit("Chr1", 1081, 1162)), "paralog")  # the adjacent tandem copy
        self.assertEqual(label(hit("Chr2", 7000, 7081), hit("Chr1", 1000, 1081)), "ambiguous")
        self.assertEqual(label(hit("Chr2", 7000, 7081), hit("Chr2", 9000, 9081)), "paralog")  # all ties elsewhere
        self.assertEqual(label(), "unmapped")
        # A diverged allele the genome-wide search missed, scored by the site realignment.
        self.assertEqual(label(hit("Chr2", 7000, 7081, score=112), site_score=132), "on_target")
        self.assertEqual(label(hit("Chr2", 7000, 7081, score=132), site_score=112), "paralog")
        result = fp.genome_label(fp.summarize_hits([hit("Chr2", 7000, 7081, score=112)]), marker, 132)
        self.assertEqual((result.target_score, result.target_from, result.off_target_score), (132, "site", 112))
        unplaced = fp.place_marker("marker7", None, self.contigs, 81)
        self.assertEqual(fp.genome_label(fp.summarize_hits([hit("Chr1", 0, 81)]), unplaced).label, "unplaced")
        multicopy = fp.place_marker(
            "marker7", fp.GenomePlacement((hit("Chr1", 1000, 1081), hit("Chr2", 5000, 5081))), self.contigs, 81
        )
        self.assertEqual(fp.genome_label(fp.summarize_hits([hit("Chr1", 1000, 1081)]), multicopy).label, "ambiguous")

    def test_target_score_is_the_better_of_both_searches(self):
        marker = fp.MarkerPlacement((hit("Chr1", 1000, 1081),), "ref_tag", 1, "agree")
        off = hit("Chr2", 7000, 7081, score=122)
        weaker = fp.genome_label(fp.summarize_hits([hit("Chr1", 1010, 1081, score=112), off]), marker, 132)
        absent = fp.genome_label(fp.summarize_hits([off]), marker, 132)
        self.assertEqual(weaker, absent)  # a weaker genome-wide hit at the site changes nothing
        self.assertEqual((weaker.label, weaker.target_score, weaker.target_from), ("on_target", 132, "site"))
        stronger = fp.genome_label(fp.summarize_hits([hit("Chr1", 1000, 1081, score=142), off]), marker, 132)
        self.assertEqual((stronger.label, stronger.target_score, stronger.target_from), ("on_target", 142, "genome"))
        tied = fp.genome_label(fp.summarize_hits([hit("Chr1", 1000, 1081, score=132), off]), marker, 132)
        self.assertEqual((tied.target_score, tied.target_from), (132, "genome"))


class FilterTests(unittest.TestCase):
    """filter_paralogs on synthetic MADCs with preset genome-wide hits.

    One marker, tag at Chr1:1001-1081.  The genome-wide search finds the
    paralog's exact copy on Chr2 and an equal-scoring copy of the ambiguous
    haplotype; it misses the diverged allele at its own site, which the
    realignment to the marker's site then scores.
    """

    clone = "Chr1_001040"

    def setUp(self):
        rng = random.Random(7)

        def other(base):
            return "A" if base != "A" else "C"

        ref = "".join(rng.choice("ACGT") for _ in range(81))
        self.ref = ref
        self.alt = mutate(ref, 40, other(ref[40]))
        self.paralog = mutate(ref, 10, other(ref[10]))
        self.ambiguous = mutate(ref, 60, other(ref[60]))
        self.true_match = mutate(ref, 70, other(ref[70]))
        self.diverged = mutate(mutate(mutate(ref, 5, other(ref[5])), 30, other(ref[30])), 75, other(ref[75]))
        filler = "".join(rng.choice("ACGT") for _ in range(12000))
        genome = {"Chr1": filler[:1000] + ref + filler[1000:2000], "Chr2": filler[2000:12000]}
        target = hit("Chr1", 1000, 1081)
        self.aligner = FakeAligner(genome, {
            self.ref: [target],
            self.alt: [hit("Chr1", 1000, 1081, score=152)],
            self.paralog: [hit("Chr2", 7000, 7081), hit("Chr1", 1000, 1081, score=152)],
            self.ambiguous: [hit("Chr2", 9000, 9081, score=152), hit("Chr1", 1000, 1081, score=152)],
            self.true_match: [hit("Chr1", 1000, 1081, score=152)],
            self.diverged: [hit("Chr2", 5000, 5081, score=112)],
        })
        self.specs = [
            ("Ref", self.ref, [10, 0, 5]),
            ("Alt", self.alt, [0, 9, 5]),
            ("RefMatch", self.paralog, [4, 4, 4]),
            ("RefMatch", self.ambiguous, [1, 2, 0]),
            ("AltMatch", self.true_match, [0, 3, 0]),
            ("Other", self.diverged, [1, 1, 1]),
        ]
        self.rows = [madc_row(self.clone, *spec) for spec in self.specs]
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def run_filter(self, madc, **options):
        output = self.dir / "out_MADC.csv"
        with contextlib.redirect_stderr(io.StringIO()):
            summary = fp.filter_paralogs(
                madc, output, self.aligner, self.dir / "labels.csv", self.dir / "markers.csv", **options
            )
        return output, summary

    def test_default_drops_paralog_and_ambiguous_rows_and_tallies_reads(self):
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        output, summary = self.run_filter(madc)
        kept = read_rows(output)
        self.assertEqual(kept[0], read_rows(madc)[0])  # header untouched in drop mode
        self.assertEqual(kept[1:], [self.rows[0], self.rows[1], self.rows[4], self.rows[5]])  # copied unchanged
        rows = summary["rows"]
        self.assertEqual(rows["input"], {"rows": 6, "reads": 50})
        self.assertEqual(rows["output"], {"rows": 4, "reads": 35})
        self.assertEqual(rows["by_label"]["paralog"], {"rows": 1, "reads": 12})
        self.assertEqual(rows["by_label"]["ambiguous"], {"rows": 1, "reads": 3})
        self.assertEqual(rows["by_outcome"]["dropped"], {"rows": 2, "reads": 15})
        self.assertEqual(rows["by_category"]["RefMatch"]["paralog"]["rows"], 1)
        self.assertEqual(
            summary["haplotypes"]["realigned_at_site"], {"haplotypes": 6, "found_at_site": 6, "target_from_site": 1}
        )
        self.assertEqual(self.aligner.windows.count(("Chr1", 980, 1101)), 1)  # the marker's site, indexed once
        labels = {row["AlleleSequence"]: row for row in read_dicts(self.dir / "labels.csv")}
        self.assertEqual(len(labels), 6)  # dropped rows keep a record
        paralog = labels[self.paralog]
        self.assertEqual((paralog["Label"], paralog["Action"]), ("paralog", "dropped"))
        self.assertEqual((paralog["TargetScore"], paralog["OffTargetScore"]), ("152", "162"))
        self.assertEqual(paralog["OffTargetHit"], "Chr2:7001-7081(+)")  # 1-based
        self.assertEqual(labels[self.ambiguous]["Label"], "ambiguous")
        diverged = labels[self.diverged]
        self.assertEqual((diverged["Label"], diverged["TargetFrom"], diverged["TargetScore"]), ("on_target", "site", "132"))
        self.assertEqual((labels[self.ref]["TargetFrom"], labels[self.ref]["OffTargetScore"]), ("genome", "."))
        markers = read_dicts(self.dir / "markers.csv")
        self.assertEqual(markers[0]["Placement"], "ref_tag")
        self.assertEqual(markers[0]["CloneIDCheck"], "agree")
        self.assertEqual(markers[0]["Sites"], "Chr1:1001-1081(+)")

    def test_weaker_genome_wide_hit_at_the_site_does_not_hide_the_realignment(self):
        # The genome-wide search aligns the diverged allele only partly at its site (112), below its Chr2
        # hit (122); realigned at the site it scores 132, so it is kept.
        self.aligner._hits[self.diverged] = [hit("Chr2", 5000, 5081, score=122), hit("Chr1", 1010, 1081, score=112)]
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        self.run_filter(madc)
        diverged = {row["AlleleSequence"]: row for row in read_dicts(self.dir / "labels.csv")}[self.diverged]
        self.assertEqual(
            (diverged["Label"], diverged["TargetScore"], diverged["TargetFrom"], diverged["OffTargetScore"]),
            ("on_target", "132", "site", "122"),
        )

    def test_ref_tag_site_missing_from_the_genome_wide_hits(self):
        # The genome-wide search returns only a weaker Chr2 site for the Ref tag.  Aligned in the window
        # around its CloneID position, the tag finds Chr1:1001-1081, which becomes the marker's site.
        self.aligner._hits[self.ref] = [hit("Chr2", 6000, 6081, score=96)]
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        output, summary = self.run_filter(madc)
        marker = read_dicts(self.dir / "markers.csv")[0]
        self.assertEqual(
            (marker["Placement"], marker["RefTagSites"], marker["CloneIDCheck"], marker["Sites"]),
            ("clone_id_site", "1", "disagree", "Chr1:1001-1081(+)"),
        )
        self.assertEqual(summary["markers"]["clone_id_disagreements"], [self.clone])
        self.assertEqual(read_rows(output)[1:], [self.rows[0], self.rows[1], self.rows[4], self.rows[5]])
        # Ignored, the CloneID cannot move the marker: it stays on Chr2, and the real AltMatch allele is dropped.
        output, _ = self.run_filter(madc, clone_positions="ignore")
        marker = read_dicts(self.dir / "markers.csv")[0]
        self.assertEqual((marker["Placement"], marker["Sites"]), ("ref_tag", "Chr2:6001-6081(+)"))
        labels = {row["AlleleSequence"]: row["Label"] for row in read_dicts(self.dir / "labels.csv")}
        self.assertEqual(labels[self.true_match], "paralog")

    def test_ninth_equal_best_ref_tag_site(self):
        # Nine exact copies of the Ref tag, and the one the CloneID names sorts after the eight a haplotype keeps.
        self.aligner._hits[self.ref] = [hit("Chr1", start, start + 81) for start in range(0, 800, 100)] + [
            hit("Chr1", 1000, 1081)
        ]
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        output, _ = self.run_filter(madc)
        marker = read_dicts(self.dir / "markers.csv")[0]
        self.assertEqual(
            (marker["Placement"], marker["RefTagSites"], marker["CloneIDCheck"], marker["Sites"]),
            ("ref_tag+clone_id", "9", "agree", "Chr1:1001-1081(+)"),
        )
        self.assertEqual(read_rows(output)[1:], [self.rows[0], self.rows[1], self.rows[4], self.rows[5]])

    def test_canonical_rows_are_never_dropped(self):
        # The genome-wide search places this marker's Alt tag elsewhere, better than at its site.
        self.aligner._hits[self.alt] = [hit("Chr2", 3000, 3081)]
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        output, summary = self.run_filter(madc)
        self.assertIn(self.alt, [row[2] for row in read_rows(output)])
        self.assertEqual(summary["rows"]["by_outcome"]["protected"], {"rows": 1, "reads": 14})
        self.assertEqual(summary["rows"]["by_category"]["Alt"], {"paralog": {"rows": 1, "reads": 14}})

    def test_flag_mode_labels_every_row_in_the_picref_column(self):
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows)
        output, summary = self.run_filter(madc, paralog_action="flag", ambiguous_action="keep")
        rows = read_rows(output)
        column = rows[0].index("GenomeLabel")
        self.assertEqual(column, mp.MADC_METADATA_COLUMNS.index("PICRef"))
        self.assertEqual(len(rows[0]), 16 + len(SAMPLES))
        self.assertEqual(
            [row[column] for row in rows[1:]],
            ["on_target", "on_target", "paralog", "ambiguous", "on_target", "on_target"],
        )
        self.assertEqual(summary["rows"]["by_outcome"]["flagged"]["rows"], 1)
        self.assertEqual(summary["rows"]["by_outcome"]["kept"]["rows"], 5)
        self.assertEqual(summary["genome_label_column"], "GenomeLabel")
        # Filtering an already-flagged table rewrites the column in place.
        again = self.dir / "again.csv"
        with contextlib.redirect_stderr(io.StringIO()):
            fp.filter_paralogs(output, again, self.aligner, self.dir / "l2.csv", self.dir / "m2.csv", paralog_action="flag")
        self.assertEqual(read_rows(again), rows[:4] + rows[5:])  # ambiguous now dropped, nothing else changed

    def test_vendor_madc_keeps_preamble_and_statistics(self):
        preamble = [
            [""] * 16 + ["DLet22"] * 3,
            ["*"] * 16 + ["1", "1", "2"],
            ["*"] * 16 + ["3077912", "3077913", "3077914"],
        ]
        vendor_rows = [madc_row(self.clone, *spec, vendor=True) for spec in self.specs]
        madc = self.dir / "vendor.csv"
        write_madc(madc, vendor_rows, vendor=True, preamble=preamble)
        output, summary = self.run_filter(madc)
        rows = read_rows(output)
        self.assertEqual(rows[:3], preamble)
        self.assertEqual(rows[3][11:13], ["PICSnp", "AvgPIC"])
        self.assertEqual(rows[4:], [vendor_rows[0], vendor_rows[1], vendor_rows[4], vendor_rows[5]])
        self.assertEqual(summary["rows"]["output"]["reads"], 35)

    def test_provenance_is_appended_to_a_madcats_preamble_only(self):
        preamble = [["* MADCATS v0.13.0 tag_length=81"] + ["*"] * 15 + ["proj"] * 3]
        madc = self.dir / "in.csv"
        write_madc(madc, self.rows, preamble=preamble)
        output, _ = self.run_filter(madc, provenance="filter-paralogs paralogs=drop")
        self.assertEqual(read_rows(output)[0][0], "* MADCATS v0.13.0 tag_length=81 | filter-paralogs paralogs=drop")

    def test_flag_mode_needs_the_picref_column(self):
        madc = self.dir / "short.csv"
        with madc.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["AlleleID", "CloneID", "AlleleSequence"] + SAMPLES)
            writer.writerow(["Chr1_001040|Ref", "Chr1_001040", self.ref, "1", "2", "3"])
        with self.assertRaises(mp.MadcError):
            self.run_filter(madc, paralog_action="flag")
        output, summary = self.run_filter(madc)  # drop mode works on the short layout
        self.assertEqual(summary["rows"]["output"], {"rows": 1, "reads": 6})

    def test_clone_position_from_another_assembly(self):
        # The Ref tag does not align, and the CloneID names Chr2:500, where this genome has no copy of the
        # locus (as when CloneIDs carry another assembly's coordinates).  Used, that position becomes the
        # marker's site and the real AltMatch allele, aligning only at Chr1:1001, is dropped as a paralog.
        self.aligner._hits[self.ref] = []
        clone = "Chr2_000500"
        rows = [madc_row(clone, *spec) for spec in self.specs]
        madc = self.dir / "in.csv"
        write_madc(madc, rows)
        _, used = self.run_filter(madc)
        self.assertEqual(read_dicts(self.dir / "markers.csv")[0]["Placement"], "clone_id")
        labels = {row["AlleleSequence"]: row for row in read_dicts(self.dir / "labels.csv")}
        self.assertEqual((labels[self.true_match]["Label"], labels[self.true_match]["Action"]), ("paralog", "dropped"))
        output, ignored = self.run_filter(madc, clone_positions="ignore")
        self.assertEqual(ignored["markers"]["placement"], {"unplaced": 1})
        self.assertEqual(ignored["markers"]["clone_id_check"], {"ignored": 1})
        self.assertEqual(read_rows(output)[1:], rows)  # an unplaced marker keeps every row
        self.assertEqual({row["Label"] for row in read_dicts(self.dir / "labels.csv")}, {"unplaced"})

    def test_hapapp_fixed_allele_id_madc(self):
        # HapApp writes a leading Code_version column, then AlleleID, CloneID, AlleleSequence and the
        # samples, with fixed IDs such as |Ref_0001; counts of rows it merged come out as floats.
        self.aligner._hits[self.alt] = [hit("Chr2", 3000, 3081)]  # the Alt tag aligns better elsewhere
        fixed = {"Ref": "Ref_0001", "Alt": "Alt_0002"}
        rows = [
            ["v1", f"{self.clone}|{fixed.get(category, f'{category}_{serial:04d}')}", self.clone, sequence]
            + [f"{count}.0" for count in counts]
            for serial, (category, sequence, counts) in enumerate(self.specs, 1)
        ]
        madc = self.dir / "hapapp.csv"
        with madc.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(["Code_version", "AlleleID", "CloneID", "AlleleSequence"] + SAMPLES)
            writer.writerows(rows)
        output, summary = self.run_filter(madc)
        self.assertEqual(read_rows(output), [read_rows(madc)[0], rows[0], rows[1], rows[4], rows[5]])
        self.assertEqual(summary["rows"]["input"], {"rows": 6, "reads": 50})
        self.assertEqual(set(summary["rows"]["by_category"]), {"Ref", "Alt", "RefMatch", "AltMatch", "Other"})
        self.assertEqual(summary["rows"]["by_outcome"]["protected"], {"rows": 1, "reads": 14})
        with self.assertRaises(mp.MadcError):  # no PICRef column to carry a flag
            self.run_filter(madc, paralog_action="flag")
        parsed = mp.read_madc(madc)  # validate reads it too
        self.assertEqual(parsed["samples"], SAMPLES)
        self.assertEqual(parsed["records"][0]["counts"], [10, 0, 5])


class CountOverflowTests(unittest.TestCase):
    """Vendor MADCs that stored counts in a signed byte: 128-255 come out negative."""

    def test_signature_is_reported(self):
        check = mp.CountOverflowCheck()
        check.add("m1", [3, 127, -128, -2])
        check.add("m2", [0, 5])
        self.assertEqual(check.summary(), {
            "negative_cells": 2, "markers_with_negative_cells": 1, "minimum": -128, "maximum": 127,
            "signed_8bit_signature": True,
        })
        self.assertIn("signed 8-bit overflow", check.warning(Path("x.csv")))
        clean = mp.CountOverflowCheck()
        clean.add("m1", [0, 127, 300])
        self.assertIsNone(clean.warning(Path("x.csv")))
        self.assertFalse(clean.summary()["signed_8bit_signature"])

    def test_both_readers_warn(self):
        with tempfile.TemporaryDirectory() as temporary:
            madc = Path(temporary) / "vendor.csv"
            write_madc(madc, [madc_row("Chr1_001040", "Ref", "ACGT" * 20 + "A", [127, -100, 4], vendor=True)], vendor=True)
            with contextlib.redirect_stderr(io.StringIO()) as log:
                parsed = mp.read_madc(madc)
            self.assertIn("negative read counts", log.getvalue())
            self.assertEqual(parsed["count_overflow"]["negative_cells"], 1)
            aligner = FakeAligner({"Chr1": "A" * 2000}, {})
            with contextlib.redirect_stderr(io.StringIO()) as log:
                summary = fp.filter_paralogs(
                    madc, Path(temporary) / "out.csv", aligner, Path(temporary) / "l.csv", Path(temporary) / "m.csv"
                )
            self.assertIn("signed 8-bit overflow", log.getvalue())
            self.assertTrue(summary["count_overflow"]["signed_8bit_signature"])

    def test_a_failed_run_leaves_no_partial_files(self):
        tag = "ACGT" * 20 + "A"
        with tempfile.TemporaryDirectory() as temporary:
            madc = Path(temporary) / "in.csv"
            write_madc(madc, [
                madc_row("Chr1_001040", "Ref", tag, [1, 2, 3]),
                madc_row("Chr1_001040", "Alt", mutate(tag, 40, "C"), [1, "n/a", 3]),
            ])
            outputs = [Path(temporary) / name for name in ("out.csv", "labels.csv", "markers.csv")]
            with self.assertRaises(mp.MadcError):
                fp.filter_paralogs(madc, outputs[0], FakeAligner({"Chr1": "A" * 2000}, {}), outputs[1], outputs[2])
            self.assertEqual([path.exists() for path in outputs], [False, False, False])

    def test_per_row_read_total_is_not_a_sample(self):
        # Processed MADCs carry a readCountSum column between the sequence and the samples.
        tag = "ACGT" * 20 + "A"
        with tempfile.TemporaryDirectory() as temporary:
            madc = Path(temporary) / "processed.csv"
            with madc.open("w", newline="") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(["AlleleID", "CloneID", "AlleleSequence", "readCountSum"] + SAMPLES)
                writer.writerow(["Chr1_001040|Ref", "Chr1_001040", tag, "50", "10", "15", "25"])
            parsed = mp.read_madc(madc)
            self.assertEqual(parsed["samples"], SAMPLES)
            self.assertEqual(parsed["records"][0]["counts"], [10, 15, 25])
            aligner = FakeAligner({"Chr1": "A" * 2000}, {})
            summary = fp.filter_paralogs(
                madc, Path(temporary) / "out.csv", aligner, Path(temporary) / "l.csv", Path(temporary) / "m.csv"
            )
            self.assertEqual(summary["rows"]["input"], {"rows": 1, "reads": 50})
            self.assertEqual(read_rows(Path(temporary) / "out.csv"), read_rows(madc))


class CommandTests(unittest.TestCase):
    def test_missing_mappy_names_the_extra(self):
        with mock.patch.dict(sys.modules, {"mappy": None}):
            with self.assertRaises(mp.MadcError) as caught:
                fp.GenomeAligner.open(Path("genome.fa"), None, 1)
        self.assertIn("madcats[genome]", str(caught.exception))

    def test_needs_a_genome_or_an_index(self):
        with tempfile.TemporaryDirectory() as temporary:
            madc = Path(temporary) / "in.csv"
            madc.write_text("AlleleID,CloneID,AlleleSequence\n")
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as caught:
                mp.main(["filter-paralogs", "--madc", str(madc), "--output", str(Path(temporary) / "out.csv")])
        self.assertEqual(caught.exception.code, 2)

    def test_script_run_shares_madc_error_with_the_filter_module(self):
        # python3 madcats.py runs the file as __main__; an error raised in madcats_paralogs must still
        # be reported as "error: ..." with exit code 2, not a traceback.
        with tempfile.TemporaryDirectory() as temporary:
            madc = Path(temporary) / "in.csv"
            madc.write_text("AlleleID,CloneID,AlleleSequence\n")
            result = subprocess.run(
                [sys.executable, str(REPOSITORY / "madcats.py"), "filter-paralogs", "--madc", str(madc),
                 "--output", str(Path(temporary) / "out.csv")],
                capture_output=True, text=True,
            )
        self.assertEqual(result.returncode, 2, result.stderr)
        self.assertIn("error: filter-paralogs needs --genome", result.stderr)

    def test_other_commands_do_not_import_mappy(self):
        with mock.patch.dict(sys.modules, {"mappy": None}):
            with contextlib.redirect_stdout(io.StringIO()):
                args = mp.build_parser().parse_args(["filter-paralogs", "--madc", "x", "--output", "y"])
        self.assertEqual((args.paralogs, args.ambiguous, args.clone_positions), ("drop", "drop", "use"))


@unittest.skipUnless(HAS_MAPPY, "mappy not installed (pip install 'madcats[genome]')")
class GenomeAlignerTests(unittest.TestCase):
    """The real command on a synthetic 3-contig genome.

    Marker A sits on chr1 and has a 2-mismatch paralogous copy on chr2.  Marker B
    sits on chr1 and has a copy on chr3 differing at two positions, so a
    haplotype carrying one base of each copy ties.  Marker C's tag is
    duplicated exactly on chr1 and chr3; its CloneID names the chr1 copy.
    """

    @classmethod
    def setUpClass(cls):
        rng = random.Random(11)
        contigs = {name: [rng.choice("ACGT") for _ in range(6000)] for name in ("chr1", "chr2", "chr3")}

        def segment(name, start):
            return "".join(contigs[name][start : start + 81])

        def other(base):
            return "A" if base != "A" else "C"

        cls.a_ref = segment("chr1", 1000)
        copy = mutate(mutate(cls.a_ref, 20, other(cls.a_ref[20])), 50, other(cls.a_ref[50]))
        contigs["chr2"][3000:3081] = copy
        cls.a_paralog = copy  # a perfect hit on chr2, two mismatches at the target
        cls.a_alt = mutate(cls.a_ref, 40, other(cls.a_ref[40]))
        cls.a_match = mutate(cls.a_ref, 65, other(cls.a_ref[65]))  # one mismatch, only at the target

        cls.b_ref = segment("chr1", 4000)
        b_copy = mutate(mutate(cls.b_ref, 30, other(cls.b_ref[30])), 60, other(cls.b_ref[60]))
        contigs["chr3"][2000:2081] = b_copy
        cls.b_alt = mutate(cls.b_ref, 45, other(cls.b_ref[45]))
        cls.b_tie = mutate(cls.b_ref, 30, other(cls.b_ref[30]))  # one mismatch against each copy

        cls.c_ref = segment("chr1", 5500)
        contigs["chr3"][4500:4581] = cls.c_ref
        cls.c_alt = mutate(cls.c_ref, 40, other(cls.c_ref[40]))

        cls.temporary = tempfile.TemporaryDirectory()
        cls.dir = Path(cls.temporary.name)
        cls.genome = cls.dir / "genome.fa"
        with cls.genome.open("w") as handle:
            for name, bases in contigs.items():
                handle.write(f">{name}\n{''.join(bases)}\n")
        rows = [
            madc_row("chr1_001041", "Ref", cls.a_ref, [10, 0, 5]),
            madc_row("chr1_001041", "Alt", cls.a_alt, [0, 10, 5]),
            madc_row("chr1_001041", "RefMatch", cls.a_paralog, [6, 6, 6]),
            madc_row("chr1_001041", "RefMatch", cls.a_match, [2, 0, 0]),
            madc_row("chr1_004046", "Ref", cls.b_ref, [8, 8, 8]),
            madc_row("chr1_004046", "Alt", cls.b_alt, [0, 4, 0]),
            madc_row("chr1_004046", "RefMatch", cls.b_tie, [3, 3, 3]),
            madc_row("chr1_005541", "Ref", cls.c_ref, [9, 9, 9]),
            madc_row("chr1_005541", "Alt", cls.c_alt, [1, 1, 1]),
        ]
        cls.madc = cls.dir / "in_MADC.csv"
        write_madc(cls.madc, rows)
        cls.index = cls.dir / "genome.mmi"
        cls.output = cls.dir / "out_MADC.csv"
        with contextlib.redirect_stderr(io.StringIO()):
            fp.run_filter_paralogs(mp.build_parser().parse_args([
                "filter-paralogs", "--madc", str(cls.madc), "--genome", str(cls.genome),
                "--genome-index", str(cls.index), "--output", str(cls.output), "--jobs", "2",
            ]))
        cls.qc = json.loads(cls.output.with_suffix(".csv.qc.json").read_text())
        cls.labels = {row["AlleleSequence"]: row for row in read_dicts(cls.output.with_suffix(".csv.genome_labels.csv"))}
        cls.markers = {row["CloneID"]: row for row in read_dicts(cls.output.with_suffix(".csv.genome_markers.csv"))}

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_labels_from_the_genome(self):
        expected = {
            self.a_ref: "on_target", self.a_alt: "on_target", self.a_paralog: "paralog", self.a_match: "on_target",
            self.b_ref: "on_target", self.b_alt: "on_target", self.b_tie: "ambiguous",
            self.c_ref: "ambiguous", self.c_alt: "ambiguous",
        }
        self.assertEqual({sequence: row["Label"] for sequence, row in self.labels.items()}, expected)
        self.assertTrue(self.labels[self.a_paralog]["OffTargetHit"].startswith("chr2:3001-3081"))
        self.assertEqual(self.labels[self.a_paralog]["OffTargetScore"], "162")

    def test_marker_placement(self):
        self.assertEqual(self.markers["chr1_001041"]["Sites"], "chr1:1001-1081(+)")
        self.assertEqual(self.markers["chr1_001041"]["CloneIDCheck"], "agree")
        multicopy = self.markers["chr1_005541"]
        self.assertEqual((multicopy["Placement"], multicopy["RefTagSites"]), ("ref_tag+clone_id", "2"))
        self.assertEqual(self.qc["markers"]["multicopy_ref_tags"], ["chr1_005541"])

    def test_output_and_qc(self):
        kept = [row[2] for row in read_rows(self.output)[1:]]
        self.assertNotIn(self.a_paralog, kept)
        self.assertNotIn(self.b_tie, kept)
        self.assertIn(self.c_ref, kept)  # canonical, protected despite its multi-copy tag
        self.assertEqual(self.qc["rows"]["by_outcome"]["protected"]["rows"], 2)
        self.assertEqual(self.qc["rows"]["by_outcome"]["dropped"], {"rows": 2, "reads": 27})
        self.assertEqual(self.qc["aligner"]["k"], 21)

    def test_saved_index_is_reused(self):
        self.assertTrue(self.index.is_file())
        output = self.dir / "from_index_MADC.csv"
        with contextlib.redirect_stderr(io.StringIO()) as log:
            fp.run_filter_paralogs(mp.build_parser().parse_args(
                ["filter-paralogs", "--madc", str(self.madc), "--genome-index", str(self.index), "--output", str(output)]
            ))
        self.assertIn("Loading index", log.getvalue())
        self.assertEqual(read_rows(output), read_rows(self.output))

    def test_clone_positions_ignored(self):
        # Marker C's Ref tag has equal-best copies on chr1 and chr3; only its CloneID picked chr1.
        output = self.dir / "ignored_MADC.csv"
        with contextlib.redirect_stderr(io.StringIO()):
            fp.run_filter_paralogs(mp.build_parser().parse_args([
                "filter-paralogs", "--madc", str(self.madc), "--genome-index", str(self.index),
                "--output", str(output), "--clone-positions", "ignore",
            ]))
        markers = {row["CloneID"]: row for row in read_dicts(output.with_suffix(".csv.genome_markers.csv"))}
        self.assertEqual((markers["chr1_005541"]["Placement"], markers["chr1_005541"]["Sites"].count(";")), ("ref_tag_multicopy", 1))
        self.assertEqual({row["CloneIDCheck"] for row in markers.values()}, {"ignored"})
        labels = {row["AlleleSequence"]: row["Label"] for row in read_dicts(output.with_suffix(".csv.genome_labels.csv"))}
        self.assertEqual(labels, {sequence: row["Label"] for sequence, row in self.labels.items()})
        qc = json.loads(output.with_suffix(".csv.qc.json").read_text())
        self.assertEqual(qc["parameters"]["clone_positions"], "ignore")


@unittest.skipUnless(HAS_MAPPY, "mappy not installed (pip install 'madcats[genome]')")
class ReviewCaseTests(unittest.TestCase):
    """Counterexamples from an external review (docs/reviews), each on a tiny genome with the real aligner."""

    def setUp(self):
        self.rng = random.Random(34)
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def dna(self, length):
        return "".join(self.rng.choice("ACGT") for _ in range(length))

    def aligner(self, contigs, name="genome"):
        fasta = self.dir / f"{name}.fa"
        fasta.write_text("".join(f">{contig}\n{sequence}\n" for contig, sequence in contigs.items()))
        with contextlib.redirect_stderr(io.StringIO()):
            return fp.GenomeAligner.open(fasta, None, 1)

    def run_filter(self, aligner, rows):
        madc = self.dir / "in.csv"
        write_madc(madc, rows)
        with contextlib.redirect_stderr(io.StringIO()):
            fp.filter_paralogs(madc, self.dir / "out.csv", aligner, self.dir / "labels.csv", self.dir / "markers.csv")
        labels = {row["AlleleSequence"]: row for row in read_dicts(self.dir / "labels.csv")}
        return labels, read_dicts(self.dir / "markers.csv")

    def test_assembly_carrying_the_alt_allele_at_the_target(self):
        # chr1 is the design site but carries the Alt allele; a paralog on chr2 carries Ref.  The Ref tag
        # aligns better on chr2, yet the marker is placed where its CloneID says.
        ref = self.dna(81)
        alt = changed(ref, 40)
        allele = changed(alt, 65)  # one SNP from the target, two from the paralog
        aligner = self.aligner({"chr1": self.dna(500) + alt + self.dna(500), "chr2": self.dna(500) + ref + self.dna(500)})
        specs = [("Ref", ref, [5, 5, 5]), ("Alt", alt, [5, 5, 5]), ("AltMatch", allele, [3, 3, 3])]
        labels, markers = self.run_filter(aligner, [madc_row("chr1_000540", *spec) for spec in specs])
        self.assertEqual(
            (markers[0]["Placement"], markers[0]["CloneIDCheck"], markers[0]["Sites"]),
            ("clone_id_site", "disagree", "chr1:501-581(+)"),
        )
        found = labels[allele]
        self.assertEqual(
            (found["Label"], found["Action"], found["TargetScore"], found["OffTargetScore"]), ("on_target", "kept", "152", "142")
        )
        # A CloneID position where the Ref tag does not align leaves the marker unplaced, keeping every row.
        labels, markers = self.run_filter(aligner, [madc_row("chr1_000900", *spec) for spec in specs])
        self.assertEqual((markers[0]["Placement"], markers[0]["CloneIDCheck"]), ("unplaced", "disagree"))
        self.assertEqual({(row["Label"], row["Action"]) for row in labels.values()}, {("unplaced", "kept")})

    def test_markers_sharing_a_site_do_not_crowd_each_other_out(self):
        # The allele has no 21-base seed at its own site, and its genome-wide hit is a 50-base partial copy.
        ref = self.dna(81)
        allele = changed(ref, 10, 30, 50, 70)
        aligner = self.aligner({"target": self.dna(500) + ref + self.dna(500), "off": self.dna(500) + allele[:50] + self.dna(500)})
        for copies in (1, 30):  # more markers on one site than the short-read preset returns hits
            rows = [
                row for number in range(copies)
                for row in (madc_row(f"m{number}", "Ref", ref, [9, 9, 9]), madc_row(f"m{number}", "RefMatch", allele, [3, 3, 3]))
            ]
            write_madc(self.dir / "in.csv", rows)
            with contextlib.redirect_stderr(io.StringIO()):
                fp.filter_paralogs(
                    self.dir / "in.csv", self.dir / "out.csv", aligner, self.dir / "labels.csv", self.dir / "markers.csv"
                )
            found = Counter(
                (row["Label"], row["TargetScore"], row["OffTargetScore"])
                for row in read_dicts(self.dir / "labels.csv") if row["AlleleID"].endswith("|RefMatch")
            )
            self.assertEqual(found, {("on_target", "122", "100"): copies})

    def test_adjacent_tandem_copies(self):
        repeat = self.dna(80)
        ref = repeat + repeat[0]
        allele = changed(ref, 40)  # one SNP from both copies
        aligner = self.aligner({"chr1": self.dna(500) + repeat + repeat + repeat[0] + self.dna(500)})
        labels, markers = self.run_filter(
            aligner, [madc_row("chr1_000540", "Ref", ref, [5, 5, 5]), madc_row("chr1_000540", "RefMatch", allele, [3, 3, 3])]
        )
        self.assertEqual(
            (markers[0]["Placement"], markers[0]["RefTagSites"], markers[0]["Sites"]), ("ref_tag+clone_id", "2", "chr1:501-581(+)")
        )
        found = labels[allele]
        self.assertEqual(
            (found["Label"], found["TargetScore"], found["OffTargetScore"], found["OffTargetHit"]),
            ("ambiguous", "152", "152", "chr1:581-661(+)"),
        )

    def test_multipart_index_is_refused(self):
        import mappy

        parts = []
        for number in (1, 2):
            fasta = self.dir / f"part{number}.fa"
            fasta.write_text(f">chr{number}\n{self.dna(1000)}\n")
            parts.append(self.dir / f"part{number}.mmi")
            mappy.Aligner(str(fasta), preset="sr", fn_idx_out=str(parts[-1]))
        joined = self.dir / "joined.mmi"  # what minimap2 -I writes for a genome larger than its batch
        joined.write_bytes(parts[0].read_bytes() + parts[1].read_bytes())
        self.assertEqual(fp._index_parts(parts[0]), [["chr1"]])
        self.assertEqual(fp._index_parts(joined), [["chr1"], ["chr2"]])
        self.assertIsNone(fp._index_parts(self.dir / "part1.fa"))
        truncated = self.dir / "truncated.mmi"
        truncated.write_bytes(parts[0].read_bytes()[:-5])
        with self.assertRaises(mp.MadcError):
            fp._index_parts(truncated)
        madc = self.dir / "in.csv"
        write_madc(madc, [madc_row("chr1_000100", "Ref", self.dna(81), [1, 1, 1])])
        output = self.dir / "out.csv"
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(mp.MadcError) as caught:
            fp.run_filter_paralogs(mp.build_parser().parse_args(
                ["filter-paralogs", "--madc", str(madc), "--genome-index", str(joined), "--output", str(output)]
            ))
        self.assertIn("in 2 parts", str(caught.exception))
        self.assertFalse(output.exists())
        # The same file is refused when --genome names it.
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(mp.MadcError) as caught:
            fp.run_filter_paralogs(mp.build_parser().parse_args(
                ["filter-paralogs", "--madc", str(madc), "--genome", str(joined), "--output", str(output)]
            ))
        self.assertIn("in 2 parts", str(caught.exception))
        self.assertFalse(output.exists())
        # A one-part index given as --genome is loaded, and no copy is saved over --genome-index.
        other = self.dir / "other.mmi"
        with contextlib.redirect_stderr(io.StringIO()) as log:
            fp.run_filter_paralogs(mp.build_parser().parse_args(
                ["filter-paralogs", "--madc", str(madc), "--genome", str(parts[0]), "--genome-index", str(other),
                 "--output", str(output)]
            ))
        self.assertIn("Loading index", log.getvalue())
        self.assertFalse(other.exists())
        self.assertTrue(output.exists())

    def test_label_change_between_assemblies_without_a_copy_change(self):
        # Both assemblies have the target and one other copy; only the other copy's sequence differs.  The
        # label changes, so a label change alone does not show that a copy was gained or lost.
        ref = self.dna(81)
        allele = changed(ref, 40)
        flanks = [self.dna(500) for _ in range(4)]
        labels = []
        for name, copy in (("a", changed(ref, 20, 50)), ("b", allele)):
            aligner = self.aligner({"target": flanks[0] + ref + flanks[1], "copy": flanks[2] + copy + flanks[3]}, name)
            ref_placement = fp.summarize_hits(aligner.hits(ref, aligner.new_buffer()), keep=None)
            marker = fp.place_marker("marker", ref_placement, aligner.contigs, 81, use_clone_position=False)
            labels.append(fp.genome_label(fp.summarize_hits(aligner.hits(allele, aligner.new_buffer())), marker).label)
        self.assertEqual(labels, ["on_target", "paralog"])


if __name__ == "__main__":
    unittest.main()
