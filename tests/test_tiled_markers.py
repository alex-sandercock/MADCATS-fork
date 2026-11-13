"""Tiled markers through ``generate``, ``filter-paralogs`` and ``call``.

Two markers that read the same tag window (tiles) share every read of their
amplicon.  ``generate`` credits each read to both.  The other two commands
were written without tiles in mind, so this checks what they do with such a
table: the paralog filter must not take one tile's rows for a copy of the
other, and ``call`` must get each tile's own variant right.

The panel is written both ways a tile's tag can carry the other tile's
variant: as an ambiguity code, the way DArTag panels write it, and as the
reference base.
"""
import contextlib
import csv
import importlib.util
import io
import json
import random
import tempfile
import unittest
from pathlib import Path

import madcats as mp

HAS_MAPPY = importlib.util.find_spec("mappy") is not None

TAG_LENGTH = 81
SNP_A, SNP_B = 30, 55                      # 0-based offsets of the two tiles' variants in the tag
AMBIGUITY_CODE = {
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("AC"): "M",
    frozenset("GT"): "K", frozenset("AT"): "W", frozenset("CG"): "S",
}
# Each sample's two copies of the amplicon, as (base at tile A's variant, base at tile B's):
# r is the reference base, a the alternate.
SAMPLES = {
    "s01": ("rr", "rr"), "s02": ("rr", "ar"), "s03": ("ar", "ar"), "s04": ("rr", "ra"), "s05": ("ra", "ra"),
    "s06": ("ar", "ra"), "s07": ("rr", "aa"), "s08": ("aa", "aa"), "s09": ("ar", "aa"), "s10": ("ra", "aa"),
}
READS_PER_COPY = 30


def with_base(sequence, offset, base):
    return sequence[:offset] + base + sequence[offset + 1:]


def another_base(base):
    return {"A": "C", "C": "T", "G": "A", "T": "G"}[base]


def run_quietly(arguments):
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return mp.main([str(argument) for argument in arguments])


def madc_rows(path):
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    header = next(index for index, row in enumerate(rows) if row and row[0] == "AlleleID")
    first_sample = rows[header].index(next(iter(SAMPLES)))
    return [(row[0], row[2], [int(value) for value in row[first_sample:]]) for row in rows[header + 1:]]


def alternate_dosages(vcf_path):
    """Records of a VCF as {(marker, position): (is the target, alternate dosage per sample)}."""
    records = {}
    with open(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            dosages = [column.split(":")[0].replace("|", "/").split("/").count("1") for column in fields[9:]]
            records[(fields[0], int(fields[1]))] = ("TARGET" in fields[7].split(";"), dosages)
    return records


class TiledPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.directory = tempfile.TemporaryDirectory()
        cls.root = Path(cls.directory.name)
        generator = random.Random(7)
        genome = "".join(generator.choice("ACGT") for _ in range(20000))
        tag = genome[5000 : 5000 + TAG_LENGTH]
        other_marker = genome[12000 : 12000 + TAG_LENGTH]
        alternate_a, alternate_b = another_base(tag[SNP_A]), another_base(tag[SNP_B])
        (cls.root / "genome.fasta").write_text(">chr1\n" + genome + "\n")
        cls.site = ("chr1", 5001, 5000 + TAG_LENGTH)

        copies = {
            "rr": tag,
            "ar": with_base(tag, SNP_A, alternate_a),
            "ra": with_base(tag, SNP_B, alternate_b),
            "aa": with_base(with_base(tag, SNP_A, alternate_a), SNP_B, alternate_b),
        }
        other_alternate = with_base(other_marker, 40, another_base(other_marker[40]))
        fastq_dir = cls.root / "fastq"
        fastq_dir.mkdir()
        for sample, pair in SAMPLES.items():
            reads = [copies[copy] for copy in pair for _ in range(READS_PER_COPY)]
            reads += [other_marker] * 20 + [other_alternate] * 20
            with (fastq_dir / f"{sample}.fastq").open("w") as handle:
                for number, read in enumerate(reads):
                    handle.write(f"@{sample}_{number}\n{read}\n+\n{'I' * len(read)}\n")
        cls.truth = {
            "tileA": [sum(copy[0] == "a" for copy in pair) for pair in SAMPLES.values()],
            "tileB": [sum(copy[1] == "a" for copy in pair) for pair in SAMPLES.values()],
        }

        cls.madc = {}
        for style in ("ambiguity codes", "reference bases"):
            # Each tile's tags hold the other tile's variant as a code for both bases, or as the reference base.
            at_b = AMBIGUITY_CODE[frozenset((tag[SNP_B], alternate_b))] if style == "ambiguity codes" else tag[SNP_B]
            at_a = AMBIGUITY_CODE[frozenset((tag[SNP_A], alternate_a))] if style == "ambiguity codes" else tag[SNP_A]
            ref_a, ref_b = with_base(tag, SNP_B, at_b), with_base(tag, SNP_A, at_a)
            loci = [
                mp.Locus(0, "tileA", tag[:20], "", tag[SNP_A], alternate_a, SNP_A, ref_a,
                         with_base(ref_a, SNP_A, alternate_a)),
                mp.Locus(1, "tileB", tag[:20], "", tag[SNP_B], alternate_b, SNP_B, ref_b,
                         with_base(ref_b, SNP_B, alternate_b)),
                mp.Locus(2, "other", other_marker[:20], "", other_marker[40], other_alternate[40], 40, other_marker,
                         other_alternate),
            ]
            prefix = cls.root / style.replace(" ", "_")
            mp.write_reference_tags(Path(str(prefix) + ".tags.csv"), loci)
            cls.madc[style] = Path(str(prefix) + ".madc.csv")
            exit_code = run_quietly([
                "generate", "--fastq-dir", fastq_dir, "--reference-tags", str(prefix) + ".tags.csv",
                "--tag-length", TAG_LENGTH, "--jobs", 1, "--reproducible-min-samples", 0, "--depth-min-reads", 1,
                "--cache-dir", str(prefix) + ".cache", "--output", cls.madc[style],
            ])
            assert exit_code == 0

    @classmethod
    def tearDownClass(cls):
        cls.directory.cleanup()

    def call(self, style, *options):
        prefix = self.root / (style.replace(" ", "_") + ".calls")
        self.assertEqual(run_quietly([
            "call", "--madc", self.madc[style], "--ploidy", 2, "--jobs", 1, "--min-dp", 0, "--min-gq", 0,
            "--output-prefix", prefix, *options,
        ]), 0)
        return prefix

    def test_generate_gives_both_tiles_every_read_of_the_amplicon(self):
        for style, madc in self.madc.items():
            with self.subTest(tags=style):
                totals = {}
                for allele_id, _, counts in madc_rows(madc):
                    marker = allele_id.split("|")[0]
                    totals[marker] = [total + count for total, count in zip(totals.get(marker, [0] * 10), counts)]
                self.assertEqual(totals["tileA"], [2 * READS_PER_COPY] * len(SAMPLES))
                self.assertEqual(totals["tileB"], totals["tileA"])
                with open(str(madc) + ".qc.json") as handle:
                    self.assertEqual(json.load(handle)["tile_families"], [["tileA", "tileB"]])

    def test_tags_with_ambiguity_codes_leave_the_ref_and_alt_rows_empty(self):
        # No read can equal a tag that holds a code, so every read is on a RefMatch or
        # AltMatch row.  A genotype taken from the Ref and Alt rows alone would have no reads.
        for allele_id, _, counts in madc_rows(self.madc["ambiguity codes"]):
            if allele_id.startswith("tile") and allele_id.split("|")[1] in ("Ref", "Alt"):
                self.assertEqual(sum(counts), 0, allele_id)
        reads_on_match_rows = sum(
            sum(counts) for allele_id, _, counts in madc_rows(self.madc["ambiguity codes"]) if "Match" in allele_id
        )
        self.assertEqual(reads_on_match_rows, 2 * 2 * READS_PER_COPY * len(SAMPLES))      # both tiles, two copies each

    @unittest.skipUnless(HAS_MAPPY, "mappy not installed (pip install 'madcats[genome]')")
    def test_paralog_filter_places_both_tiles_at_one_site_and_drops_nothing(self):
        for style, madc in self.madc.items():
            with self.subTest(tags=style):
                output = Path(str(madc) + ".filtered.csv")
                self.assertEqual(run_quietly([
                    "filter-paralogs", "--madc", madc, "--genome", self.root / "genome.fasta", "--output", output,
                    "--jobs", 1,
                ]), 0)
                with open(str(output) + ".genome_markers.csv", newline="") as handle:
                    sites = {row["CloneID"]: row["Sites"] for row in csv.DictReader(handle)}
                contig, first, last = self.site
                self.assertEqual(sites["tileA"], f"{contig}:{first}-{last}(+)")
                self.assertEqual(sites["tileB"], sites["tileA"])
                with open(str(output) + ".genome_labels.csv", newline="") as handle:
                    labels = list(csv.DictReader(handle))
                self.assertEqual({row["Label"] for row in labels}, {"on_target"})
                self.assertEqual({row["Action"] for row in labels}, {"kept"})
                self.assertEqual(madc_rows(output), madc_rows(madc))

    def test_call_gets_each_tile_its_own_variant(self):
        for style in self.madc:
            with self.subTest(tags=style):
                targets = alternate_dosages(str(self.call(style)) + ".targets.vcf")
                self.assertEqual(targets[("tileA", SNP_A + 1)], (True, self.truth["tileA"]))
                self.assertEqual(targets[("tileB", SNP_B + 1)], (True, self.truth["tileB"]))

    def test_the_other_tiles_variant_is_an_off_target_record_only_where_the_tag_has_a_base(self):
        # With the reference base in the tag, each tile also reports the other tile's variant,
        # from the same reads and so with the same genotypes.
        records = alternate_dosages(str(self.call("reference bases", "--vcf", "all")) + ".all.vcf")
        self.assertEqual(records[("tileA", SNP_B + 1)], (False, self.truth["tileB"]))
        self.assertEqual(records[("tileB", SNP_A + 1)], (False, self.truth["tileA"]))
        self.assertEqual(len(records), 5)
        # A position where the Ref tag holds an ambiguity code gets no off-target record:
        # the variant is reported once, as the other tile's target.
        records = alternate_dosages(str(self.call("ambiguity codes", "--vcf", "all")) + ".all.vcf")
        self.assertEqual(sorted(records), [("other", 41), ("tileA", SNP_A + 1), ("tileB", SNP_B + 1)])


if __name__ == "__main__":
    unittest.main()
