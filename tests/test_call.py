"""Tests of ``madcats call`` (madcats_call.py): the genotype model, then the command.

The parity tests compare against values computed by polyRAD and stored in
tests/fixtures/polyrad_parity.json (made by scripts/make_polyrad_fixture.R),
so they run without R.
"""

import contextlib
import csv
import io
import itertools
import json
import math
import random
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import madcats  # noqa: E402
import madcats_call as call  # noqa: E402

EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "madcats-demo"
FIXTURE_PATH = Path(__file__).resolve().parent / "fixtures" / "polyrad_parity.json"


def load_fixture():
    with open(FIXTURE_PATH) as handle:
        fixture = json.load(handle)
    fixture["sample_index"] = {name: index for index, name in enumerate(fixture["samples"])}
    fixture["allele_index"] = {name: index for index, name in enumerate(fixture["alleles"])}
    # Columns of the depth table that belong to each locus, in fixture order.
    fixture["locus_columns"] = {}
    for column, locus in enumerate(fixture["allele_locus"]):
        fixture["locus_columns"].setdefault(locus, []).append(column)
    return fixture


class LikelihoodTests(unittest.TestCase):
    def test_likelihoods_of_all_possible_counts_sum_to_one(self):
        settings = call.ModelSettings()
        marker_depth = 25
        for dosage in range(5):
            total = sum(
                math.exp(call.allele_log_likelihoods(allele_depth, marker_depth, 4, 0.3, settings)[dosage])
                for allele_depth in range(marker_depth + 1)
            )
            self.assertAlmostEqual(total, 1.0, places=10)

    def test_no_reads_gives_equal_likelihoods(self):
        self.assertEqual(call.allele_log_likelihoods(0, 0, 4, 0.3, call.ModelSettings()), [0.0] * 5)

    def test_very_deep_coverage_stays_finite(self):
        log_likelihoods = call.allele_log_likelihoods(2000000, 2000000, 6, 0.5, call.ModelSettings())
        self.assertTrue(all(math.isfinite(value) for value in log_likelihoods))
        self.assertEqual(log_likelihoods, sorted(log_likelihoods))
        self.assertLess(log_likelihoods[0], -100)

    def test_without_contamination_impossible_dosages_have_zero_likelihood(self):
        settings = call.ModelSettings(contamination_rate=0.0)
        log_likelihoods = call.allele_log_likelihoods(3, 10, 2, 0.5, settings)
        self.assertEqual(log_likelihoods[0], call.NEGATIVE_INFINITY)
        self.assertEqual(log_likelihoods[2], call.NEGATIVE_INFINITY)
        self.assertTrue(math.isfinite(log_likelihoods[1]))


class PriorAndPosteriorTests(unittest.TestCase):
    def test_hwe_prior_is_binomial(self):
        self.assertEqual(call.hwe_prior(2, 0.25), [0.5625, 0.375, 0.0625])
        self.assertAlmostEqual(sum(call.hwe_prior(6, 0.37)), 1.0)

    def test_posterior_falls_back_to_prior_when_nothing_is_possible(self):
        prior = [1.0, 0.0, 0.0]
        self.assertEqual(call.dosage_posterior([call.NEGATIVE_INFINITY, -1.0, -2.0], prior), prior)

    def test_posterior_survives_extreme_log_likelihoods(self):
        posterior = call.dosage_posterior([-5000.0, -5002.0, -9000.0], [1 / 3] * 3)
        self.assertAlmostEqual(sum(posterior), 1.0)
        self.assertAlmostEqual(posterior[0] / posterior[1], math.exp(2.0))


class MarkerTests(unittest.TestCase):
    def test_marker_without_reads_is_not_callable(self):
        result = call.call_marker([[0, 0], [0, 0]], [2, 2])
        self.assertFalse(result.callable)

    def test_two_allele_marginals_mirror_each_other(self):
        depths = [[30, 0], [14, 17], [0, 25], [22, 9], [5, 31]]
        ploidies = [2, 4, 4, 6, 3]
        result = call.call_marker(depths, ploidies)
        for sample_posteriors in result.posteriors:
            first, second = sample_posteriors
            for value, mirrored in zip(first, reversed(second)):
                self.assertAlmostEqual(value, mirrored, places=12)

    def test_clear_diploid_genotypes_are_called(self):
        depths = [[40, 0], [21, 19], [0, 38], [18, 22], [35, 0], [0, 41]]
        result = call.call_marker(depths, [2] * 6)
        self.assertTrue(result.converged)
        called = [call.most_probable_dosage(sample[1]) for sample in result.posteriors]
        self.assertEqual(called, [0, 1, 2, 1, 0, 2])

    def test_iteration_limit_is_reported(self):
        depths = [[9, 3], [2, 8], [5, 5], [1, 0]]
        result = call.call_marker(depths, [4] * 4, call.ModelSettings(tolerance=1e-15, max_iterations=2))
        self.assertEqual(result.iterations, 2)
        self.assertFalse(result.converged)
        self.assertGreater(result.final_change, 1e-15)

    def test_uniform_prior_is_a_single_pass(self):
        depths = [[9, 3], [2, 8], [0, 0]]
        result = call.call_marker(depths, [4] * 3, call.ModelSettings(prior=call.PRIOR_UNIFORM))
        self.assertEqual(result.iterations, 1)
        # A sample without reads keeps the prior, so no dosage stands out.
        self.assertEqual(result.posteriors[2][0], [0.2] * 5)
        self.assertIsNone(call.most_probable_dosage(result.posteriors[2][0]))


class PloidyCorrectionTests(unittest.TestCase):
    def test_dosages_that_already_fit_are_unchanged(self):
        posteriors = [[0.1, 0.8, 0.1], [0.1, 0.8, 0.1]]
        self.assertEqual(call.correct_dosages_to_ploidy([1, 1], posteriors, 2), [1, 1])
        self.assertEqual(call.best_dosages_summing_to_ploidy(posteriors, 2), [1, 1])

    def test_copy_is_removed_where_it_costs_least(self):
        # Sum is 3 in a diploid.  Dropping a copy of the second allele keeps the most probability.
        posteriors = [[0.01, 0.09, 0.90], [0.45, 0.55, 0.00]]
        self.assertEqual(call.correct_dosages_to_ploidy([2, 1], posteriors, 2), [2, 0])

    def test_copy_is_added_where_it_costs_least(self):
        posteriors = [[0.05, 0.60, 0.30, 0.03, 0.02], [0.10, 0.80, 0.05, 0.03, 0.02], [0.70, 0.25, 0.03, 0.01, 0.01]]
        # Sum is 2 in a tetraploid, so two copies are added one at a time.
        self.assertEqual(call.correct_dosages_to_ploidy([1, 1, 0], posteriors, 4), [2, 1, 1])

    def test_tied_dosage_gives_no_call(self):
        self.assertIsNone(call.correct_dosages_to_ploidy([None, 1], [[0.5, 0.5, 0.0], [0.1, 0.8, 0.1]], 2))

    def test_no_possible_correction_gives_no_call(self):
        # Both alleles are certainly absent, so no copy can be added.
        posteriors = [[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        self.assertIsNone(call.correct_dosages_to_ploidy([0, 0], posteriors, 2))
        self.assertIsNone(call.best_dosages_summing_to_ploidy(posteriors, 2))

    def test_exact_search_matches_brute_force(self):
        posteriors = [
            [0.02, 0.10, 0.50, 0.30, 0.05, 0.02, 0.01],
            [0.30, 0.40, 0.20, 0.05, 0.03, 0.01, 0.01],
            [0.05, 0.15, 0.45, 0.25, 0.05, 0.03, 0.02],
            [0.60, 0.30, 0.05, 0.02, 0.01, 0.01, 0.01],
        ]
        ploidy = 6

        def log_probability(dosages):
            return sum(math.log(posteriors[allele][dosage]) for allele, dosage in enumerate(dosages))

        every_genotype = [
            dosages for dosages in itertools.product(range(ploidy + 1), repeat=4) if sum(dosages) == ploidy
        ]
        self.assertEqual(len(every_genotype), 84)
        brute_force = max(every_genotype, key=log_probability)
        exact = call.best_dosages_summing_to_ploidy(posteriors, ploidy)
        self.assertEqual(exact, list(brute_force))
        greedy = call.correct_dosages_to_ploidy(
            [call.most_probable_dosage(posterior) for posterior in posteriors], posteriors, ploidy
        )
        self.assertLessEqual(log_probability(greedy), log_probability(exact) + 1e-12)


def simulated_two_allele_markers(overdispersion, ploidy, markers, samples, genotypes="hwe", seed=11):
    """Read counts ``[marker][sample][allele]`` drawn from the beta-binomial model.

    ``genotypes`` is ``hwe`` (random mating) or ``inbred`` (every sample
    homozygous).  The seed is fixed, so the counts are the same on every run.
    """
    generator = random.Random(seed)
    table = []
    for _ in range(markers):
        frequency = generator.uniform(0.15, 0.85)
        depths = []
        for _ in range(samples):
            if genotypes == "hwe":
                dosage = sum(generator.random() < frequency for _ in range(ploidy))
            else:
                dosage = ploidy if generator.random() < frequency else 0
            read_fraction = dosage / ploidy
            if 0 < read_fraction < 1:
                read_fraction = generator.betavariate(
                    read_fraction * overdispersion, (1 - read_fraction) * overdispersion
                )
            marker_depth = generator.randint(20, 120)
            alternate_depth = sum(generator.random() < read_fraction for _ in range(marker_depth))
            depths.append([marker_depth - alternate_depth, alternate_depth])
        table.append(depths)
    return table


class OverdispersionEstimateTests(unittest.TestCase):
    def test_likelihood_of_one_sample_is_the_prior_weighted_average(self):
        settings = call.ModelSettings(overdispersion=12)
        frequency = 0.3
        by_dosage = call.allele_log_likelihoods(7, 20, 4, frequency, settings)
        expected = math.log(
            sum(weight * math.exp(value) for weight, value in zip(call.hwe_prior(4, frequency), by_dosage))
        )
        # A second sample without reads adds nothing.
        observed = call.log_likelihood_of_counts([[13, 7], [0, 0]], [4, 4], frequency, settings)
        self.assertAlmostEqual(observed, expected, places=12)

    def test_recovers_the_simulated_value_in_diploids_and_polyploids(self):
        for ploidy in (2, 4, 6):
            with self.subTest(ploidy=ploidy):
                table = simulated_two_allele_markers(20, ploidy, markers=40, samples=60)
                estimate = call.estimate_overdispersion(table, [ploidy] * 60)
                self.assertIn(estimate.value, (16, 18, 20, 25))
                self.assertIsNone(estimate.at_limit)
                self.assertEqual(estimate.markers_used, 40)
                self.assertEqual(max(estimate.log_likelihood, key=estimate.log_likelihood.get), estimate.value)

    def test_inbred_samples_end_at_the_lowest_value(self):
        # Without heterozygotes the counts say nothing about scatter, and the
        # Hardy-Weinberg weights of the estimate do not fit.
        table = simulated_two_allele_markers(20, 2, markers=30, samples=60, genotypes="inbred")
        estimate = call.estimate_overdispersion(table, [2] * 60)
        self.assertEqual((estimate.value, estimate.at_limit), (min(call.OVERDISPERSION_GRID), "lowest"))

    def test_nothing_to_estimate_from(self):
        estimate = call.estimate_overdispersion([[[0, 0], [0, 0]]], [2, 2])
        self.assertIsNone(estimate.value)
        self.assertIsNone(call.estimate_overdispersion([], [2, 2]).value)

    def test_marker_subset_is_limited_and_evenly_spaced(self):
        self.assertEqual(call._evenly_spaced(list(range(10)), 20), list(range(10)))
        self.assertEqual(call._evenly_spaced(list(range(10)), 5), [0, 2, 4, 6, 8])


class PolyradParityTests(unittest.TestCase):
    """madcats_call against polyRAD on the same read counts."""

    @classmethod
    def setUpClass(cls):
        cls.fixture = load_fixture()
        cls.settings = call.ModelSettings(
            overdispersion=cls.fixture["overdispersion"],
            contamination_rate=cls.fixture["contamination_rate"],
            tolerance=cls.fixture["iteration_tolerance"],
            max_iterations=1000,
        )

    def locus_depths(self, columns):
        return [[row[column] for column in columns] for row in self.fixture["depth"]]

    def position(self, record, columns):
        """Sample index and index of the record's allele within its locus, or None for another locus."""
        column = self.fixture["allele_index"][record["allele"]]
        if column not in columns:
            return None
        return self.fixture["sample_index"][record["sample"]], columns.index(column)

    def test_single_pass_matches_polyrad(self):
        fixture = self.fixture
        sample_count = len(fixture["samples"])
        compared = 0
        for columns in fixture["locus_columns"].values():
            depths = self.locus_depths(columns)
            frequencies = call.depth_ratio_frequencies(depths)
            for column, frequency in zip(columns, frequencies):
                self.assertAlmostEqual(frequency, fixture["single_pass"]["allele_frequency"][column], places=14)
            for record in fixture["single_pass"]["records"]:
                position = self.position(record, columns)
                if position is None:
                    continue
                sample, allele = position
                ploidy = fixture["ploidy"][sample]
                log_likelihoods = call.allele_log_likelihoods(
                    depths[sample][allele],
                    sum(depths[sample]),
                    ploidy,
                    call.frequency_for_likelihood(frequencies[allele], sample_count),
                    self.settings,
                )
                prior = call.hwe_prior(ploidy, frequencies[allele])
                posterior = call.dosage_posterior(log_likelihoods, prior)
                for dosage in range(ploidy + 1):
                    expected_likelihood = record["likelihood"][dosage]
                    self.assertAlmostEqual(
                        math.exp(log_likelihoods[dosage]) / expected_likelihood, 1.0, places=9,
                        msg=f"likelihood {record['sample']} {record['allele']} dosage {dosage}",
                    )
                    self.assertAlmostEqual(prior[dosage], record["prior"][dosage], places=13)
                    self.assertAlmostEqual(posterior[dosage], record["posterior"][dosage], places=11)
                    compared += 1
        self.assertEqual(compared, sum(len(record["posterior"]) for record in fixture["single_pass"]["records"]))

    def test_iterated_frequencies_posteriors_and_genotypes_match_polyrad(self):
        fixture = self.fixture
        for columns in fixture["locus_columns"].values():
            depths = self.locus_depths(columns)
            result = call.call_marker(depths, fixture["ploidy"], self.settings)
            self.assertTrue(result.converged)
            for column, frequency in zip(columns, result.allele_frequencies):
                self.assertAlmostEqual(frequency, fixture["iterated"]["allele_frequency"][column], places=7)
            for record in fixture["iterated"]["records"]:
                position = self.position(record, columns)
                if position is None:
                    continue
                sample, allele = position
                for value, expected in zip(result.posteriors[sample][allele], record["posterior"]):
                    self.assertAlmostEqual(value, expected, places=7)
            for sample, ploidy in enumerate(fixture["ploidy"]):
                most_probable = [call.most_probable_dosage(posterior) for posterior in result.posteriors[sample]]
                corrected = call.correct_dosages_to_ploidy(most_probable, result.posteriors[sample], ploidy)
                expected = [fixture["iterated"]["genotypes"][sample][column] for column in columns]
                self.assertEqual(corrected, expected, msg=f"{fixture['samples'][sample]} columns {columns}")



# ---------------------------------------------------------------------------
# The command
# ---------------------------------------------------------------------------

# A 30-base Ref tag and tags derived from it.  Positions below are 1-based.
REF_TAG = "ACGTACGTACGGATCCTTAGCATGCAAGTC"


def with_base(tag, position, base):
    return tag[: position - 1] + base + tag[position:]


ALT_TAG = with_base(REF_TAG, 12, "C")           # the designed target: 12 G>C
REF_TAG_SITE_20 = with_base(REF_TAG, 20, "T")   # Ref at the target, G>T at position 20
ALT_TAG_SITE_20 = with_base(ALT_TAG, 20, "T")
OTHER_TAG = with_base(REF_TAG, 12, "T")         # a third base at the target


def run_madcats(arguments):
    """Run a madcats command quietly; return what it printed to stderr."""
    captured = io.StringIO()
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(captured):
        madcats.main([str(argument) for argument in arguments])
    return captured.getvalue()


def read_vcf(path):
    """Parse a VCF into (header lines, sample names, records).

    Each record is a dict of the fixed columns plus ``samples``: one dict of
    FORMAT key to value per sample.
    """
    header, records, sample_names = [], [], []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#"):
                sample_names = line.split("\t")[9:]
            else:
                fields = line.split("\t")
                keys = fields[8].split(":")
                records.append(
                    {
                        "CHROM": fields[0], "POS": int(fields[1]), "ID": fields[2], "REF": fields[3],
                        "ALT": [] if fields[4] == "." else fields[4].split(","), "INFO": fields[7], "FORMAT": keys,
                        "samples": [dict(zip(keys, column.split(":"))) for column in fields[9:]],
                    }
                )
    return header, sample_names, records


class CommandTestCase(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.dir = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def write_madc(self, rows, samples, name="madc.csv", leading_columns=()):
        """Write a short-layout MADC.  ``rows`` are (marker, category, sequence, counts)."""
        path = self.dir / name
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([*leading_columns, "AlleleID", "CloneID", "AlleleSequence", *samples])
            for marker, category, sequence, counts in rows:
                writer.writerow([*["x" for _ in leading_columns], f"{marker}|{category}", marker, sequence, *counts])
        return path

    def call(self, madc, *options):
        prefix = self.dir / "calls"
        self.stderr = run_madcats(["call", "--madc", madc, "--output-prefix", prefix, "--jobs", "1", *options])
        with open(str(prefix) + ".qc.json") as handle:
            self.qc = json.load(handle)
        return prefix

    def vcf(self, prefix, vcf_type):
        return read_vcf(str(prefix) + call.VCF_FILE_SUFFIX[vcf_type])


class TargetVcfTests(CommandTestCase):
    def test_default_writes_only_the_target_vcf(self):
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [30, 15, 0]), ("M1", "Alt", ALT_TAG, [0, 16, 28])], ["s1", "s2", "s3"]
        )
        prefix = self.call(madc)
        self.assertTrue(Path(str(prefix) + ".targets.vcf").is_file())
        self.assertFalse(Path(str(prefix) + ".all.vcf").exists())
        self.assertFalse(Path(str(prefix) + ".haplotypes.vcf").exists())
        header, samples, records = self.vcf(prefix, "target")
        self.assertEqual(samples, ["s1", "s2", "s3"])
        self.assertIn("##contig=<ID=M1,length=30>", header)
        self.assertTrue(any(line.startswith("##madcatsCoordinates=marker-local") for line in header))
        (record,) = records
        self.assertEqual((record["CHROM"], record["POS"], record["REF"], record["ALT"]), ("M1", 12, "G", ["C"]))
        self.assertEqual([sample["GT"] for sample in record["samples"]], ["0/0", "0/1", "1/1"])
        self.assertEqual(record["FORMAT"], ["GT", "AD", "DP", "XD", "GQ", "DS", "GP", "PL"])

    def test_match_rows_count_for_their_allele_and_other_rows_are_excluded(self):
        madc = self.write_madc(
            [
                ("M1", "Ref", REF_TAG, [20, 10, 0]),
                ("M1", "Alt", ALT_TAG, [0, 10, 20]),
                ("M1", "RefMatch", REF_TAG_SITE_20, [5, 0, 0]),
                ("M1", "AltMatch", ALT_TAG_SITE_20, [0, 4, 6]),
                ("M1", "Other", OTHER_TAG, [3, 0, 7]),
            ],
            ["s1", "s2", "s3"],
        )
        prefix = self.call(madc)
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual([sample["AD"] for sample in record["samples"]], ["25,0", "10,14", "0,26"])
        self.assertEqual([sample["DP"] for sample in record["samples"]], ["25", "24", "26"])
        self.assertEqual([sample["XD"] for sample in record["samples"]], ["3", "0", "7"])
        self.assertEqual(self.qc["vcfs"]["target"]["reads_excluded"], 10)
        self.assertEqual(self.qc["vcfs"]["target"]["records_with_many_excluded_reads"], 1)

    def test_excluded_read_share_filter(self):
        madc = self.write_madc(
            [
                ("M1", "Ref", REF_TAG, [20, 20]),
                ("M1", "Alt", ALT_TAG, [0, 20]),
                ("M1", "Other", OTHER_TAG, [30, 2]),
            ],
            ["s1", "s2"],
        )
        prefix = self.call(madc, "--max-other-fraction", "0.5")
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual([sample["GT"] for sample in record["samples"]], ["./.", "0/1"])
        self.assertEqual(self.qc["vcfs"]["target"]["genotypes_missing_excluded_read_share"], 1)

    def test_low_depth_and_no_reads_are_missing_but_keep_their_counts(self):
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [30, 4, 0, 12]), ("M1", "Alt", ALT_TAG, [0, 3, 0, 14])], ["s1", "s2", "s3", "s4"]
        )
        prefix = self.call(madc)
        (record,) = self.vcf(prefix, "target")[2]
        low_depth, no_reads = record["samples"][1], record["samples"][2]
        self.assertEqual((low_depth["GT"], low_depth["AD"], low_depth["DP"]), ("./.", "4,3", "7"))
        # The posterior of a low-depth sample is still reported; only the call is withheld.
        self.assertEqual(len(low_depth["GP"].split(",")), 3)
        self.assertEqual(
            no_reads, {"GT": "./.", "AD": "0,0", "DP": "0", "XD": "0", "GQ": ".", "DS": ".", "GP": ".", "PL": "."}
        )
        self.assertIn("NS=2", record["INFO"])
        self.assertEqual(self.qc["vcfs"]["target"]["genotypes_missing_low_depth"], 1)
        self.assertEqual(self.qc["vcfs"]["target"]["genotypes_missing_no_reads"], 1)

        prefix = self.call(madc, "--min-dp", "0")
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual(record["samples"][1]["GT"], "0/1")

    def test_low_genotype_quality_is_missing(self):
        # 16 of 20 reads are Ref: between a heterozygote and a Ref homozygote.
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [16, 30, 14, 0]), ("M1", "Alt", ALT_TAG, [4, 0, 15, 30])], ["s1", "s2", "s3", "s4"]
        )
        prefix = self.call(madc, "--min-gq", "0")
        relaxed = self.vcf(prefix, "target")[2][0]["samples"][0]
        self.assertNotEqual(relaxed["GT"], "./.")
        prefix = self.call(madc, "--min-gq", str(int(relaxed["GQ"]) + 1))
        strict = self.vcf(prefix, "target")[2][0]["samples"][0]
        self.assertEqual((strict["GT"], strict["GQ"]), ("./.", relaxed["GQ"]))
        self.assertEqual(self.qc["vcfs"]["target"]["genotypes_missing_low_confidence"], 1)

    def test_mixed_ploidy_from_a_ploidy_file(self):
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [40, 30, 20, 61]), ("M1", "Alt", ALT_TAG, [0, 10, 40, 59])],
            ["dip", "tet", "hex", "tri"],
        )
        ploidy_file = self.dir / "ploidy.tsv"
        ploidy_file.write_text("sample\tploidy\ntet\t4\nhex\t6\ntri\t3\nnot_in_madc\t8\n")
        # The uniform prior keeps the expected calls a matter of read fractions alone.
        prefix = self.call(madc, "--ploidy-file", ploidy_file, "--min-gq", "0", "--prior", "uniform")
        (record,) = self.vcf(prefix, "target")[2]
        genotypes = {name: sample for name, sample in zip(["dip", "tet", "hex", "tri"], record["samples"])}
        self.assertEqual(genotypes["dip"]["GT"], "0/0")
        self.assertEqual(genotypes["tet"]["GT"], "0/0/0/1")
        self.assertEqual(genotypes["hex"]["GT"], "0/0/1/1/1/1")
        self.assertEqual(len(genotypes["tri"]["GT"].split("/")), 3)
        for name, ploidy in (("dip", 2), ("tet", 4), ("hex", 6), ("tri", 3)):
            self.assertEqual(len(genotypes[name]["GP"].split(",")), ploidy + 1)
            self.assertEqual(len(genotypes[name]["PL"].split(",")), ploidy + 1)
            self.assertAlmostEqual(sum(map(float, genotypes[name]["GP"].split(","))), 1.0, places=3)
        self.assertEqual(self.qc["samples"]["by_ploidy"], {"2": 1, "3": 1, "4": 1, "6": 1})
        self.assertEqual(self.qc["samples"]["ploidy_file_samples_not_in_madc"], 1)
        self.assertIn("not in the MADC", self.stderr)

    def test_pl_comes_from_the_likelihood_and_gp_from_the_posterior(self):
        # Nearly every sample is a Ref homozygote, so the HWE prior favours Ref strongly.
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [40] * 11 + [10]), ("M1", "Alt", ALT_TAG, [0] * 11 + [1])],
            [f"s{number}" for number in range(12)],
        )
        prefix = self.call(madc, "--prior", "uniform")
        uniform = self.vcf(prefix, "target")[2][0]["samples"][-1]
        prefix = self.call(madc, "--prior", "hwe")
        hwe = self.vcf(prefix, "target")[2][0]["samples"][-1]
        # 10 Ref reads and 1 Alt read are most likely from a heterozygote, whatever the prior ...
        for sample in (uniform, hwe):
            self.assertEqual(sample["PL"].split(",")[1], "0")
        # ... but the prior decides how much posterior probability the heterozygote gets.
        heterozygote_with_uniform = float(uniform["GP"].split(",")[1])
        heterozygote_with_hwe = float(hwe["GP"].split(",")[1])
        self.assertGreater(heterozygote_with_uniform, heterozygote_with_hwe + 0.05)

    def test_target_is_not_guessed_when_the_design_is_ambiguous(self):
        two_differences = with_base(ALT_TAG, 25, "G")
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [30, 15]), ("M1", "Alt", two_differences, [0, 16])], ["s1", "s2"]
        )
        prefix = self.call(madc, "--vcf", "target", "haplotype")
        self.assertEqual(self.vcf(prefix, "target")[2], [])
        self.assertEqual(self.qc["vcfs"]["target"]["targets_skipped_ambiguous_design"], 1)
        self.assertIn("no target record", self.stderr)
        # The haplotype record does not depend on knowing the target.
        self.assertEqual(len(self.vcf(prefix, "haplotype")[2]), 1)

        tags = self.dir / "tags.csv"
        with open(tags, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(
                ["MarkerName", "LeftPrimerSequence", "RefSequence", "AltSequence", "VariantPositionInTagSequence"]
            )
            writer.writerow(["M1", REF_TAG[:8], REF_TAG, two_differences, 25])
        prefix = self.call(madc, "--reference-tags", tags)
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual((record["POS"], record["REF"], record["ALT"]), (25, "C", ["G"]))

    def test_marker_without_an_alt_row_has_no_target_record(self):
        madc = self.write_madc([("M1", "Ref", REF_TAG, [30, 15])], ["s1", "s2"])
        prefix = self.call(madc)
        self.assertEqual(self.vcf(prefix, "target")[2], [])
        self.assertEqual(self.qc["vcfs"]["target"]["targets_skipped_no_alt_row"], 1)

    def test_designed_deletion_is_written_as_an_indel(self):
        # Alt deletes bases 13-14 (AT); the fixed window pulls two further bases in at the end.
        deletion_alt = REF_TAG[:12] + REF_TAG[14:] + "TT"
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [30, 15, 0]), ("M1", "Alt", deletion_alt, [0, 16, 28])], ["s1", "s2", "s3"]
        )
        prefix = self.call(madc)
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual(len(record["REF"]) - len(record["ALT"][0]), 2)
        self.assertEqual(record["REF"], REF_TAG[record["POS"] - 1 : record["POS"] + 2])
        self.assertEqual([sample["GT"] for sample in record["samples"]], ["0/0", "0/1", "1/1"])

    def test_hapapp_layout_and_fixed_allele_ids(self):
        path = self.dir / "hapapp.csv"
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Code_version", "AlleleID", "CloneID", "AlleleSequence", "s1", "s2"])
            writer.writerow(["v1", "M1|Ref_0001", "M1", REF_TAG, "30.0", "14.0"])
            writer.writerow(["v1", "M1|Alt_0002", "M1", ALT_TAG, "0.0", "15.0"])
            writer.writerow(["v1", "M1|RefMatch_0001", "M1", REF_TAG_SITE_20, "6.0", "0.0"])
        prefix = self.call(path)
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual([sample["AD"] for sample in record["samples"]], ["36,0", "14,15"])

    def test_negative_counts_clear_the_sample_at_that_marker(self):
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [-100, 15, 40]), ("M1", "Alt", ALT_TAG, [20, 16, 0])], ["s1", "s2", "s3"]
        )
        prefix = self.call(madc)
        (record,) = self.vcf(prefix, "target")[2]
        self.assertEqual((record["samples"][0]["GT"], record["samples"][0]["DP"]), ("./.", "0"))
        self.assertEqual(self.qc["markers"]["sample_cells_cleared_for_negative_counts"], 1)
        self.assertEqual(self.qc["count_overflow"]["negative_cells"], 1)

    def test_iteration_limit_is_reported(self):
        madc = self.write_madc(
            [("M1", "Ref", REF_TAG, [9, 2, 15, 11]), ("M1", "Alt", ALT_TAG, [3, 8, 5, 10])], ["s1", "s2", "s3", "s4"]
        )
        self.call(madc, "--ploidy", "4", "--max-iterations", "1", "--tol", "1e-12")
        self.assertEqual(self.qc["vcfs"]["target"]["records_not_converged"], 1)
        self.assertIn("--max-iterations", self.stderr)


class AllVcfTests(CommandTestCase):
    def rows(self):
        return [
            ("M1", "Ref", REF_TAG, [30, 12, 0, 20]),
            ("M1", "Alt", ALT_TAG, [0, 0, 14, 0]),
            ("M1", "RefMatch", REF_TAG_SITE_20, [0, 13, 0, 0]),
            ("M1", "AltMatch", ALT_TAG_SITE_20, [0, 0, 15, 0]),
        ]

    def test_target_and_off_target_substitution(self):
        prefix = self.call(self.write_madc(self.rows(), ["s1", "s2", "s3", "s4"]), "--vcf", "all", "target")
        records = self.vcf(prefix, "all")[2]
        self.assertEqual([(record["POS"], record["ID"]) for record in records], [(12, "M1"), (20, "M1_20")])
        target, site = records
        self.assertIn("TARGET", target["INFO"].split(";"))
        self.assertNotIn("TARGET", site["INFO"].split(";"))
        self.assertEqual((site["REF"], site["ALT"]), ("G", ["T"]))
        # Reads are grouped by the base at position 20, whatever the row's base at the target.
        self.assertEqual([sample["AD"] for sample in site["samples"]], ["30,0", "12,13", "14,15", "20,0"])
        self.assertEqual([sample["GT"] for sample in site["samples"]], ["0/0", "0/1", "0/1", "0/0"])
        self.assertEqual(site["FORMAT"], ["GT", "AD", "DP", "XD", "GQ", "DS", "GP", "PL"])
        # The target record is the same line in both files.
        self.assertEqual(self.vcf(prefix, "target")[2], [target])

    def test_third_base_at_a_site_gives_a_multiallelic_record_with_mpc(self):
        rows = self.rows() + [("M1", "RefMatch", with_base(REF_TAG, 20, "A"), [0, 0, 0, 22])]
        prefix = self.call(self.write_madc(rows, ["s1", "s2", "s3", "s4"]), "--vcf", "all")
        site = self.vcf(prefix, "all")[2][1]
        # Alternate bases are ordered by read count: T has 28 reads, A has 22.
        self.assertEqual((site["REF"], site["ALT"]), ("G", ["T", "A"]))
        self.assertEqual(site["FORMAT"], ["GT", "AD", "DP", "XD", "DS", "MPC"])
        self.assertEqual(site["samples"][3]["AD"], "20,0,22")
        self.assertEqual(site["samples"][3]["GT"], "0/2")

    def test_row_that_deletes_a_site_is_excluded_there(self):
        # A RefMatch row that deletes bases 19-21, which removes position 20.
        deletion = REF_TAG[:18] + REF_TAG[21:] + "TTT"
        rows = self.rows() + [("M1", "RefMatch", deletion, [0, 0, 0, 18])]
        prefix = self.call(self.write_madc(rows, ["s1", "s2", "s3", "s4"]), "--vcf", "all")
        records = self.vcf(prefix, "all")[2]
        # The deletion itself is not a record: only substitutions are, besides the target.
        self.assertEqual([record["POS"] for record in records], [12, 20])
        target, site = records
        self.assertEqual((target["samples"][3]["AD"], target["samples"][3]["XD"]), ("38,0", "0"))
        self.assertEqual((site["samples"][3]["AD"], site["samples"][3]["XD"]), ("20,0", "18"))

    def test_base_substitution_at_the_target_is_not_a_second_record(self):
        rows = self.rows() + [("M1", "Other", OTHER_TAG, [0, 0, 0, 9])]
        prefix = self.call(self.write_madc(rows, ["s1", "s2", "s3", "s4"]), "--vcf", "all")
        records = self.vcf(prefix, "all")[2]
        self.assertEqual([record["POS"] for record in records], [12, 20])
        self.assertEqual(records[0]["samples"][3]["XD"], "9")


class HaplotypeVcfTests(CommandTestCase):
    def test_alleles_are_whole_tags_and_identical_rows_are_added_together(self):
        rows = [
            ("M1", "Ref", REF_TAG, [30, 12, 0]),
            ("M1", "Alt", ALT_TAG, [0, 0, 14]),
            ("M1", "RefMatch_0001", REF_TAG_SITE_20, [0, 6, 0]),
            ("M1", "RefMatch_0002", REF_TAG_SITE_20, [0, 7, 0]),
            ("M1", "AltMatch", ALT_TAG_SITE_20, [0, 0, 0]),
            ("M2", "Ref", REF_TAG, [25, 25, 25]),
            ("M2", "Alt", ALT_TAG, [0, 0, 0]),
        ]
        prefix = self.call(self.write_madc(rows, ["s1", "s2", "s3"]), "--vcf", "haplotype")
        header, _, records = self.vcf(prefix, "haplotype")
        first, second = records
        self.assertEqual((first["POS"], first["REF"]), (1, REF_TAG))
        # The AltMatch row has no reads in any sample, so it is not an allele.
        self.assertEqual(first["ALT"], [ALT_TAG, REF_TAG_SITE_20])
        self.assertEqual(first["FORMAT"], ["GT", "AD", "DP", "DS", "MPC"])
        self.assertEqual([sample["AD"] for sample in first["samples"]], ["30,0,0", "12,0,13", "0,14,0"])
        self.assertEqual([sample["GT"] for sample in first["samples"]], ["0/0", "0/2", "1/1"])
        # A marker where only the Ref tag has reads has no alternate allele.
        self.assertEqual((second["ALT"], [sample["GT"] for sample in second["samples"]]), ([], ["0/0"] * 3))
        self.assertFalse(any("ID=GQ" in line or "ID=GP" in line or "ID=PL" in line for line in header))

        with open(str(prefix) + ".haplotypes.alleles.csv", newline="") as handle:
            allele_table = list(csv.DictReader(handle))
        self.assertEqual(
            [(row["CloneID"], row["AlleleIndex"], row["AlleleIDs"], row["TotalReads"]) for row in allele_table],
            [
                ("M1", "0", "M1|Ref", "42"),
                ("M1", "1", "M1|Alt", "14"),
                ("M1", "2", "M1|RefMatch_0001;M1|RefMatch_0002", "13"),
                ("M2", "0", "M2|Ref", "75"),
            ],
        )

    def test_low_marginal_posterior_is_missing(self):
        rows = [
            ("M1", "Ref", REF_TAG, [16, 30, 14, 0]),
            ("M1", "Alt", ALT_TAG, [4, 0, 15, 30]),
        ]
        madc = self.write_madc(rows, ["s1", "s2", "s3", "s4"])
        prefix = self.call(madc, "--vcf", "haplotype", "--min-marginal-posterior", "0")
        relaxed = self.vcf(prefix, "haplotype")[2][0]["samples"][0]
        self.assertNotEqual(relaxed["GT"], "./.")
        self.assertLess(float(relaxed["MPC"]), 0.999)
        prefix = self.call(madc, "--vcf", "haplotype", "--min-marginal-posterior", "0.999")
        strict = self.vcf(prefix, "haplotype")[2][0]["samples"][0]
        self.assertEqual((strict["GT"], strict["MPC"]), ("./.", relaxed["MPC"]))

    def test_hexaploid_dosages_add_up_to_the_ploidy(self):
        third = with_base(REF_TAG, 5, "C")
        fourth = with_base(REF_TAG, 27, "A")
        rows = [
            ("M1", "Ref", REF_TAG, [40, 21, 33, 10, 62]),
            ("M1", "Alt", ALT_TAG, [38, 19, 0, 31, 0]),
            ("M1", "RefMatch_1", third, [22, 18, 30, 9, 0]),
            ("M1", "RefMatch_2", fourth, [19, 0, 0, 11, 58]),
        ]
        prefix = self.call(self.write_madc(rows, ["s1", "s2", "s3", "s4", "s5"]), "--vcf", "haplotype", "--ploidy", "6",
                           "--min-marginal-posterior", "0")
        (record,) = self.vcf(prefix, "haplotype")[2]
        for sample in record["samples"]:
            alleles = sample["GT"].split("/")
            self.assertEqual(len(alleles), 6)
            self.assertEqual(alleles, sorted(alleles))
        counts = self.qc["vcfs"]["haplotype"]
        self.assertGreaterEqual(counts["genotypes_called"], 4)
        self.assertLessEqual(counts.get("corrections_where_exact_search_differs", 0),
                             counts.get("genotypes_corrected_to_ploidy", 0))


class CommandBehaviourTests(CommandTestCase):
    def madc_from_simulation(self, table):
        rows = []
        for number, depths in enumerate(table):
            marker = f"M{number:03d}"
            rows.append((marker, "Ref", REF_TAG, [pair[0] for pair in depths]))
            rows.append((marker, "Alt", ALT_TAG, [pair[1] for pair in depths]))
        return self.write_madc(rows, [f"s{sample:02d}" for sample in range(len(table[0]))])

    def test_overdispersion_auto_is_estimated_and_recorded(self):
        madc = self.madc_from_simulation(simulated_two_allele_markers(20, 4, markers=40, samples=60))
        prefix = self.call(madc, "--ploidy", "4", "--overdispersion", "auto")
        report = self.qc["overdispersion"]
        self.assertEqual(report["requested"], "auto")
        self.assertIn(report["used"], (16, 18, 20, 25))
        self.assertEqual(report["used"], report["estimate"])
        self.assertEqual(report["markers_used"], 40)
        self.assertIn("Estimated overdispersion", self.stderr)
        header = self.vcf(prefix, "target")[0]
        model_line = next(line for line in header if line.startswith("##madcatsModel="))
        self.assertIn(f"overdispersion {report['used']:g}, estimated from the read counts", model_line)

        self.call(madc, "--ploidy", "4", "--overdispersion", "30")
        self.assertEqual(self.qc["overdispersion"], {"requested": 30.0, "used": 30.0})

    def test_overdispersion_auto_falls_back_for_inbred_samples(self):
        table = simulated_two_allele_markers(20, 2, markers=30, samples=60, genotypes="inbred")
        prefix = self.call(self.madc_from_simulation(table), "--overdispersion", "auto")
        report = self.qc["overdispersion"]
        self.assertEqual((report["at_limit"], report["used"]), ("lowest", call.DEFAULT_OVERDISPERSION))
        self.assertIn("Hardy-Weinberg", self.stderr)
        model_line = next(line for line in self.vcf(prefix, "target")[0] if line.startswith("##madcatsModel="))
        self.assertNotIn("estimated from the read counts", model_line)

    def test_parallel_run_writes_the_same_files(self):
        rows = []
        for number in range(12):
            marker = f"M{number:02d}"
            rows += [
                (marker, "Ref", REF_TAG, [30 + number, 12, 0, 20]),
                (marker, "Alt", ALT_TAG, [0, 9 + number, 14, 0]),
                (marker, "RefMatch", REF_TAG_SITE_20, [0, 13, 0, number]),
            ]
        madc = self.write_madc(rows, ["s1", "s2", "s3", "s4"])
        outputs = {}
        for jobs in ("1", "3"):
            prefix = self.dir / f"jobs{jobs}"
            run_madcats(["call", "--madc", madc, "--output-prefix", prefix, "--vcf", "target", "all", "haplotype",
                         "--jobs", jobs])
            outputs[jobs] = [Path(str(prefix) + suffix).read_text() for suffix in call.VCF_FILE_SUFFIX.values()]
        self.assertEqual(outputs["1"], outputs["3"])

    def test_repeated_sample_names_are_an_error(self):
        madc = self.write_madc([("M1", "Ref", REF_TAG, [1, 2]), ("M1", "Alt", ALT_TAG, [1, 2])], ["s1", "s1"])
        with self.assertRaises(SystemExit), contextlib.redirect_stderr(io.StringIO()) as stderr:
            madcats.main(["call", "--madc", str(madc), "--output-prefix", str(self.dir / "x")])
        self.assertIn("repeated sample names", stderr.getvalue())

    def test_ploidy_file_errors(self):
        for content, message in (("s1\tfour\n", "not a whole number"), ("s1\t4\ns1\t2\n", "more than once")):
            path = self.dir / "ploidy.tsv"
            path.write_text("sample\tploidy\n" + content)
            with self.assertRaises(madcats.MadcError) as raised:
                call.read_ploidy_table(path)
            self.assertIn(message, str(raised.exception))
        path.write_text("s1,4\ns2,6\n")
        self.assertEqual(call.read_ploidy_table(path), {"s1": 4, "s2": 6})


@unittest.skipUnless(EXAMPLE.is_dir(), "example data not present")
class ExampleProjectCallTests(unittest.TestCase):
    """Genotypes called from the bundled example against its simulated truth."""

    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory()
        out = Path(cls.temporary.name)
        tags, madc, cls.prefix = out / "reference_tags.csv", out / "demo_MADC.csv", out / "demo"
        run_madcats(["design", "--primers", EXAMPLE / "primers.csv", "--probe-design", EXAMPLE / "probe_design.tsv",
                     "--output", tags])
        run_madcats(["generate", "--fastq-dir", EXAMPLE / "fastq", "--reference-tags", tags,
                     "--sample-lookup", EXAMPLE / "sample_names.csv", "--jobs", "1", "--output", madc])
        run_madcats(["call", "--madc", madc, "--output-prefix", cls.prefix, "--vcf", "target", "all", "haplotype",
                     "--jobs", "1"])
        with open(EXAMPLE / "truth_haplotypes.csv", newline="") as handle:
            cls.truth_sequence = {
                (row["MarkerName"], row["Haplotype"]): row["TagSequence"] for row in csv.DictReader(handle)
            }
        with open(EXAMPLE / "truth_genotypes.csv", newline="") as handle:
            cls.truth_genotypes = list(csv.DictReader(handle))
        with open(madc, newline="") as handle:
            cls.category_of_sequence = {
                (row["CloneID"], row["AlleleSequence"]): row["AlleleID"].split("|")[1] for row in csv.DictReader(handle)
            }

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def called_genotypes(self, vcf_type):
        _, samples, records = read_vcf(str(self.prefix) + call.VCF_FILE_SUFFIX[vcf_type])
        return samples, {record["ID"]: record for record in records}

    def test_haplotype_genotypes_match_the_truth(self):
        samples, records = self.called_genotypes("haplotype")
        self.assertEqual(len(self.truth_genotypes), 96)
        for truth in self.truth_genotypes:
            record = records[truth["MarkerName"]]
            alleles = [record["REF"]] + record["ALT"]
            expected = sorted(
                alleles.index(self.truth_sequence[(truth["MarkerName"], haplotype)])
                for haplotype in truth["Genotype"].split("/")
            )
            called = record["samples"][samples.index(truth["Sample"])]["GT"]
            self.assertEqual(called, "/".join(map(str, expected)), msg=f"{truth['Sample']} {truth['MarkerName']}")

    def test_target_genotypes_match_the_truth(self):
        samples, records = self.called_genotypes("target")
        for truth in self.truth_genotypes:
            marker = truth["MarkerName"]
            # A haplotype carries the alternate target allele when its MADC row is Alt or AltMatch.
            expected = sorted(
                int(self.category_of_sequence[(marker, self.truth_sequence[(marker, haplotype)])].startswith("Alt"))
                for haplotype in truth["Genotype"].split("/")
            )
            called = records[marker]["samples"][samples.index(truth["Sample"])]["GT"]
            self.assertEqual(called, "/".join(map(str, expected)), msg=f"{truth['Sample']} {marker}")

    def test_all_vcf_adds_the_two_simulated_second_snps(self):
        _, records = self.called_genotypes("all")
        off_target = sorted(
            record_id for record_id, record in records.items() if "TARGET" not in record["INFO"].split(";")
        )
        self.assertEqual(off_target, ["Chr02_001257730_67", "Chr04_002318006_25"])
        self.assertEqual(len(records), 10)


if __name__ == "__main__":
    unittest.main()
