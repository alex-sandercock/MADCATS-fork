"""Tests of the genotype-calling engine in madcats_call.py.

The parity tests compare against values computed by polyRAD and stored in
tests/fixtures/polyrad_parity.json (made by scripts/make_polyrad_fixture.R),
so they run without R.
"""

import itertools
import json
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import madcats_call as call  # noqa: E402

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


if __name__ == "__main__":
    unittest.main()
