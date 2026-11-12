"""``madcats call``: genotype calls and VCFs from a MADC count table.

The command reads a MADC (from ``generate``, a vendor delivery or HapApp),
estimates each sample's genotype at every marker, and writes one VCF per
requested type:

``target``     one record per marker for the variant the panel was designed
               for.  Reads of the Ref and RefMatch rows count for the
               reference allele, reads of the Alt and AltMatch rows for the
               alternate allele, and reads of Other rows are left out.
``all``        the target records plus one record for every other position
               of the tag where a row shows a base substitution.
``haplotype``  one record per marker whose alleles are the complete tag
               sequences (microhaplotypes).

Coordinates are marker-local: CHROM is the marker (CloneID) and POS is a
1-based position on its Ref tag.

The statistical model is adapted from polyRAD (Clark, Lipka and Sacks 2019,
G3 9:663-673; https://github.com/lvclark/polyRAD, GPL >= 2).  The ``norm``
prior is the one of updog (Gerard and Ferrao 2020, Bioinformatics
36:1795-1800; https://github.com/dcgerard/updog, GPL >= 3).  In a record,
each allele is considered on its own against all the other alleles together,
which turns a multiallelic record into one two-outcome problem per allele:

1. Likelihood.  For a sample of ploidy P, an allele can be present in
   0, 1, ..., P copies (its *dosage*).  The probability of the observed number
   of reads of the allele, out of all the sample's reads in the record, is
   beta-binomial.  The expected read fraction for a dosage is dosage / P,
   moved slightly toward the population allele frequency to allow for
   contamination between samples.
2. Prior.  Either every dosage is equally probable (``uniform``), or dosages
   follow Hardy-Weinberg proportions for the allele's frequency (``hwe``):
   each of the P chromosome copies carries the allele independently.  That is
   the expectation for a diploid or an autopolyploid in a randomly mating
   population.  It is not the expectation under disomic inheritance in an
   allopolyploid, or in a family.  The ``norm`` prior makes no such
   assumption: the dosage probabilities follow a normal shape whose mean and
   spread are fitted to the samples of each ploidy.
3. Posterior.  Prior times likelihood, normalised over the dosages.
4. Iteration (``hwe`` and ``norm``).  The allele frequencies, and for
   ``norm`` the prior's two parameters, are re-estimated from the posteriors
   and steps 1-3 repeat until they stop changing.  Records do not share any parameter, so each converges on
   its own (polyRAD stops all loci together, on their mean change).

The posteriors are *marginal*: one distribution per allele, not a joint
distribution over whole genotypes.  With exactly two alleles the two
marginals mirror each other, so the posterior of the alternate allele is the
full genotype posterior and the record carries GQ, GP and PL.  With more
alleles the most probable dosages need not add up to the ploidy;
``correct_dosages_to_ploidy`` (polyRAD's method) reconciles them, and the
record carries MPC, the lowest marginal posterior among the called dosages,
in place of GQ, GP and PL.

Likelihoods are kept as logarithms and only normalised posteriors are
exponentiated, so a dosage that the reads make very improbable keeps a
usable value at any read depth.

Reading map: the model (``allele_log_likelihoods``, ``hwe_prior``,
``normal_prior``, ``fit_normal_prior``, ``dosage_posterior``, ``call_marker``); from marginal posteriors to one
genotype (``correct_dosages_to_ploidy``, ``best_dosages_summing_to_ploidy``);
the overdispersion estimate (``estimate_overdispersion``); input
(``read_count_table``, ``read_ploidy_table``); what each VCF type counts
(``resolve_target``, ``RowReading``, ``process_marker``); VCF text
(``format_record``, ``vcf_header``); the command line (``run_call``,
``add_call_command``).

This module uses only the Python standard library.  The MADC readers and the
tag alignment it shares with ``generate`` stay in ``madcats``.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import datetime as dt
import functools
import math
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

from madcats import (
    GAP_BAND_DEFAULT,
    IUPAC_BASES,
    CountOverflowCheck,
    Locus,
    MadcError,
    __version__,
    _default_jobs,
    _designed_events,
    _write_json,
    allele_type,
    find_madc_header,
    load_reference_tags,
    madc_identity_columns,
    madc_read_counts,
    madc_sample_start,
    nonnegative_int,
    open_text,
    positive_int,
    variant_events,
)

PRIOR_UNIFORM = "uniform"
PRIOR_HWE = "hwe"
PRIOR_NORM = "norm"
PRIOR_CHOICES = (PRIOR_UNIFORM, PRIOR_HWE, PRIOR_NORM)

# polyRAD's defaults: AddGenotypeLikelihood(overdispersion = 9), readDArTag(contamRate = 0.001)
# and IterateHWE(tol = 1e-5).  polyRAD has no iteration limit.  The command estimates the
# overdispersion unless one is given, and uses this default only when it cannot.
DEFAULT_OVERDISPERSION = 9.0
DEFAULT_CONTAMINATION_RATE = 0.001
DEFAULT_TOLERANCE = 1e-5
DEFAULT_MAX_ITERATIONS = 100

NEGATIVE_INFINITY = float("-inf")


@dataclasses.dataclass(frozen=True)
class ModelSettings:
    """Parameters of the genotype model.

    overdispersion: how tightly read fractions follow the expected fraction.
        The two beta-binomial shape parameters always add up to this value;
        larger values mean less spread, and the distribution approaches the
        binomial as it grows.
    contamination_rate: fraction of a sample's reads expected to come from
        other samples.  Those reads carry an allele at its population frequency.
    prior: ``uniform``, ``hwe`` or ``norm`` (see ``normal_prior``).
    selfing_rate: with the ``hwe`` prior, the share of individuals that come
        from self-fertilisation, from 0 (random mating) to 1 (fully inbred).
    tolerance: the ``hwe`` iteration stops once the mean absolute change in
        the marker's allele frequencies is at or below this.
    max_iterations: upper limit on ``hwe`` iterations.
    """

    overdispersion: float = DEFAULT_OVERDISPERSION
    contamination_rate: float = DEFAULT_CONTAMINATION_RATE
    prior: str = PRIOR_HWE
    selfing_rate: float = 0.0
    tolerance: float = DEFAULT_TOLERANCE
    max_iterations: int = DEFAULT_MAX_ITERATIONS

    def __post_init__(self) -> None:
        if self.overdispersion <= 0:
            raise ValueError("overdispersion must be positive")
        if not 0 <= self.contamination_rate < 1:
            raise ValueError("contamination_rate must be at least 0 and below 1")
        if self.prior not in PRIOR_CHOICES:
            raise ValueError(f"prior must be one of {', '.join(PRIOR_CHOICES)}")
        if not 0 <= self.selfing_rate <= 1:
            raise ValueError("selfing_rate must be between 0 and 1")
        if self.tolerance <= 0:
            raise ValueError("tolerance must be positive")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be at least 1")


@dataclasses.dataclass
class MarkerCall:
    """The model's result for one marker.

    ``log_likelihoods`` and ``posteriors`` are indexed
    ``[sample][allele][dosage]``, with dosage running from 0 to the sample's
    ploidy.  ``allele_frequencies`` are the final estimates; the likelihoods
    and posteriors were computed from the estimates of the round before, which
    differ from the final ones by ``final_change`` on average.

    ``priors`` holds the prior the posteriors were computed with, as
    ``[ploidy][allele][dosage]``.  With the ``norm`` prior,
    ``normal_prior_parameters`` holds the last fitted
    ``(mean, standard deviation)`` as ``[ploidy][allele]``, for the ploidy
    groups large enough to be fitted.

    ``callable`` is False when no sample has a read at the marker: there is
    then no allele frequency to start from, and the other fields are empty.
    """

    callable: bool
    allele_frequencies: List[float]
    log_likelihoods: List[List[List[float]]]
    posteriors: List[List[List[float]]]
    iterations: int = 0
    final_change: Optional[float] = None
    converged: bool = True
    priors: Dict[int, List[List[float]]] = dataclasses.field(default_factory=dict)
    normal_prior_parameters: Dict[int, List[Tuple[float, float]]] = dataclasses.field(default_factory=dict)


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------


def log_binomial_coefficient(total: int, chosen: int) -> float:
    """Natural log of the number of ways to choose ``chosen`` items out of ``total``."""
    return math.lgamma(total + 1) - math.lgamma(chosen + 1) - math.lgamma(total - chosen + 1)


def minimum_allele_frequency(sample_count: int) -> float:
    """The smallest allele frequency the likelihood is given.

    An estimated frequency of exactly 0 or 1 would make the contamination
    term vanish for that allele.  polyRAD replaces 0 with one allele copy in
    a diploid population of this many samples, and 1 with its complement.
    """
    return 1.0 / (2.0 * sample_count)


def frequency_for_likelihood(allele_frequency: float, sample_count: int) -> float:
    """``allele_frequency`` with exact 0 and 1 moved inside the open interval."""
    floor = minimum_allele_frequency(sample_count)
    if allele_frequency == 0.0:
        return floor
    if allele_frequency == 1.0:
        return 1.0 - floor
    return allele_frequency


def expected_read_fraction(dosage: int, ploidy: int, allele_frequency: float, contamination_rate: float) -> float:
    """Probability that one read of a sample shows the allele.

    A read comes from the sample itself with probability
    1 - contamination_rate, and then shows the allele in proportion to its
    dosage.  Otherwise it comes from another sample and shows the allele at
    its population frequency.
    """
    return (dosage / ploidy) * (1.0 - contamination_rate) + contamination_rate * allele_frequency


def allele_log_likelihoods(
    allele_depth: int,
    marker_depth: int,
    ploidy: int,
    allele_frequency: float,
    settings: ModelSettings,
) -> List[float]:
    """Log-likelihood of each dosage (0 to ``ploidy``) of one allele in one sample.

    ``allele_depth`` reads out of ``marker_depth`` showed the allele.  The
    likelihood is the beta-binomial probability of that count, with shape
    parameters ``fraction * overdispersion`` and
    ``(1 - fraction) * overdispersion``, where ``fraction`` is the expected
    read fraction of the dosage.  ``allele_frequency`` must already be inside
    the open interval (see ``frequency_for_likelihood``).

    With no reads at the marker every dosage has likelihood 1 (log 0.0).
    """
    if marker_depth == 0:
        return [0.0] * (ploidy + 1)
    other_depth = marker_depth - allele_depth
    overdispersion = settings.overdispersion
    # The two shape parameters add up to the overdispersion for every dosage, so these
    # terms of the beta-binomial formula are the same for all dosages.
    constant = (
        log_binomial_coefficient(marker_depth, allele_depth)
        + math.lgamma(overdispersion)
        - math.lgamma(marker_depth + overdispersion)
    )
    log_likelihoods = []
    for dosage in range(ploidy + 1):
        fraction = expected_read_fraction(dosage, ploidy, allele_frequency, settings.contamination_rate)
        if fraction <= 0.0 or fraction >= 1.0:
            # Only without contamination: a dosage of 0 cannot produce a read of the
            # allele, and a dosage equal to the ploidy cannot produce any other read.
            impossible = allele_depth > 0 if fraction <= 0.0 else other_depth > 0
            log_likelihoods.append(NEGATIVE_INFINITY if impossible else 0.0)
            continue
        shape_allele = fraction * overdispersion
        shape_other = (1.0 - fraction) * overdispersion
        log_likelihoods.append(
            constant
            + math.lgamma(allele_depth + shape_allele)
            + math.lgamma(other_depth + shape_other)
            - math.lgamma(shape_allele)
            - math.lgamma(shape_other)
        )
    return log_likelihoods


# ---------------------------------------------------------------------------
# Priors and posteriors
# ---------------------------------------------------------------------------


def uniform_prior(ploidy: int) -> List[float]:
    """Every dosage from 0 to ``ploidy`` equally probable."""
    return [1.0 / (ploidy + 1)] * (ploidy + 1)


def gamete_dosage_probabilities(parent_dosage: int, ploidy: int) -> List[float]:
    """Probability that a gamete carries 0, 1, ..., ploidy / 2 copies of an allele.

    The gamete receives half of the parent's chromosome copies, drawn without
    replacement from the ``parent_dosage`` copies that carry the allele and
    the rest that do not (polysomic inheritance without double reduction).
    """
    gamete_size = ploidy // 2
    return [
        math.comb(parent_dosage, copies) * math.comb(ploidy - parent_dosage, gamete_size - copies)
        / math.comb(ploidy, gamete_size)
        for copies in range(gamete_size + 1)
    ]


def _solve_linear_system(matrix: Sequence[Sequence[float]], right_hand_side: Sequence[float]) -> List[float]:
    """Solve ``matrix @ x = right_hand_side`` by Gaussian elimination with row pivoting."""
    size = len(right_hand_side)
    rows = [list(matrix[index]) + [right_hand_side[index]] for index in range(size)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda index: abs(rows[index][column]))
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for index in range(column + 1, size):
            factor = rows[index][column] / rows[column][column]
            for position in range(column, size + 1):
                rows[index][position] -= factor * rows[column][position]
    solution = [0.0] * size
    for index in range(size - 1, -1, -1):
        already_solved = sum(rows[index][position] * solution[position] for position in range(index + 1, size))
        solution[index] = (rows[index][size] - already_solved) / rows[index][index]
    return solution


@functools.lru_cache(maxsize=None)
def _selfing_system(ploidy: int, selfing_rate: float) -> Tuple[Tuple[float, ...], ...]:
    """The matrix ``I - selfing_rate * S`` of ``hwe_prior``'s selfing adjustment.

    ``S[offspring][parent]`` is the probability that self-fertilising a
    parent of one dosage gives an offspring of another: the offspring's
    dosage is the sum of two gametes of the same parent.
    """
    size = ploidy + 1
    selfed = [[0.0] * size for _ in range(size)]
    for parent_dosage in range(size):
        gamete = gamete_dosage_probabilities(parent_dosage, ploidy)
        for first_copies, first_probability in enumerate(gamete):
            for second_copies, second_probability in enumerate(gamete):
                selfed[first_copies + second_copies][parent_dosage] += first_probability * second_probability
    return tuple(
        tuple((1.0 if row == column else 0.0) - selfing_rate * selfed[row][column] for column in range(size))
        for row in range(size)
    )


def hwe_prior(ploidy: int, allele_frequency: float, selfing_rate: float = 0.0) -> List[float]:
    """Dosage probabilities under Hardy-Weinberg proportions, optionally with self-fertilisation.

    With random mating each of the ``ploidy`` chromosome copies carries the
    allele independently with probability ``allele_frequency``, so the
    dosage is binomial.

    With a ``selfing_rate`` between 0 and 1, a share of each generation comes
    from self-fertilisation and the rest from random mating.  The genotype
    proportions that such a population settles at are those of equation 6 of
    de Silva et al. (2005, Heredity 95:327-334), as in polyRAD:
    ``(1 - s) * (I - s * S)^-1 * binomial``.  Selfing moves probability from
    the intermediate dosages to 0 and ``ploidy``.  With a rate of 1 only
    those two remain.  As in polyRAD, an odd ploidy is left at random-mating
    proportions unless the rate is 1, since its gametes cannot carry half
    the copies.
    """
    other_frequency = 1.0 - allele_frequency
    if selfing_rate >= 1.0:
        return [other_frequency] + [0.0] * (ploidy - 1) + [allele_frequency]
    random_mating = [
        math.comb(ploidy, dosage) * allele_frequency**dosage * other_frequency ** (ploidy - dosage)
        for dosage in range(ploidy + 1)
    ]
    if selfing_rate <= 0.0 or ploidy % 2 != 0:
        return random_mating
    settled = _solve_linear_system(_selfing_system(ploidy, selfing_rate), random_mating)
    return [(1.0 - selfing_rate) * probability for probability in settled]


# Limits on the two parameters of the normal-shaped prior, as in updog.  The largest
# standard deviation makes the prior practically uniform.
NORMAL_PRIOR_MIN_SD = 1e-8
NORMAL_PRIOR_MAX_SD = 1000.0
NORMAL_PRIOR_MEAN_MARGIN = 1.0     # the mean may lie this far below 0 or above the ploidy
# A ploidy group smaller than this uses the Hardy-Weinberg prior instead: two prior
# parameters cannot be estimated from a handful of samples.
NORMAL_PRIOR_MIN_SAMPLES = 10


def _log_normal_prior(ploidy: int, mean: float, standard_deviation: float) -> List[float]:
    """Logarithms of ``normal_prior``; finite even where the probability itself rounds to zero."""
    exponents = [-((dosage - mean) ** 2) / (2.0 * standard_deviation**2) for dosage in range(ploidy + 1)]
    largest = max(exponents)
    log_total = largest + math.log(sum(math.exp(exponent - largest) for exponent in exponents))
    return [exponent - log_total for exponent in exponents]


def normal_prior(ploidy: int, mean: float, standard_deviation: float) -> List[float]:
    """Dosage probabilities proportional to a normal density at 0, 1, ..., ``ploidy``.

    This is updog's ``norm`` genotype distribution (Gerard and Ferrao 2020,
    Bioinformatics 36:1795-1800).  With two parameters it can be narrower or
    wider than Hardy-Weinberg proportions, so it does not assume random mating.
    """
    return [math.exp(log_probability) for log_probability in _log_normal_prior(ploidy, mean, standard_deviation)]


def weighted_log_normal_prior(weights: Sequence[float], ploidy: int, mean: float, standard_deviation: float) -> float:
    """Sum of ``weight * log(prior)`` over the dosages: what ``fit_normal_prior`` maximises."""
    log_prior = _log_normal_prior(ploidy, mean, standard_deviation)
    return sum(weight * log_probability for weight, log_probability in zip(weights, log_prior) if weight > 0.0)


def _golden_section_maximum(function, low: float, high: float, tolerance: float = 1e-10) -> float:
    """The argument between ``low`` and ``high`` at which a single-peaked function is highest."""
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    inner_low, inner_high = high - ratio * (high - low), low + ratio * (high - low)
    value_low, value_high = function(inner_low), function(inner_high)
    while high - low > tolerance:
        if value_low < value_high:
            low, inner_low, value_low = inner_low, inner_high, value_high
            inner_high = low + ratio * (high - low)
            value_high = function(inner_high)
        else:
            high, inner_high, value_high = inner_high, inner_low, value_low
            inner_low = high - ratio * (high - low)
            value_low = function(inner_low)
    return (low + high) / 2.0


def _normal_prior_newton_fit(shares: Sequence[float], ploidy: int, mean: float, standard_deviation: float):
    """Fit the normal-shaped prior without limits on its parameters, by Newton's method.

    The prior is written as ``log probability = a * dosage + b * dosage**2 -
    constant``, with ``a = mean / sd**2`` and ``b = -1 / (2 * sd**2)``.  In
    ``a`` and ``b`` the weighted log prior has a single peak, reached where
    the prior's mean dosage and mean squared dosage equal those of
    ``shares``.  Returns ``(a, b)``, or None when the iteration does not
    settle (the shares sit on one or two dosages and the peak is at infinity).
    A step that cannot improve the fit any further, to within rounding, means
    the peak has been reached.
    """
    dosages = range(ploidy + 1)
    target_first = sum(share * dosage for share, dosage in zip(shares, dosages))
    target_second = sum(share * dosage**2 for share, dosage in zip(shares, dosages))
    linear, quadratic = mean / standard_deviation**2, -1.0 / (2.0 * standard_deviation**2)

    def prior_and_objective(linear_term: float, quadratic_term: float):
        exponents = [linear_term * dosage + quadratic_term * dosage**2 for dosage in dosages]
        largest = max(exponents)
        unnormalised = [math.exp(exponent - largest) for exponent in exponents]
        total = sum(unnormalised)
        objective = linear_term * target_first + quadratic_term * target_second - largest - math.log(total)
        return [value / total for value in unnormalised], objective

    prior, objective = prior_and_objective(linear, quadratic)
    for _ in range(100):
        moments = [
            sum(probability * dosage**power for probability, dosage in zip(prior, dosages)) for power in (1, 2, 3, 4)
        ]
        gradient_first, gradient_second = target_first - moments[0], target_second - moments[1]
        # Variances and covariance of dosage and squared dosage under the current prior.
        variance_first = moments[1] - moments[0] ** 2
        variance_second = moments[3] - moments[1] ** 2
        covariance = moments[2] - moments[0] * moments[1]
        determinant = variance_first * variance_second - covariance**2
        if determinant <= 1e-300:
            return None
        step_linear = (variance_second * gradient_first - covariance * gradient_second) / determinant
        step_quadratic = (variance_first * gradient_second - covariance * gradient_first) / determinant
        scale = 1.0
        while scale > 1e-8:
            candidate = prior_and_objective(linear + scale * step_linear, quadratic + scale * step_quadratic)
            if candidate[1] >= objective:
                break
            scale /= 2.0
        else:
            return linear, quadratic
        linear, quadratic = linear + scale * step_linear, quadratic + scale * step_quadratic
        improvement = candidate[1] - objective
        prior, objective = candidate
        if improvement < 1e-14:
            return linear, quadratic
    return None


def fit_normal_prior(
    dosage_weights: Sequence[float], ploidy: int, mean: float, standard_deviation: float
) -> Tuple[float, float]:
    """The mean and standard deviation of the normal-shaped prior that best fit dosage weights.

    ``dosage_weights[d]`` is the summed posterior probability of dosage ``d``
    over the samples: the expected number of samples with that dosage.  The
    fit maximises ``sum(weight * log(prior))``, which is the prior's update
    in an EM iteration, and is exact because the prior depends on the samples
    only through these sums.  ``mean`` and ``standard_deviation`` are the
    starting values, usually the previous fit.

    The parameters are kept within updog's limits: a mean from -1 to
    ploidy + 1 and a standard deviation up to 1000.  The search first ignores
    the limits.  If the best fit then lies outside them, the best fit within
    them lies on their edge, because the function being maximised has a
    single peak; the two edges for the mean are searched along the standard
    deviation.  (The third edge, the largest standard deviation, is a
    practically uniform prior and is where those searches end if nothing
    fits better.)
    """
    total = sum(dosage_weights)
    shares = [weight / total for weight in dosage_weights]
    lowest_mean, highest_mean = -NORMAL_PRIOR_MEAN_MARGIN, ploidy + NORMAL_PRIOR_MEAN_MARGIN

    def objective(candidate_mean: float, candidate_sd: float) -> float:
        return weighted_log_normal_prior(shares, ploidy, candidate_mean, candidate_sd)

    unlimited = _normal_prior_newton_fit(shares, ploidy, mean, standard_deviation)
    if unlimited is not None and unlimited[1] < 0.0:
        linear, quadratic = unlimited
        fitted_sd = math.sqrt(-1.0 / (2.0 * quadratic))
        fitted_mean = -linear / (2.0 * quadratic)
        if lowest_mean <= fitted_mean <= highest_mean and NORMAL_PRIOR_MIN_SD <= fitted_sd <= NORMAL_PRIOR_MAX_SD:
            return fitted_mean, fitted_sd

    # The starting values are a candidate too, so that a fit is never worse than the one before.
    start_mean = min(highest_mean, max(lowest_mean, mean))
    start_sd = min(NORMAL_PRIOR_MAX_SD, max(NORMAL_PRIOR_MIN_SD, standard_deviation))
    candidates = [(objective(start_mean, start_sd), start_mean, start_sd)]
    for edge_mean in (lowest_mean, highest_mean):
        log_sd = _golden_section_maximum(
            lambda value: objective(edge_mean, math.exp(value)),
            math.log(NORMAL_PRIOR_MIN_SD), math.log(NORMAL_PRIOR_MAX_SD),
        )
        candidates.append((objective(edge_mean, math.exp(log_sd)), edge_mean, math.exp(log_sd)))
    _, best_mean, best_sd = max(candidates)
    return best_mean, best_sd


def dosage_posterior(log_likelihoods: Sequence[float], prior: Sequence[float]) -> List[float]:
    """Posterior probability of each dosage: prior times likelihood, normalised.

    When the prior gives zero probability to every dosage the reads allow,
    the product is zero everywhere and cannot be normalised; the prior is
    returned unchanged, as polyRAD does.
    """
    log_products = [
        log_likelihood + math.log(probability) if probability > 0.0 else NEGATIVE_INFINITY
        for log_likelihood, probability in zip(log_likelihoods, prior)
    ]
    largest = max(log_products)
    if largest == NEGATIVE_INFINITY:
        return list(prior)
    # Subtract the largest term before exponentiating so that it becomes exp(0) = 1 and
    # the sum cannot underflow to zero.
    scaled = [math.exp(log_product - largest) for log_product in log_products]
    total = sum(scaled)
    return [value / total for value in scaled]


def posterior_mean_dosage(posterior: Sequence[float]) -> float:
    """Expected number of copies of the allele."""
    return sum(dosage * probability for dosage, probability in enumerate(posterior))


# ---------------------------------------------------------------------------
# One marker
# ---------------------------------------------------------------------------


def depth_ratio_frequencies(depths: Sequence[Sequence[int]]) -> Optional[List[float]]:
    """Starting allele frequencies: each allele's mean share of a sample's reads.

    ``depths`` is ``[sample][allele]``.  Samples with no reads at the marker
    are left out.  Returns None when no sample has a read.
    """
    allele_count = len(depths[0])
    share_totals = [0.0] * allele_count
    samples_with_reads = 0
    for sample_depths in depths:
        marker_depth = sum(sample_depths)
        if marker_depth == 0:
            continue
        samples_with_reads += 1
        for allele, allele_depth in enumerate(sample_depths):
            share_totals[allele] += allele_depth / marker_depth
    if samples_with_reads == 0:
        return None
    return [total / samples_with_reads for total in share_totals]


def _frequency_priors(
    ploidy_groups: Sequence[int], allele_frequencies: Sequence[float], settings: ModelSettings
) -> Dict[int, List[List[float]]]:
    """Priors that follow from the allele frequencies alone, as ``[ploidy][allele][dosage]``.

    Uniform for the ``uniform`` prior and Hardy-Weinberg proportions
    otherwise.  The ``norm`` prior starts from the latter and is then fitted
    (``_fitted_normal_priors``).  Samples of the same ploidy share a prior.
    """
    if settings.prior == PRIOR_UNIFORM:
        return {ploidy: [uniform_prior(ploidy)] * len(allele_frequencies) for ploidy in ploidy_groups}
    return {
        ploidy: [hwe_prior(ploidy, frequency, settings.selfing_rate) for frequency in allele_frequencies]
        for ploidy in ploidy_groups
    }


def _fitted_normal_priors(
    posteriors: Sequence[Sequence[Sequence[float]]],
    ploidies: Sequence[int],
    allele_frequencies: Sequence[float],
    previous_parameters: Dict[int, List[Tuple[float, float]]],
) -> Tuple[Dict[int, List[List[float]]], Dict[int, List[Tuple[float, float]]]]:
    """Normal-shaped priors fitted to the current posteriors, one per ploidy group and allele.

    Returns the priors as ``[ploidy][allele][dosage]`` and their
    ``(mean, standard deviation)`` as ``[ploidy][allele]``.  Each ploidy
    group is fitted on its own samples, because a normal shape over the
    dosages of one ploidy says nothing about another.  A group with fewer
    than ``NORMAL_PRIOR_MIN_SAMPLES`` samples gets the Hardy-Weinberg prior
    for the marker's allele frequencies instead, and no parameters.
    """
    priors: Dict[int, List[List[float]]] = {}
    parameters: Dict[int, List[Tuple[float, float]]] = {}
    for ploidy in set(ploidies):
        group = [
            sample_posteriors for sample_posteriors, sample_ploidy in zip(posteriors, ploidies) if sample_ploidy == ploidy
        ]
        if len(group) < NORMAL_PRIOR_MIN_SAMPLES:
            priors[ploidy] = [hwe_prior(ploidy, frequency) for frequency in allele_frequencies]
            continue
        priors[ploidy], parameters[ploidy] = [], []
        for allele, frequency in enumerate(allele_frequencies):
            dosage_weights = [
                sum(sample_posteriors[allele][dosage] for sample_posteriors in group) for dosage in range(ploidy + 1)
            ]
            if ploidy in previous_parameters:
                start_mean, start_sd = previous_parameters[ploidy][allele]
            else:
                # The first fit starts from the binomial's mean and spread.
                start_mean = ploidy * frequency
                start_sd = max(0.1, math.sqrt(ploidy * frequency * (1.0 - frequency)))
            mean, standard_deviation = fit_normal_prior(dosage_weights, ploidy, start_mean, start_sd)
            parameters[ploidy].append((mean, standard_deviation))
            priors[ploidy].append(normal_prior(ploidy, mean, standard_deviation))
    return priors, parameters


def _marker_pass(
    depths: Sequence[Sequence[int]],
    ploidies: Sequence[int],
    allele_frequencies: Sequence[float],
    priors: Dict[int, List[List[float]]],
    settings: ModelSettings,
) -> Tuple[List[List[List[float]]], List[List[List[float]]]]:
    """Likelihoods and posteriors of every sample and allele for given allele frequencies and priors."""
    sample_count = len(depths)
    likelihood_frequencies = [frequency_for_likelihood(frequency, sample_count) for frequency in allele_frequencies]
    log_likelihoods = []
    posteriors = []
    for sample_depths, ploidy in zip(depths, ploidies):
        marker_depth = sum(sample_depths)
        sample_log_likelihoods = [
            allele_log_likelihoods(allele_depth, marker_depth, ploidy, likelihood_frequencies[allele], settings)
            for allele, allele_depth in enumerate(sample_depths)
        ]
        log_likelihoods.append(sample_log_likelihoods)
        posteriors.append(
            [
                dosage_posterior(sample_log_likelihoods[allele], priors[ploidy][allele])
                for allele in range(len(sample_depths))
            ]
        )
    return log_likelihoods, posteriors


def frequencies_from_posteriors(
    posteriors: Sequence[Sequence[Sequence[float]]], ploidies: Sequence[int]
) -> List[float]:
    """Allele frequencies re-estimated from the posteriors.

    A sample's estimate of an allele's frequency is its posterior mean dosage
    divided by its ploidy; the marker's estimate is the mean over samples.
    Every sample counts equally, whatever its ploidy or read depth, as in
    polyRAD.
    """
    allele_count = len(posteriors[0])
    totals = [0.0] * allele_count
    for sample_posteriors, ploidy in zip(posteriors, ploidies):
        for allele, posterior in enumerate(sample_posteriors):
            totals[allele] += posterior_mean_dosage(posterior) / ploidy
    return [total / len(posteriors) for total in totals]


def call_marker(
    depths: Sequence[Sequence[int]],
    ploidies: Sequence[int],
    settings: ModelSettings = ModelSettings(),
    start_frequencies: Optional[Sequence[float]] = None,
) -> MarkerCall:
    """Estimate dosage posteriors for every sample and allele of one marker.

    depths: read counts, ``[sample][allele]``.
    ploidies: one ploidy per sample.
    start_frequencies: allele frequencies to start from, in place of the
        depth-ratio estimates.

    With the ``uniform`` prior there is one pass, and the allele frequencies
    are used only by the contamination term.  With the ``hwe`` prior the
    passes repeat until the frequencies settle (see ``ModelSettings``).  With
    the ``norm`` prior each pass also re-fits the prior to the posteriors
    (the two steps of an EM iteration), and the passes repeat until neither
    the frequencies nor any prior probability changes by more than the
    tolerance.
    """
    if len(depths) != len(ploidies):
        raise ValueError("depths and ploidies must have one entry per sample")
    frequencies = list(start_frequencies) if start_frequencies is not None else depth_ratio_frequencies(depths)
    if frequencies is None:
        return MarkerCall(callable=False, allele_frequencies=[], log_likelihoods=[], posteriors=[])

    ploidy_groups = sorted(set(ploidies))
    priors = _frequency_priors(ploidy_groups, frequencies, settings)
    if settings.prior == PRIOR_UNIFORM:
        log_likelihoods, posteriors = _marker_pass(depths, ploidies, frequencies, priors, settings)
        return MarkerCall(True, frequencies, log_likelihoods, posteriors, iterations=1, priors=priors)

    normal_parameters: Dict[int, List[Tuple[float, float]]] = {}
    iterations = 0
    while True:
        iterations += 1
        log_likelihoods, posteriors = _marker_pass(depths, ploidies, frequencies, priors, settings)
        new_frequencies = frequencies_from_posteriors(posteriors, ploidies)
        change = sum(abs(new - old) for new, old in zip(new_frequencies, frequencies)) / len(frequencies)
        used_priors = priors
        if settings.prior == PRIOR_NORM:
            priors, normal_parameters = _fitted_normal_priors(
                posteriors, ploidies, new_frequencies, normal_parameters
            )
            largest_prior_change = max(
                abs(new - old)
                for ploidy in ploidy_groups
                for new_prior, old_prior in zip(priors[ploidy], used_priors[ploidy])
                for new, old in zip(new_prior, old_prior)
            )
            change = max(change, largest_prior_change)
        else:
            priors = _frequency_priors(ploidy_groups, new_frequencies, settings)
        frequencies = new_frequencies
        if change <= settings.tolerance or iterations >= settings.max_iterations:
            break
    return MarkerCall(
        True,
        frequencies,
        log_likelihoods,
        posteriors,
        iterations=iterations,
        final_change=change,
        converged=change <= settings.tolerance,
        priors=used_priors,
        normal_prior_parameters=normal_parameters,
    )


# ---------------------------------------------------------------------------
# From marginal posteriors to one genotype per sample
# ---------------------------------------------------------------------------


def most_probable_dosage(posterior: Sequence[float]) -> Optional[int]:
    """The dosage with the highest posterior probability, or None when two tie for it."""
    highest = max(posterior)
    if sum(1 for probability in posterior if probability == highest) > 1:
        return None
    return list(posterior).index(highest)


def _log_probability(probability: float) -> float:
    return math.log(probability) if probability > 0.0 else NEGATIVE_INFINITY


def correct_dosages_to_ploidy(
    dosages: Sequence[Optional[int]],
    posteriors: Sequence[Sequence[float]],
    ploidy: int,
) -> Optional[List[int]]:
    """Adjust the dosages of one sample's alleles until they add up to ``ploidy``.

    ``dosages`` holds each allele's most probable dosage and ``posteriors``
    its dosage posterior.  This is polyRAD's ``CorrectGenos``: while the sum
    is short, add one copy to the allele for which doing so leaves the product
    of the posteriors of all the dosages highest; while it is over, remove one
    copy by the same rule.  One copy at a time is fast but not guaranteed to
    reach the best possible set of dosages; ``best_dosages_summing_to_ploidy``
    is.

    Returns None (no call) when any allele's dosage is None, or when every
    possible adjustment has zero probability.
    """
    if any(dosage is None for dosage in dosages):
        return None
    corrected = [int(dosage) for dosage in dosages]
    while sum(corrected) != ploidy:
        step = 1 if sum(corrected) < ploidy else -1
        best_allele = None
        best_log_probability = NEGATIVE_INFINITY
        for allele, dosage in enumerate(corrected):
            if not 0 <= dosage + step <= ploidy:
                continue
            log_probability = sum(
                _log_probability(posteriors[other][corrected[other] + (step if other == allele else 0)])
                for other in range(len(corrected))
            )
            if log_probability > best_log_probability:
                best_allele, best_log_probability = allele, log_probability
        if best_allele is None:
            return None
        corrected[best_allele] += step
    return corrected


def best_dosages_summing_to_ploidy(posteriors: Sequence[Sequence[float]], ploidy: int) -> Optional[List[int]]:
    """The dosages that add up to ``ploidy`` with the highest product of marginal posteriors.

    It maximises the same quantity as ``correct_dosages_to_ploidy`` but over
    every combination, so its result is never less probable.  The search goes
    through the alleles in order, keeping for each number of copies used so
    far the best assignment that uses that many; the work grows with the
    number of alleles times the square of the ploidy.

    Returns None when no combination has a probability above zero.
    """
    # best_by_copies_used[copies] = (log probability, dosages of the alleles handled so far)
    best_by_copies_used: Dict[int, Tuple[float, List[int]]] = {0: (0.0, [])}
    for posterior in posteriors:
        extended: Dict[int, Tuple[float, List[int]]] = {}
        for copies_used, (log_probability, chosen) in best_by_copies_used.items():
            for dosage in range(ploidy - copies_used + 1):
                if posterior[dosage] <= 0.0:
                    continue
                candidate = log_probability + math.log(posterior[dosage])
                total_copies = copies_used + dosage
                if total_copies not in extended or candidate > extended[total_copies][0]:
                    extended[total_copies] = (candidate, chosen + [dosage])
        best_by_copies_used = extended
    if ploidy not in best_by_copies_used:
        return None
    return best_by_copies_used[ploidy][1]


# ---------------------------------------------------------------------------
# Estimating the overdispersion from the data
# ---------------------------------------------------------------------------

# Values tried by ``estimate_overdispersion``.  Steps widen as the values grow, because
# the beta-binomial changes less and less as it approaches the binomial.
OVERDISPERSION_GRID = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500)
# The estimate uses at most this many sample-by-marker cells, from markers taken at even
# spacing through the table, so that it stays quick and the same table always gives the
# same value.
OVERDISPERSION_MAX_CELLS = 20000
OVERDISPERSION_MIN_MARKERS = 20
# With fewer sample-by-marker cells that have reads, the estimate would be mostly noise.
OVERDISPERSION_MIN_CELLS = 100
OVERDISPERSION_MAX_ROUNDS = 4


@dataclasses.dataclass
class OverdispersionEstimate:
    """Result of ``estimate_overdispersion``.

    ``value`` is the chosen overdispersion, or None when fewer than
    ``OVERDISPERSION_MIN_CELLS`` sample-by-marker cells had reads
    (``cells_with_reads``), which is too little to estimate from.
    ``log_likelihood`` maps every value tried to the log-likelihood of the
    read counts under it, in the last round.  ``at_limit`` is ``"lowest"`` or
    ``"highest"`` when the chosen value is an end of the grid, so the best
    value may lie beyond it.  ``values_by_round`` is the value chosen in each
    round, and ``markers_used`` the number of markers the estimate is based on.
    """

    value: Optional[float]
    log_likelihood: Dict[float, float]
    markers_used: int
    values_by_round: List[float]
    at_limit: Optional[str] = None
    cells_with_reads: int = 0


def log_likelihood_of_counts(
    depths: Sequence[Sequence[int]],
    ploidies: Sequence[int],
    alternate_frequency: float,
    settings: ModelSettings,
) -> float:
    """Log-likelihood of one two-allele marker's read counts, with the genotypes summed out.

    depths: ``[sample][allele]`` read counts of the reference and alternate allele.

    A sample's genotype is unknown, so the probability of its counts is the
    average over the possible dosages of the alternate allele, weighted by
    their Hardy-Weinberg probabilities (with ``settings.selfing_rate``).  The
    samples' logs are added.
    """
    likelihood_frequency = frequency_for_likelihood(alternate_frequency, len(depths))
    priors = {ploidy: hwe_prior(ploidy, alternate_frequency, settings.selfing_rate) for ploidy in set(ploidies)}
    total = 0.0
    for (reference_depth, alternate_depth), ploidy in zip(depths, ploidies):
        marker_depth = reference_depth + alternate_depth
        if marker_depth == 0:
            continue
        log_terms = [
            log_likelihood + math.log(probability)
            for log_likelihood, probability in zip(
                allele_log_likelihoods(alternate_depth, marker_depth, ploidy, likelihood_frequency, settings),
                priors[ploidy],
            )
            if probability > 0.0
        ]
        largest = max(log_terms)
        if largest == NEGATIVE_INFINITY:
            continue
        total += largest + math.log(sum(math.exp(term - largest) for term in log_terms))
    return total


def _evenly_spaced(items: Sequence, limit: int) -> List:
    """At most ``limit`` of the items, taken at even spacing."""
    if len(items) <= limit:
        return list(items)
    return [items[index * len(items) // limit] for index in range(limit)]


def estimate_overdispersion(
    depths_by_marker: Sequence[Sequence[Sequence[int]]],
    ploidies: Sequence[int],
    contamination_rate: float = DEFAULT_CONTAMINATION_RATE,
    grid: Sequence[float] = OVERDISPERSION_GRID,
    start: float = DEFAULT_OVERDISPERSION,
    selfing_rate: float = 0.0,
) -> OverdispersionEstimate:
    """Choose the overdispersion that makes the observed read counts most probable.

    depths_by_marker: ``[marker][sample][allele]`` read counts of two alleles
        per marker (the reference-type and alternate-type rows).

    For every value in ``grid`` the log-likelihood of the counts is added up
    over the markers (``log_likelihood_of_counts``), and the value with the
    highest total is chosen.  The allele frequencies that the likelihood
    needs are estimated by the genotype model under the current value, so
    the search is repeated with the chosen value until the choice no longer
    changes; the frequencies depend little on the overdispersion, and this
    usually takes two rounds.

    This is a maximum-likelihood estimate.  polyRAD's ``TestOverdispersion``
    instead tests the read counts of samples that a rough genotype estimate
    assigns one allele copy.  That works in diploids, but in polyploids the
    rough estimate assigns the wrong dosage often enough that the counts look
    more scattered than they are (see docs/DEVELOPMENT.md).
    """
    samples_with_reads = [
        sum(1 for sample_depths in depths if sum(sample_depths) > 0) for depths in depths_by_marker
    ]
    usable = [(depths, count) for depths, count in zip(depths_by_marker, samples_with_reads) if count > 0]
    marker_limit = max(OVERDISPERSION_MIN_MARKERS, OVERDISPERSION_MAX_CELLS // max(1, len(ploidies)))
    chosen_markers = _evenly_spaced(usable, marker_limit)
    cells_with_reads = sum(count for _, count in chosen_markers)
    if cells_with_reads < OVERDISPERSION_MIN_CELLS:
        return OverdispersionEstimate(None, {}, len(chosen_markers), [], cells_with_reads=cells_with_reads)
    markers = [depths for depths, _ in chosen_markers]

    current = start
    values_by_round: List[float] = []
    log_likelihood: Dict[float, float] = {}
    for _ in range(OVERDISPERSION_MAX_ROUNDS):
        frequency_settings = ModelSettings(
            overdispersion=current, contamination_rate=contamination_rate, selfing_rate=selfing_rate
        )
        alternate_frequencies = [
            call_marker(depths, ploidies, frequency_settings).allele_frequencies[1] for depths in markers
        ]
        log_likelihood = {}
        for value in grid:
            settings = ModelSettings(
                overdispersion=value, contamination_rate=contamination_rate, selfing_rate=selfing_rate
            )
            log_likelihood[value] = sum(
                log_likelihood_of_counts(depths, ploidies, frequency, settings)
                for depths, frequency in zip(markers, alternate_frequencies)
            )
        best = max(grid, key=lambda value: log_likelihood[value])
        values_by_round.append(best)
        if best == current:
            break
        current = best
    chosen = values_by_round[-1]
    at_limit = "lowest" if chosen == min(grid) else "highest" if chosen == max(grid) else None
    return OverdispersionEstimate(chosen, log_likelihood, len(markers), values_by_round, at_limit, cells_with_reads)


# ---------------------------------------------------------------------------
# Input: the MADC count table and the ploidy table
# ---------------------------------------------------------------------------

VCF_TARGET = "target"
VCF_ALL = "all"
VCF_HAPLOTYPE = "haplotype"
VCF_TYPES = (VCF_TARGET, VCF_ALL, VCF_HAPLOTYPE)
VCF_FILE_SUFFIX = {VCF_TARGET: ".targets.vcf", VCF_ALL: ".all.vcf", VCF_HAPLOTYPE: ".haplotypes.vcf"}

READ_TOTAL_COLUMN = "readCountSum"
REFERENCE_CATEGORIES = ("Ref", "RefMatch")
ALTERNATE_CATEGORIES = ("Alt", "AltMatch")
CONCRETE_BASES = frozenset("ACGT")

DEFAULT_PLOIDY = 2
DEFAULT_MIN_DEPTH = 10
DEFAULT_MIN_GENOTYPE_QUALITY = 10
DEFAULT_MIN_MARGINAL_POSTERIOR = 0.9
MAX_GENOTYPE_QUALITY = 99
MAX_PHRED_LIKELIHOOD = 99999
# A record is listed in the QC report when at least this share of its reads was left out.
NOTABLE_EXCLUDED_READ_SHARE = 0.10


class MadcRow(NamedTuple):
    """One allele row of a MADC: its ID, category (Ref, Alt, RefMatch, AltMatch or Other),
    tag sequence and read count in each sample."""

    allele_id: str
    category: str
    sequence: str
    counts: Tuple[int, ...]


@dataclasses.dataclass
class MarkerTask:
    """Everything needed to call one marker.

    ``design_ref_tag``, ``design_alt_tag`` and ``declared_target_position``
    (1-based) come from ``--reference-tags`` and are empty or None when the
    marker has no entry there; the Ref and Alt rows are used instead.
    """

    clone_id: str
    rows: List[MadcRow]
    design_ref_tag: str = ""
    design_alt_tag: str = ""
    declared_target_position: Optional[int] = None


def read_count_table(path: Path) -> Tuple[List[str], Dict[str, List[MadcRow]], Dict[str, object]]:
    """Read a MADC into sample names, rows grouped by marker, and a count-overflow summary.

    Accepts the layouts ``filter-paralogs`` does: a MADC from ``generate``, a
    vendor MADC with or without its starred preamble, and HapApp's
    fixed-allele-ID table.  Sample names are the header row's.  Markers keep
    the order in which they first appear.
    """
    with open_text(path, "rt") as handle:
        table = list(csv.reader(handle))
    header_index = find_madc_header(table)
    header = [cell.strip() for cell in table[header_index]]
    allele_column, clone_column, sequence_column = madc_identity_columns(header, path)
    sample_start = madc_sample_start(header)
    # Processed MADCs can carry a per-row read total between the sequence and the samples.
    if header[sample_start : sample_start + 1] == [READ_TOTAL_COLUMN]:
        sample_start += 1
    sample_names = header[sample_start:]
    if not sample_names or any(not name for name in sample_names):
        raise MadcError(f"MADC has no sample columns, or a sample column without a name: {path}")
    repeated = sorted(name for name, occurrences in Counter(sample_names).items() if occurrences > 1)
    if repeated:
        raise MadcError(f"MADC has repeated sample names ({', '.join(repeated[:5])}): {path}")

    rows_by_marker: Dict[str, List[MadcRow]] = {}
    overflow = CountOverflowCheck()
    for cells in table[header_index + 1 :]:
        if len(cells) <= allele_column or not cells[allele_column].strip():
            continue
        padded = list(cells) + [""] * max(0, len(header) - len(cells))
        clone_id = padded[clone_column].strip()
        counts = madc_read_counts(padded[sample_start : sample_start + len(sample_names)], path)
        overflow.add(clone_id, counts)
        allele_id = padded[allele_column].strip()
        rows_by_marker.setdefault(clone_id, []).append(
            MadcRow(allele_id, allele_type(allele_id), padded[sequence_column].strip().upper(), tuple(counts))
        )
    if not rows_by_marker:
        raise MadcError(f"MADC has no allele rows: {path}")
    warning = overflow.warning(path)
    if warning:
        print(warning, file=sys.stderr)
    return sample_names, rows_by_marker, overflow.summary()


def read_ploidy_table(path: Path) -> Dict[str, int]:
    """Read sample ploidies from a two-column table: sample name, ploidy.

    Columns are separated by a tab or a comma.  A first line whose second
    column is not a whole number is taken as a header and skipped.
    """
    ploidy_by_sample: Dict[str, int] = {}
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = [field.strip() for field in re.split(r"[\t,]", line.rstrip("\r\n"))]
            if not any(fields):
                continue
            if len(fields) < 2:
                raise MadcError(f"{path} line {line_number}: expected a sample name and a ploidy")
            name, value = fields[0], fields[1]
            if not re.fullmatch(r"\d+", value):
                if line_number == 1:
                    continue
                raise MadcError(f"{path} line {line_number}: ploidy {value!r} is not a whole number")
            if int(value) < 1:
                raise MadcError(f"{path} line {line_number}: ploidy must be at least 1")
            if name in ploidy_by_sample:
                raise MadcError(f"{path} lists sample {name!r} more than once")
            ploidy_by_sample[name] = int(value)
    return ploidy_by_sample


# ---------------------------------------------------------------------------
# What each VCF type counts
# ---------------------------------------------------------------------------


class TargetVariant(NamedTuple):
    """The designed variant of a marker: 1-based position on the Ref tag and its two alleles.

    A substitution has single-base alleles.  An insertion or deletion is
    written the VCF way, starting at the base before it (``C`` > ``CAAAA``).
    """

    position: int
    ref_allele: str
    alt_allele: str


def _is_concrete(sequence: str) -> bool:
    return bool(sequence) and set(sequence) <= CONCRETE_BASES


def resolve_target(
    ref_tag: str, alt_tag: str, declared_position: Optional[int], gap_band: int
) -> Tuple[Optional[TargetVariant], Optional[str]]:
    """Find the designed variant from the marker's Ref and Alt tags.

    Returns ``(target, None)``, or ``(None, reason)`` when the target cannot
    be identified with certainty.  It is certain when the two tags differ by
    exactly one substitution, insertion or deletion, or when they differ by
    more but ``declared_position`` (``VariantPositionInTagSequence`` of the
    reference tags) names one of the substitutions.  The target is never
    guessed among several differences.
    """
    if not alt_tag:
        return None, "no_alt_row"
    differences = [
        event for event in _designed_events(ref_tag, alt_tag, gap_band) if event[0] in ("snp", "ins", "del")
    ]
    if len(differences) != 1 and declared_position is not None:
        differences = [event for event in differences if event[0] == "snp" and event[1] == declared_position]
    if len(differences) != 1:
        return None, "ambiguous_design"
    _, position, ref_allele, alt_allele = differences[0]
    if not _is_concrete(ref_allele) or not _is_concrete(alt_allele):
        return None, "ambiguous_design"
    return TargetVariant(position, ref_allele, alt_allele), None


class RowReading:
    """Which base a row's tag shows at each position of the marker's Ref tag.

    Built from the row's events against the Ref tag (``madcats.variant_events``,
    the same reading that fills the MADC's Variants column).  A position is
    not shown when the row deleted it, when an insertion pushed it out of the
    fixed-length window, or when the row's tag is not the Ref tag's length.
    """

    def __init__(self, ref_tag: str, sequence: str, events) -> None:
        self.ref_tag = ref_tag
        self.comparable = len(sequence) == len(ref_tag)
        self.substitutions: Dict[int, str] = {}
        self.positions_not_shown: Set[int] = set()
        for kind, position, ref_allele, alt_allele in events:
            if kind == "snp":
                self.substitutions[position] = alt_allele
            elif kind == "del":
                # Written from the base before the deleted ones: that base is still shown.
                # A deletion at the very start of the tag is written from the base after it
                # instead; the two cannot always be told apart there, so the whole span is
                # treated as not shown.
                first_missing = position if position == 1 else position + 1
                self.positions_not_shown.update(range(first_missing, position + len(ref_allele)))
            elif kind == "clip":
                self.positions_not_shown.update(range(position, len(ref_tag) + 1))

    def base_at(self, position: int) -> Optional[str]:
        """The base at a 1-based Ref-tag position, or None when the row does not show one
        (including an N or other ambiguous base in the read)."""
        if not self.comparable or position in self.positions_not_shown:
            return None
        base = self.substitutions.get(position, self.ref_tag[position - 1])
        return base if base in CONCRETE_BASES else None


def vcf_reference_sequence(tag: str) -> str:
    """A tag as a VCF REF allele.

    VCF allows only A, C, G, T and N in REF, and says to replace any other
    ambiguity code by the alphabetically first base it stands for.
    """
    return "".join(
        base if base in CONCRETE_BASES or base == "N" else min(IUPAC_BASES.get(base, "N")) for base in tag
    )


def _summed_counts(counts_by_row: Sequence[Sequence[int]], row_indices: Sequence[int], sample_count: int) -> List[int]:
    """Per-sample read totals over the given rows."""
    totals = [0] * sample_count
    for row_index in row_indices:
        for sample, count in enumerate(counts_by_row[row_index]):
            totals[sample] += count
    return totals


def target_allele_depths(
    rows: Sequence[MadcRow], counts_by_row: Sequence[Sequence[int]], sample_count: int
) -> Tuple[List[List[int]], List[int]]:
    """Reads for and against the designed variant, and the reads that show neither.

    Returns ``[sample][allele]`` depths of the reference allele (Ref and
    RefMatch rows) and the alternate allele (Alt and AltMatch rows), and the
    per-sample depth of every other row.
    """
    reference_rows = [index for index, row in enumerate(rows) if row.category in REFERENCE_CATEGORIES]
    alternate_rows = [index for index, row in enumerate(rows) if row.category in ALTERNATE_CATEGORIES]
    other_rows = [
        index for index, row in enumerate(rows) if row.category not in REFERENCE_CATEGORIES + ALTERNATE_CATEGORIES
    ]
    reference_depth = _summed_counts(counts_by_row, reference_rows, sample_count)
    alternate_depth = _summed_counts(counts_by_row, alternate_rows, sample_count)
    depths = [list(pair) for pair in zip(reference_depth, alternate_depth)]
    return depths, _summed_counts(counts_by_row, other_rows, sample_count)


def _without_negative_cells(rows: Sequence[MadcRow], sample_count: int) -> Tuple[List[List[int]], int]:
    """The rows' counts with every sample that has a negative count at the marker set to zero.

    A negative read count (signed 8-bit overflow in some vendor MADCs) makes
    the sample's depth at the marker unknown, so the sample is treated as
    having no reads there.  Returns the counts and the number of samples
    cleared.
    """
    counts_by_row = [list(row.counts) for row in rows]
    samples_with_negative = {
        sample for counts in counts_by_row for sample, count in enumerate(counts) if count < 0
    }
    for counts in counts_by_row:
        for sample in samples_with_negative:
            counts[sample] = 0
    return counts_by_row, len(samples_with_negative)


# ---------------------------------------------------------------------------
# VCF records
# ---------------------------------------------------------------------------


@dataclasses.dataclass(frozen=True)
class CallOptions:
    """Settings of one ``madcats call`` run that the per-marker work needs."""

    settings: ModelSettings
    vcf_types: Tuple[str, ...]
    min_depth: int = DEFAULT_MIN_DEPTH
    min_genotype_quality: int = DEFAULT_MIN_GENOTYPE_QUALITY
    min_marginal_posterior: float = DEFAULT_MIN_MARGINAL_POSTERIOR
    max_other_fraction: Optional[float] = None
    gap_band: int = GAP_BAND_DEFAULT
    overdispersion_estimated: bool = False


def _decimal(value: float, places: int = 4) -> str:
    """A number with at most ``places`` decimals and no trailing zeros."""
    text = f"{value:.{places}f}".rstrip("0").rstrip(".")
    return text if text not in ("", "-0") else "0"


def genotype_text(dosages: Optional[Sequence[int]], ploidy: int) -> str:
    """The VCF GT value for allele dosages (``[1, 0, 2]`` gives ``0/2/2``), or all-missing.

    Alleles are joined by ``/`` because the copies are not phased.
    """
    if dosages is None:
        return "/".join(["."] * ploidy)
    return "/".join(str(allele) for allele, dosage in enumerate(dosages) for _ in range(dosage))


def genotype_quality(posterior: Sequence[float], called_dosage: int) -> int:
    """Phred-scaled probability that the called dosage is wrong, capped at 99."""
    chance_wrong = sum(probability for dosage, probability in enumerate(posterior) if dosage != called_dosage)
    if chance_wrong <= 0.0:
        return MAX_GENOTYPE_QUALITY
    return min(MAX_GENOTYPE_QUALITY, int(round(-10.0 * math.log10(chance_wrong))))


def phred_scaled_likelihoods(log_likelihoods: Sequence[float]) -> List[int]:
    """VCF PL values: each dosage's likelihood relative to the most likely one, Phred-scaled."""
    best = max(log_likelihoods)
    scaled = []
    for log_likelihood in log_likelihoods:
        if log_likelihood == NEGATIVE_INFINITY:
            scaled.append(MAX_PHRED_LIKELIHOOD)
        else:
            scaled.append(min(MAX_PHRED_LIKELIHOOD, int(round(-10.0 * (log_likelihood - best) / math.log(10.0)))))
    return scaled


def _two_allele_genotype(
    posteriors: Sequence[Sequence[float]], log_likelihoods: Sequence[Sequence[float]], ploidy: int
) -> Tuple[Optional[List[int]], Optional[float], List[str]]:
    """Call a sample at a record with a reference and one alternate allele.

    Returns the allele dosages (None on a tie), the genotype quality, and the
    GQ, DS, GP and PL values.  The alternate allele's posterior is the
    genotype posterior, because the reference dosage is the ploidy minus it.
    """
    alternate_posterior = posteriors[1]
    alternate_dosage = most_probable_dosage(alternate_posterior)
    quality = None if alternate_dosage is None else genotype_quality(alternate_posterior, alternate_dosage)
    values = [
        "." if quality is None else str(quality),
        _decimal(posterior_mean_dosage(alternate_posterior), 3),
        ",".join(_decimal(probability) for probability in alternate_posterior),
        ",".join(str(value) for value in phred_scaled_likelihoods(log_likelihoods[1])),
    ]
    dosages = None if alternate_dosage is None else [ploidy - alternate_dosage, alternate_dosage]
    return dosages, quality, values


def _multiallelic_genotype(
    posteriors: Sequence[Sequence[float]], ploidy: int, tallies: Counter
) -> Tuple[Optional[List[int]], Optional[float], List[str]]:
    """Call a sample at a record with any number of alleles from the marginal posteriors.

    Returns the allele dosages (None when there is no call), the lowest
    marginal posterior among the called dosages (MPC), and the DS and MPC
    values.  ``tallies`` counts how often the dosages had to be corrected to
    the ploidy, and how often an exact search would have chosen differently.
    """
    most_probable = [most_probable_dosage(posterior) for posterior in posteriors]
    dosages = correct_dosages_to_ploidy(most_probable, posteriors, ploidy)
    lowest_posterior = None
    if dosages is not None:
        if sum(dosage for dosage in most_probable if dosage is not None) != ploidy:
            tallies["genotypes_corrected_to_ploidy"] += 1
            if best_dosages_summing_to_ploidy(posteriors, ploidy) != dosages:
                tallies["corrections_where_exact_search_differs"] += 1
        lowest_posterior = min(posterior[dosage] for posterior, dosage in zip(posteriors, dosages))
    mean_dosages = ",".join(_decimal(posterior_mean_dosage(posterior), 3) for posterior in posteriors[1:])
    values = [mean_dosages or ".", "." if lowest_posterior is None else _decimal(lowest_posterior)]
    return dosages, lowest_posterior, values


def format_record(
    marker: str,
    position: int,
    record_id: str,
    alleles: Sequence[str],
    depths: Sequence[Sequence[int]],
    excluded_depths: Optional[Sequence[int]],
    ploidies: Sequence[int],
    options: CallOptions,
    two_allele_model: bool,
    is_target: bool = False,
) -> Tuple[str, Counter]:
    """Call every sample at one record and format its VCF line.

    alleles: the reference allele, then the alternate alleles.
    depths: ``[sample][allele]`` read counts, in the order of ``alleles``.
    excluded_depths: reads per sample left out of ``depths`` (written as XD),
        or None for a record that uses every read.
    two_allele_model: write GQ, GP and PL from the alternate allele's
        posterior (requires exactly two alleles); otherwise write MPC.

    A genotype is written as missing, with its AD and DP kept, when the
    sample has fewer than ``min_depth`` reads, when no single dosage is most
    probable, when its confidence is below the threshold, or when too large a
    share of its reads was excluded.  Returns the line and the tallies for
    the QC report.
    """
    result = call_marker(depths, ploidies, options.settings)
    tallies: Counter = Counter(records=1)
    if result.callable and not result.converged:
        tallies["records_not_converged"] += 1
    used_reads = sum(sum(sample_depths) for sample_depths in depths)
    excluded_reads = sum(excluded_depths) if excluded_depths is not None else 0
    tallies["reads_used"] += used_reads
    tallies["reads_excluded"] += excluded_reads
    if excluded_reads and excluded_reads >= NOTABLE_EXCLUDED_READ_SHARE * (used_reads + excluded_reads):
        tallies["records_with_many_excluded_reads"] += 1

    value_count = 4 if two_allele_model else 2
    sample_columns = []
    samples_called = 0
    for sample, (sample_depths, ploidy) in enumerate(zip(depths, ploidies)):
        depth = sum(sample_depths)
        excluded = excluded_depths[sample] if excluded_depths is not None else 0
        dosages = None
        if not result.callable or depth == 0:
            missing_reason = "no_reads"
            values = ["."] * value_count
        else:
            if two_allele_model:
                dosages, confidence, values = _two_allele_genotype(
                    result.posteriors[sample], result.log_likelihoods[sample], ploidy
                )
                confidence_threshold = options.min_genotype_quality
            else:
                dosages, confidence, values = _multiallelic_genotype(result.posteriors[sample], ploidy, tallies)
                confidence_threshold = options.min_marginal_posterior
            if depth < options.min_depth:
                missing_reason = "low_depth"
            elif dosages is None:
                missing_reason = "no_single_best_dosage"
            elif confidence < confidence_threshold:
                missing_reason = "low_confidence"
            elif options.max_other_fraction is not None and excluded > options.max_other_fraction * (depth + excluded):
                missing_reason = "excluded_read_share"
            else:
                missing_reason = None
        tallies["genotypes"] += 1
        if missing_reason is None:
            samples_called += 1
            tallies["genotypes_called"] += 1
        else:
            dosages = None
            tallies[f"genotypes_missing_{missing_reason}"] += 1
        column = [genotype_text(dosages, ploidy), ",".join(str(count) for count in sample_depths), str(depth)]
        if excluded_depths is not None:
            column.append(str(excluded))
        sample_columns.append(":".join(column + values))

    info = ["TARGET"] if is_target else []
    info.append(f"NS={samples_called}")
    if result.callable and len(alleles) > 1:
        info.append("AF=" + ",".join(_decimal(frequency) for frequency in result.allele_frequencies[1:]))
    format_keys = ["GT", "AD", "DP"] + (["XD"] if excluded_depths is not None else [])
    format_keys += ["GQ", "DS", "GP", "PL"] if two_allele_model else ["DS", "MPC"]
    line = "\t".join(
        [marker, str(position), record_id, alleles[0], ",".join(alleles[1:]) or ".", ".", ".", ";".join(info),
         ":".join(format_keys)]
        + sample_columns
    )
    return line, tallies


@dataclasses.dataclass
class MarkerOutput:
    """What ``process_marker`` produced for one marker.

    ``lines`` holds the VCF lines per VCF type, ``tallies`` the QC counts
    keyed by ``(vcf type or "markers", name)``, and ``haplotype_alleles`` one
    ``(allele index, MADC allele IDs, sequence, total reads)`` entry per
    allele of the haplotype record.  ``tag_length`` is 0 for a marker without
    a Ref tag, which gets no record.
    """

    clone_id: str
    tag_length: int = 0
    lines: Dict[str, List[str]] = dataclasses.field(default_factory=dict)
    tallies: Counter = dataclasses.field(default_factory=Counter)
    haplotype_alleles: List[Tuple[int, str, str, int]] = dataclasses.field(default_factory=list)


def _first_sequence(rows: Sequence[MadcRow], category: str) -> str:
    return next((row.sequence for row in rows if row.category == category and row.sequence), "")


def process_marker(task: MarkerTask, ploidies: Sequence[int], options: CallOptions) -> MarkerOutput:
    """Call one marker and format its records for every requested VCF type."""
    output = MarkerOutput(task.clone_id, lines={vcf_type: [] for vcf_type in options.vcf_types})
    sample_count = len(ploidies)
    rows = [row for row in task.rows if row.sequence]
    ref_tag = task.design_ref_tag or _first_sequence(rows, "Ref")
    alt_tag = task.design_alt_tag or _first_sequence(rows, "Alt")
    if not ref_tag:
        output.tallies[("markers", "skipped_no_ref_row")] += 1
        return output
    output.tag_length = len(ref_tag)
    counts_by_row, samples_cleared = _without_negative_cells(rows, sample_count)
    output.tallies[("markers", "sample_cells_cleared_for_negative_counts")] += samples_cleared

    def add(vcf_type: str, line: str, tallies: Counter) -> None:
        output.lines[vcf_type].append(line)
        for name, count in tallies.items():
            output.tallies[(vcf_type, name)] += count

    site_types = [vcf_type for vcf_type in (VCF_TARGET, VCF_ALL) if vcf_type in options.vcf_types]
    target = None
    if site_types:
        target, skip_reason = resolve_target(ref_tag, alt_tag, task.declared_target_position, options.gap_band)
        if target is None:
            for vcf_type in site_types:
                output.tallies[(vcf_type, f"targets_skipped_{skip_reason}")] += 1
        else:
            depths, other_depths = target_allele_depths(rows, counts_by_row, sample_count)
            line, tallies = format_record(
                task.clone_id, target.position, task.clone_id, [target.ref_allele, target.alt_allele],
                depths, other_depths, ploidies, options, two_allele_model=True, is_target=True,
            )
            for vcf_type in site_types:
                add(vcf_type, line, tallies)

    if VCF_ALL in options.vcf_types:
        for position, line, tallies in _off_target_records(
            task.clone_id, rows, counts_by_row, ref_tag, alt_tag, target, ploidies, options
        ):
            add(VCF_ALL, line, tallies)
        # The target was added first; put the marker's records in position order.
        output.lines[VCF_ALL].sort(key=lambda record: int(record.split("\t", 2)[1]))

    if VCF_HAPLOTYPE in options.vcf_types:
        line, tallies, output.haplotype_alleles = _haplotype_record(
            task.clone_id, rows, counts_by_row, ref_tag, ploidies, options
        )
        add(VCF_HAPLOTYPE, line, tallies)
    return output


def _off_target_records(
    clone_id: str,
    rows: Sequence[MadcRow],
    counts_by_row: Sequence[Sequence[int]],
    ref_tag: str,
    alt_tag: str,
    target: Optional[TargetVariant],
    ploidies: Sequence[int],
    options: CallOptions,
):
    """Yield ``(position, line, tallies)`` for each substitution site away from the target.

    A site is a Ref-tag position where at least one row with reads shows a
    base other than the Ref tag's.  Reads are grouped by the base their row
    shows there; reads of rows that do not show the position are excluded.
    Insertions and deletions are not written as records.
    """
    sample_count = len(ploidies)
    # variant_events needs a Locus but reads only its two tags.
    locus = Locus(
        index=0, marker=clone_id, left_primer="", right_primer="", ref_base="", alt_base="", variant_pos=-1,
        ref_tag=ref_tag, alt_tag=alt_tag or ref_tag,
    )
    readings = [
        RowReading(ref_tag, row.sequence, variant_events(locus, row.sequence, options.gap_band)) for row in rows
    ]
    rows_with_reads = [index for index, counts in enumerate(counts_by_row) if any(counts)]
    target_is_substitution = target is not None and len(target.ref_allele) == len(target.alt_allele) == 1
    site_positions = sorted(
        {
            position
            for index in rows_with_reads
            if readings[index].comparable
            for position, base in readings[index].substitutions.items()
            if base in CONCRETE_BASES
            and ref_tag[position - 1] in CONCRETE_BASES
            and not (target_is_substitution and position == target.position)
        }
    )
    for position in site_positions:
        rows_by_base: Dict[str, List[int]] = {}
        rows_not_showing = []
        for index in rows_with_reads:
            base = readings[index].base_at(position)
            if base is None:
                rows_not_showing.append(index)
            else:
                rows_by_base.setdefault(base, []).append(index)
        depth_by_base = {
            base: _summed_counts(counts_by_row, indices, sample_count) for base, indices in rows_by_base.items()
        }
        reference_base = ref_tag[position - 1]
        depth_by_base.setdefault(reference_base, [0] * sample_count)
        # Alternate bases with the most reads first, then alphabetically.
        alternate_bases = sorted(
            (base for base in depth_by_base if base != reference_base),
            key=lambda base: (-sum(depth_by_base[base]), base),
        )
        alleles = [reference_base] + alternate_bases
        line, tallies = format_record(
            clone_id, position, f"{clone_id}_{position}", alleles,
            [[depth_by_base[base][sample] for base in alleles] for sample in range(sample_count)],
            _summed_counts(counts_by_row, rows_not_showing, sample_count),
            ploidies, options, two_allele_model=len(alleles) == 2,
        )
        yield position, line, tallies


def _haplotype_record(
    clone_id: str,
    rows: Sequence[MadcRow],
    counts_by_row: Sequence[Sequence[int]],
    ref_tag: str,
    ploidies: Sequence[int],
    options: CallOptions,
) -> Tuple[str, Counter, List[Tuple[int, str, str, int]]]:
    """The marker's haplotype record: the Ref tag, then every other tag sequence with reads.

    Rows with the same sequence are added together.  Alternate alleles keep
    the order of their first row in the MADC.
    """
    sample_count = len(ploidies)
    reference_sequence = vcf_reference_sequence(ref_tag)
    rows_by_sequence: Dict[str, List[int]] = {reference_sequence: []}
    for index, row in enumerate(rows):
        sequence = reference_sequence if row.sequence == ref_tag else row.sequence
        rows_by_sequence.setdefault(sequence, []).append(index)
    depth_by_sequence = {
        sequence: _summed_counts(counts_by_row, indices, sample_count) for sequence, indices in rows_by_sequence.items()
    }
    alleles = [reference_sequence] + [
        sequence for sequence in rows_by_sequence if sequence != reference_sequence and any(depth_by_sequence[sequence])
    ]
    line, tallies = format_record(
        clone_id, 1, clone_id, alleles,
        [[depth_by_sequence[sequence][sample] for sequence in alleles] for sample in range(sample_count)],
        None, ploidies, options, two_allele_model=False,
    )
    allele_table = [
        (
            allele_index,
            ";".join(rows[index].allele_id for index in rows_by_sequence[sequence]) or ".",
            sequence,
            sum(depth_by_sequence[sequence]),
        )
        for allele_index, sequence in enumerate(alleles)
    ]
    return line, tallies, allele_table


def model_description(settings: ModelSettings, overdispersion_estimated: bool = False) -> str:
    """One sentence stating the genotype model and its inheritance assumption, for the VCF header."""
    source = ", estimated from the read counts" if overdispersion_estimated else ""
    if settings.prior == PRIOR_HWE:
        mating = (
            f"a selfing rate of {settings.selfing_rate:g}" if settings.selfing_rate > 0 else "random mating"
        )
        prior = (
            "Hardy-Weinberg prior with one allele frequency per record estimated over all samples together; "
            f"it assumes {mating} and that each chromosome copy is drawn independently "
            "(polysomic inheritance in polyploids)"
        )
    elif settings.prior == PRIOR_NORM:
        prior = (
            "normal-shaped prior over dosages (updog's norm; Gerard and Ferrao 2020), with a mean and standard "
            f"deviation fitted per record and ploidy; a ploidy with fewer than {NORMAL_PRIOR_MIN_SAMPLES} samples "
            "uses Hardy-Weinberg proportions"
        )
    else:
        prior = "uniform prior over dosages"
    return (
        f"beta-binomial read-count likelihood (overdispersion {settings.overdispersion:g}{source}, contamination rate "
        f"{settings.contamination_rate:g}) adapted from polyRAD; {prior}"
    )


def vcf_header(
    vcf_type: str, sample_names: Sequence[str], contigs: Sequence[Tuple[str, int]], options: CallOptions
) -> List[str]:
    """Header lines of one VCF.  It carries no date, so the same input gives the same file."""
    multiallelic_scope = "this record's alleles" if vcf_type == VCF_HAPLOTYPE else "records with more than two alleles"
    other_fraction = "none" if options.max_other_fraction is None else format(options.max_other_fraction, "g")
    lines = [
        "##fileformat=VCFv4.3",
        f"##source=MADCATS call v{__version__}",
        f"##madcatsVcfType={vcf_type}",
        "##madcatsCoordinates=marker-local: CHROM is a marker (CloneID) and POS is a 1-based position on its "
        "Ref tag. These are not genome coordinates.",
        f"##madcatsModel={model_description(options.settings, options.overdispersion_estimated)}",
        f"##madcatsFilters=min_dp={options.min_depth}, min_gq={options.min_genotype_quality}, "
        f"min_marginal_posterior={options.min_marginal_posterior:g}, "
        f"max_other_fraction={other_fraction}"
        "; a genotype that fails is written as missing with its AD and DP kept",
    ]
    lines += [f"##contig=<ID={marker},length={length}>" for marker, length in contigs]
    if vcf_type != VCF_HAPLOTYPE:
        lines.append('##INFO=<ID=TARGET,Number=0,Type=Flag,Description="The variant the marker was designed for">')
    lines += [
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Samples with a called genotype">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated frequency of each alternate allele">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype; allele copies are unphased">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Reads of each allele">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Reads used for the call (the sum of AD)">',
    ]
    if vcf_type != VCF_HAPLOTYPE:
        lines += [
            '##FORMAT=<ID=XD,Number=1,Type=Integer,Description="Reads left out of AD and DP: Other rows at a target, '
            'rows that do not show the position elsewhere">',
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Phred-scaled posterior probability that GT is wrong '
            '(two-allele records)">',
        ]
    lines.append(
        '##FORMAT=<ID=DS,Number=A,Type=Float,Description="Posterior mean number of copies of each alternate allele">'
    )
    if vcf_type != VCF_HAPLOTYPE:
        lines += [
            '##FORMAT=<ID=GP,Number=G,Type=Float,Description="Posterior probability of each genotype '
            '(two-allele records)">',
            '##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihood of each genotype, before the '
            'prior (two-allele records)">',
        ]
    if vcf_type != VCF_TARGET:
        lines.append(
            '##FORMAT=<ID=MPC,Number=1,Type=Float,Description="Lowest marginal posterior probability among the called '
            f"dosages of {multiallelic_scope}. Each allele's dosage is estimated separately, so this is an upper bound "
            'on the probability that the whole genotype is right, not an estimate of it">'
        )
    columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
    lines.append("\t".join(columns + list(sample_names)))
    return lines


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------


def _fraction(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed <= 1:
        raise argparse.ArgumentTypeError("must be between 0 and 1")
    return parsed


def _positive_float(value: str) -> float:
    parsed = float(value)
    if parsed <= 0:
        raise argparse.ArgumentTypeError("must be greater than 0")
    return parsed


OVERDISPERSION_AUTO = "auto"


def _overdispersion_argument(value: str):
    return OVERDISPERSION_AUTO if value.strip().lower() == OVERDISPERSION_AUTO else _positive_float(value)


def choose_overdispersion(
    tasks: Sequence[MarkerTask], ploidies: Sequence[int], contamination_rate: float, selfing_rate: float = 0.0
) -> Tuple[float, Dict[str, object]]:
    """Estimate the overdispersion for ``--overdispersion auto``.

    Returns the value to use and a report for the QC JSON.  The estimate is
    made from each marker's reference-type and alternate-type read counts.
    It falls back to the default when there is too little to estimate from, or
    when the best value is the lowest one tried: that is what a table of
    inbred lines or another sample set far from Hardy-Weinberg proportions
    produces, whatever its real read-count scatter.
    """
    sample_count = len(ploidies)
    depths_by_marker = []
    for task in tasks:
        rows = [row for row in task.rows if row.sequence]
        if not any(row.category in ALTERNATE_CATEGORIES for row in rows):
            continue
        counts_by_row, _ = _without_negative_cells(rows, sample_count)
        depths_by_marker.append(target_allele_depths(rows, counts_by_row, sample_count)[0])
    estimate = estimate_overdispersion(depths_by_marker, ploidies, contamination_rate, selfing_rate=selfing_rate)
    report: Dict[str, object] = {
        "estimate": estimate.value,
        "markers_used": estimate.markers_used,
        "cells_with_reads": estimate.cells_with_reads,
        "values_by_round": estimate.values_by_round,
        "at_limit": estimate.at_limit,
        "log_likelihood_by_value": {
            format(value, "g"): round(total, 3) for value, total in estimate.log_likelihood.items()
        },
    }
    if estimate.value is None:
        print(
            f"note: too few read counts to estimate the overdispersion ({estimate.cells_with_reads:,} "
            f"sample-by-marker cells with reads; {OVERDISPERSION_MIN_CELLS} needed). Using "
            f"{DEFAULT_OVERDISPERSION:g}.",
            file=sys.stderr,
        )
        return DEFAULT_OVERDISPERSION, report
    if estimate.at_limit == "lowest":
        print(
            f"warning: the overdispersion estimate is the lowest value tried ({estimate.value:g}). The samples "
            "probably do not follow the genotype proportions the estimate assumes (Hardy-Weinberg, with "
            f"--selfing-rate if given), as with inbred lines or one family. Using {DEFAULT_OVERDISPERSION:g}; "
            "set --overdispersion to a number to choose it yourself.",
            file=sys.stderr,
        )
        return DEFAULT_OVERDISPERSION, report
    limit_note = ""
    if estimate.at_limit == "highest":
        limit_note = "; this is the highest value tried, so the counts scatter no more than random sampling of reads"
    print(
        f"Estimated overdispersion: {estimate.value:g} (from {estimate.markers_used:,} markers{limit_note})",
        file=sys.stderr,
    )
    return float(estimate.value), report


def _contamination_rate(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 0 and below 1")
    return parsed


def run_call(args: argparse.Namespace) -> int:
    """Call genotypes from a MADC and write the requested VCFs (see the module docstring)."""
    madc = Path(args.madc)
    if not madc.is_file():
        raise MadcError(f"MADC not found: {madc}")
    started = time.time()
    vcf_types = tuple(vcf_type for vcf_type in VCF_TYPES if vcf_type in args.vcf)
    if args.selfing_rate > 0 and args.prior != PRIOR_HWE:
        raise MadcError("--selfing-rate applies to --prior hwe only")

    sample_names, rows_by_marker, count_overflow = read_count_table(madc)
    ploidy_by_sample = read_ploidy_table(Path(args.ploidy_file)) if args.ploidy_file else {}
    ploidies = [ploidy_by_sample.get(name, args.ploidy) for name in sample_names]
    samples_not_in_madc = sorted(set(ploidy_by_sample) - set(sample_names))
    if samples_not_in_madc:
        print(
            f"note: {len(samples_not_in_madc):,} samples in {args.ploidy_file} are not in the MADC "
            f"(first: {samples_not_in_madc[0]})",
            file=sys.stderr,
        )

    if 0 < args.selfing_rate < 1 and any(ploidy % 2 for ploidy in ploidies):
        print(
            "note: --selfing-rate is not applied to samples of odd ploidy; they keep random-mating proportions",
            file=sys.stderr,
        )

    small_ploidy_groups = sorted(
        ploidy for ploidy, count in Counter(ploidies).items() if count < NORMAL_PRIOR_MIN_SAMPLES
    )
    if args.prior == PRIOR_NORM and small_ploidy_groups:
        print(
            f"note: fewer than {NORMAL_PRIOR_MIN_SAMPLES} samples have ploidy "
            f"{', '.join(str(ploidy) for ploidy in small_ploidy_groups)}; --prior norm uses Hardy-Weinberg "
            "proportions for them",
            file=sys.stderr,
        )

    designs = {}
    if args.reference_tags:
        designs = {locus.marker: locus for locus in load_reference_tags(Path(args.reference_tags))}
    tasks = []
    for clone_id, rows in rows_by_marker.items():
        design = designs.get(clone_id)
        if design is None:
            tasks.append(MarkerTask(clone_id, rows))
        else:
            tasks.append(MarkerTask(clone_id, rows, design.ref_tag, design.alt_tag, design.variant_pos + 1))

    overdispersion_report: Dict[str, object] = {"requested": args.overdispersion}
    if args.overdispersion == OVERDISPERSION_AUTO:
        overdispersion, estimate_report = choose_overdispersion(
            tasks, ploidies, args.contamination_rate, args.selfing_rate
        )
        overdispersion_report.update(estimate_report)
    else:
        overdispersion = args.overdispersion
    overdispersion_report["used"] = overdispersion
    # False when --overdispersion auto fell back to the default.
    overdispersion_estimated = (
        overdispersion_report.get("estimate") is not None and overdispersion_report.get("at_limit") != "lowest"
    )

    options = CallOptions(
        settings=ModelSettings(
            overdispersion=overdispersion,
            contamination_rate=args.contamination_rate,
            prior=args.prior,
            selfing_rate=args.selfing_rate,
            tolerance=args.tol,
            max_iterations=args.max_iterations,
        ),
        vcf_types=vcf_types,
        min_depth=args.min_dp,
        min_genotype_quality=args.min_gq,
        min_marginal_posterior=args.min_marginal_posterior,
        max_other_fraction=args.max_other_fraction,
        gap_band=args.max_indel_length,
        overdispersion_estimated=overdispersion_estimated,
    )

    worker = functools.partial(process_marker, ploidies=ploidies, options=options)
    if args.jobs > 1 and len(tasks) > 1:
        chunk_size = max(1, len(tasks) // (args.jobs * 8))
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            outputs = list(executor.map(worker, tasks, chunksize=chunk_size))
    else:
        outputs = [worker(task) for task in tasks]

    contigs = [(output.clone_id, output.tag_length) for output in outputs if output.tag_length]
    tallies: Counter = Counter()
    for output in outputs:
        tallies.update(output.tallies)

    prefix = str(args.output_prefix)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    written = {}
    for vcf_type in vcf_types:
        vcf_path = Path(prefix + VCF_FILE_SUFFIX[vcf_type])
        with vcf_path.open("w") as handle:
            for line in vcf_header(vcf_type, sample_names, contigs, options):
                handle.write(line + "\n")
            for output in outputs:
                for line in output.lines[vcf_type]:
                    handle.write(line + "\n")
        written[vcf_type] = vcf_path
    alleles_path = None
    if VCF_HAPLOTYPE in vcf_types:
        # Ties each allele number of the haplotype VCF to the MADC rows it was counted from.
        alleles_path = Path(prefix + ".haplotypes.alleles.csv")
        with alleles_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["CloneID", "AlleleIndex", "AlleleIDs", "AlleleSequence", "TotalReads"])
            for output in outputs:
                for allele_index, allele_ids, sequence, total_reads in output.haplotype_alleles:
                    writer.writerow([output.clone_id, allele_index, allele_ids, sequence, total_reads])

    per_type = {
        vcf_type: {name: count for (group, name), count in sorted(tallies.items()) if group == vcf_type}
        for vcf_type in vcf_types
    }
    qc_path = Path(args.qc_output) if args.qc_output else Path(prefix + ".qc.json")
    qc = {
        "pipeline": "MADCATS",
        "command": "call",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "madc": str(madc.resolve()),
            "reference_tags": str(Path(args.reference_tags).resolve()) if args.reference_tags else None,
            "ploidy_file": str(Path(args.ploidy_file).resolve()) if args.ploidy_file else None,
        },
        "parameters": {
            "vcf": list(vcf_types),
            "ploidy": args.ploidy,
            "prior": args.prior,
            "selfing_rate": args.selfing_rate,
            "overdispersion": args.overdispersion,
            "contamination_rate": args.contamination_rate,
            "tol": args.tol,
            "max_iterations": args.max_iterations,
            "min_dp": args.min_dp,
            "min_gq": args.min_gq,
            "min_marginal_posterior": args.min_marginal_posterior,
            "max_other_fraction": args.max_other_fraction,
            "max_indel_length": args.max_indel_length,
            "jobs": args.jobs,
        },
        "model": model_description(options.settings, options.overdispersion_estimated),
        "overdispersion": overdispersion_report,
        "coordinates": "marker-local",
        "samples": {
            "count": len(sample_names),
            "by_ploidy": {str(ploidy): count for ploidy, count in sorted(Counter(ploidies).items())},
            "ploidy_file_samples_not_in_madc": len(samples_not_in_madc),
            "ploidies_using_hwe_under_norm_prior": small_ploidy_groups if args.prior == PRIOR_NORM else [],
        },
        "markers": {
            "in_madc": len(tasks),
            "with_reference_tags": sum(1 for task in tasks if task.design_ref_tag),
            **{name: count for (group, name), count in sorted(tallies.items()) if group == "markers"},
        },
        "count_overflow": count_overflow,
        "vcfs": {
            vcf_type: {"path": str(written[vcf_type]), **per_type[vcf_type]} for vcf_type in vcf_types
        },
        "outputs": {"haplotype_alleles": str(alleles_path) if alleles_path else None},
        "elapsed_seconds": round(time.time() - started, 1),
    }
    _write_json(qc_path, qc)

    for vcf_type in vcf_types:
        counts = per_type[vcf_type]
        genotypes = counts.get("genotypes", 0)
        called = counts.get("genotypes_called", 0)
        share = called / genotypes if genotypes else 0.0
        print(
            f"{written[vcf_type]}: {counts.get('records', 0):,} records, {called:,} of {genotypes:,} genotypes "
            f"called ({share:.1%})",
            file=sys.stderr,
        )
        skipped = sum(count for name, count in counts.items() if name.startswith("targets_skipped_"))
        if skipped:
            print(
                f"  {skipped:,} markers have no target record because the designed variant could not be "
                "identified; see the QC JSON, and give --reference-tags if you have them",
                file=sys.stderr,
            )
        if counts.get("records_not_converged"):
            print(
                f"  {counts['records_not_converged']:,} records reached --max-iterations before their allele "
                "frequencies settled",
                file=sys.stderr,
            )
    print(f"Wrote QC to {qc_path}", file=sys.stderr)
    return 0


def add_call_command(subparsers) -> None:
    calling = subparsers.add_parser(
        "call",
        help="call genotypes from a MADC and write target, all-variant or haplotype VCFs",
        description=(
            "Estimate each sample's genotype from the read counts in a MADC and write VCFs. --vcf chooses the "
            "record types: target (the designed variant of each marker; the default), all (targets plus other "
            "base substitutions in the tags) and haplotype (whole tag sequences as alleles). Each type goes to "
            "its own file. Ploidy can differ between samples. The model is polyRAD's: a beta-binomial "
            "likelihood of the read counts with a uniform or Hardy-Weinberg prior. Coordinates are "
            "marker-local: CHROM is the marker and POS a position on its Ref tag. Reads MADCs from generate, "
            "vendor DArT MADCs and HapApp fixed-allele-ID MADCs."
        ),
    )
    calling.add_argument(
        "--madc",
        required=True,
        help=(
            "MADC to call from: from generate or filter-paralogs, a vendor DArT MADC (with or without its "
            "starred preamble), or a HapApp fixed-allele-ID MADC"
        ),
    )
    calling.add_argument(
        "--output-prefix",
        required=True,
        help="path prefix of the outputs: PREFIX.targets.vcf, PREFIX.all.vcf, PREFIX.haplotypes.vcf and PREFIX.qc.json",
    )
    calling.add_argument(
        "--vcf",
        nargs="+",
        choices=VCF_TYPES,
        default=[VCF_TARGET],
        help="VCF types to write, one file each; give one or more (default: target)",
    )
    calling.add_argument(
        "--reference-tags",
        help=(
            "reference-tag CSV of the panel (from design or infer-design). Its Ref and Alt tags and "
            "VariantPositionInTagSequence identify each marker's target; without it the MADC's Ref and Alt "
            "rows are used, and a marker whose two tags differ in more than one place gets no target record"
        ),
    )
    calling.add_argument(
        "--ploidy", type=positive_int, default=DEFAULT_PLOIDY,
        help=f"ploidy of every sample not listed in --ploidy-file (default: {DEFAULT_PLOIDY})",
    )
    calling.add_argument(
        "--ploidy-file",
        help="two-column table of sample name and ploidy, tab- or comma-separated, with or without a header line",
    )
    calling.add_argument(
        "--prior",
        choices=PRIOR_CHOICES,
        default=PRIOR_HWE,
        help=(
            "genotype prior: uniform (every dosage equally probable), hwe (Hardy-Weinberg proportions for "
            "allele frequencies estimated from all samples; assumes random mating and polysomic inheritance) "
            "or norm (a normal shape over the dosages fitted per record and ploidy, as in updog; does not "
            f"assume random mating) (default: {PRIOR_HWE})"
        ),
    )
    calling.add_argument(
        "--selfing-rate", type=_fraction, default=0.0,
        help=(
            "with --prior hwe, the share of individuals that come from self-fertilisation, from 0 (random "
            "mating) to 1 (fully inbred lines); it moves prior probability from the intermediate dosages to "
            "the two homozygous ones (default: 0)"
        ),
    )
    calling.add_argument(
        "--overdispersion", type=_overdispersion_argument, default=OVERDISPERSION_AUTO,
        help=(
            "beta-binomial overdispersion parameter, as in polyRAD: lower values allow read fractions to stray "
            "further from the expected fraction. 'auto' estimates it from the read counts by maximum "
            f"likelihood, and uses {DEFAULT_OVERDISPERSION:g} when it cannot; or give a number "
            f"(default: {OVERDISPERSION_AUTO})"
        ),
    )
    calling.add_argument(
        "--contamination-rate", type=_contamination_rate, default=DEFAULT_CONTAMINATION_RATE,
        help=(
            "expected share of a sample's reads that come from other samples "
            f"(default: {DEFAULT_CONTAMINATION_RATE:g})"
        ),
    )
    calling.add_argument(
        "--tol", type=_positive_float, default=DEFAULT_TOLERANCE,
        help=(
            "with --prior hwe or norm, stop iterating a record once its allele frequencies (and, with norm, "
            f"its prior probabilities) change by no more than this (default: {DEFAULT_TOLERANCE:g})"
        ),
    )
    calling.add_argument(
        "--max-iterations", type=positive_int, default=DEFAULT_MAX_ITERATIONS,
        help=f"with --prior hwe or norm, the most iterations a record may take (default: {DEFAULT_MAX_ITERATIONS})",
    )
    calling.add_argument(
        "--min-dp", type=nonnegative_int, default=DEFAULT_MIN_DEPTH,
        help=(
            "write a genotype as missing when the sample has fewer reads than this in the record "
            f"(default: {DEFAULT_MIN_DEPTH})"
        ),
    )
    calling.add_argument(
        "--min-gq", type=nonnegative_int, default=DEFAULT_MIN_GENOTYPE_QUALITY,
        help=(
            "two-allele records: write a genotype as missing when its genotype quality is below this "
            f"(default: {DEFAULT_MIN_GENOTYPE_QUALITY})"
        ),
    )
    calling.add_argument(
        "--min-marginal-posterior", type=_fraction, default=DEFAULT_MIN_MARGINAL_POSTERIOR,
        help=(
            "haplotype records and other records with more than two alleles: write a genotype as missing when "
            "MPC, the lowest marginal posterior among its called dosages, is below this "
            f"(default: {DEFAULT_MIN_MARGINAL_POSTERIOR:g})"
        ),
    )
    calling.add_argument(
        "--max-other-fraction", type=_fraction, default=None,
        help=(
            "target and all records: write a genotype as missing when more than this share of the sample's "
            "reads was left out of the call (XD) (default: no limit)"
        ),
    )
    calling.add_argument(
        "--max-indel-length", type=nonnegative_int, default=GAP_BAND_DEFAULT,
        help=(
            "widest insertion or deletion considered when a row is compared with the Ref tag; use the value "
            f"given to generate (default: {GAP_BAND_DEFAULT})"
        ),
    )
    calling.add_argument("--qc-output", help="QC JSON (default: PREFIX.qc.json)")
    calling.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(),
        help="worker processes; markers are called independently (default: up to 4)",
    )
    calling.set_defaults(func=run_call)
