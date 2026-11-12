"""Genotype calling from allele read counts (the engine behind ``madcats call``).

The statistical model is adapted from polyRAD (Clark, Lipka and Sacks 2019,
G3 9:663-673; https://github.com/lvclark/polyRAD, GPL >= 2).  At one marker,
each allele is considered on its own against all the other alleles of that
marker together, which turns a multiallelic marker into one two-outcome
problem per allele:

1. Likelihood.  For a sample of ploidy P, an allele can be present in
   0, 1, ..., P copies (its *dosage*).  The probability of the observed number
   of reads of the allele, out of all the sample's reads at the marker, is
   beta-binomial.  The expected read fraction for a dosage is dosage / P,
   moved slightly toward the population allele frequency to allow for
   contamination between samples.
2. Prior.  Either every dosage is equally probable (``uniform``), or dosages
   follow Hardy-Weinberg proportions for the allele's frequency (``hwe``):
   each of the P chromosome copies carries the allele independently.  That is
   the expectation for a diploid or an autopolyploid in a randomly mating
   population.  It is not the expectation under disomic inheritance in an
   allopolyploid, or in a family.
3. Posterior.  Prior times likelihood, normalised over the dosages.
4. Iteration (``hwe`` only).  The allele frequencies are re-estimated from
   the posterior mean dosages and steps 1-3 repeat until the frequencies
   stop changing.  Markers do not share any parameter, so each marker
   converges on its own (polyRAD stops all markers together, on their mean
   change).

The posteriors are *marginal*: one distribution per allele, not a joint
distribution over whole genotypes.  The most probable dosages of a marker's
alleles therefore need not add up to the ploidy, and
``correct_dosages_to_ploidy`` (polyRAD's method) or
``best_dosages_summing_to_ploidy`` (an exact search) reconciles them.  At a
marker with exactly two alleles the two marginals mirror each other, and the
posterior of either allele is the full genotype posterior.

Likelihoods are kept as logarithms and only normalised posteriors are
exponentiated, so a dosage that the reads make very improbable keeps a
usable value at any read depth.

This module uses only the Python standard library.
"""

from __future__ import annotations

import dataclasses
import math
from typing import Dict, List, Optional, Sequence, Tuple

PRIOR_UNIFORM = "uniform"
PRIOR_HWE = "hwe"
PRIOR_CHOICES = (PRIOR_UNIFORM, PRIOR_HWE)

# polyRAD's defaults: AddGenotypeLikelihood(overdispersion = 9), readDArTag(contamRate = 0.001)
# and IterateHWE(tol = 1e-5).  polyRAD has no iteration limit.
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
    prior: ``uniform`` or ``hwe``.
    tolerance: the ``hwe`` iteration stops once the mean absolute change in
        the marker's allele frequencies is at or below this.
    max_iterations: upper limit on ``hwe`` iterations.
    """

    overdispersion: float = DEFAULT_OVERDISPERSION
    contamination_rate: float = DEFAULT_CONTAMINATION_RATE
    prior: str = PRIOR_HWE
    tolerance: float = DEFAULT_TOLERANCE
    max_iterations: int = DEFAULT_MAX_ITERATIONS

    def __post_init__(self) -> None:
        if self.overdispersion <= 0:
            raise ValueError("overdispersion must be positive")
        if not 0 <= self.contamination_rate < 1:
            raise ValueError("contamination_rate must be at least 0 and below 1")
        if self.prior not in PRIOR_CHOICES:
            raise ValueError(f"prior must be one of {', '.join(PRIOR_CHOICES)}")
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


def hwe_prior(ploidy: int, allele_frequency: float) -> List[float]:
    """Dosage probabilities under Hardy-Weinberg proportions.

    Each of the ``ploidy`` chromosome copies carries the allele independently
    with probability ``allele_frequency``, so the dosage is binomial.
    """
    other_frequency = 1.0 - allele_frequency
    return [
        math.comb(ploidy, dosage) * allele_frequency**dosage * other_frequency ** (ploidy - dosage)
        for dosage in range(ploidy + 1)
    ]


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


def _marker_pass(
    depths: Sequence[Sequence[int]],
    ploidies: Sequence[int],
    allele_frequencies: Sequence[float],
    settings: ModelSettings,
) -> Tuple[List[List[List[float]]], List[List[List[float]]]]:
    """Likelihoods and posteriors of every sample and allele for given allele frequencies."""
    sample_count = len(depths)
    likelihood_frequencies = [frequency_for_likelihood(frequency, sample_count) for frequency in allele_frequencies]
    # Samples of the same ploidy share their priors.
    priors: Dict[int, List[List[float]]] = {}
    for ploidy in set(ploidies):
        if settings.prior == PRIOR_HWE:
            priors[ploidy] = [hwe_prior(ploidy, frequency) for frequency in allele_frequencies]
        else:
            priors[ploidy] = [uniform_prior(ploidy)] * len(allele_frequencies)

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


def frequencies_from_posteriors(posteriors: Sequence[Sequence[Sequence[float]]], ploidies: Sequence[int]) -> List[float]:
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
    passes repeat until the frequencies settle (see ``ModelSettings``).
    """
    if len(depths) != len(ploidies):
        raise ValueError("depths and ploidies must have one entry per sample")
    frequencies = list(start_frequencies) if start_frequencies is not None else depth_ratio_frequencies(depths)
    if frequencies is None:
        return MarkerCall(callable=False, allele_frequencies=[], log_likelihoods=[], posteriors=[])

    if settings.prior == PRIOR_UNIFORM:
        log_likelihoods, posteriors = _marker_pass(depths, ploidies, frequencies, settings)
        return MarkerCall(True, frequencies, log_likelihoods, posteriors, iterations=1)

    iterations = 0
    while True:
        iterations += 1
        log_likelihoods, posteriors = _marker_pass(depths, ploidies, frequencies, settings)
        new_frequencies = frequencies_from_posteriors(posteriors, ploidies)
        change = sum(abs(new - old) for new, old in zip(new_frequencies, frequencies)) / len(frequencies)
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
