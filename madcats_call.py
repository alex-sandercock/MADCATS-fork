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

The statistical model is adapted from two R packages.  Each function below
that follows one of them names the function it follows, and says where it
differs.

polyRAD 2.0.1 (GPL >= 2; https://github.com/lvclark/polyRAD).  Clark LV,
    Lipka AE, Sacks EJ (2019) polyRAD: genotype calling with uncertainty from
    sequencing data in polyploids and diploids.  G3 9(3):663-673.
    doi:10.1534/g3.118.200913.
    From it: the beta-binomial likelihood of each allele against the rest of
    the record with its contamination term, the Hardy-Weinberg prior and its
    selfing adjustment (after de Silva HN, Hall AJ, Rikkerink E, McNeilage
    MA, Fraser LG (2005) Heredity 95:327-334, doi:10.1038/sj.hdy.6800728),
    the iteration of the allele frequencies, the marginal posteriors, the
    correction of dosages to the ploidy, and the form of a prior with one
    allele frequency per sample.
updog 2.1.7 (GPL >= 3; https://github.com/dcgerard/updog).  Gerard D, Ferrão
    LFV, Garcia AAF, Stephens M (2018) Genotyping polyploids from messy
    sequencing data.  Genetics 210(3):789-807.
    doi:10.1534/genetics.118.301468.  Gerard D, Ferrão LFV (2020) Priors for
    genotyping polyploids.  Bioinformatics 36(6):1795-1800.
    doi:10.1093/bioinformatics/btz852.
    From it: the allelic bias with its prior and starting values, the form
    of the term for reads that show the wrong allele, the ``norm`` prior with
    its limits, and the free-weights prior (updog's ``flex``).

MADCATS's own, with no published counterpart: the estimate of the
overdispersion and of one background-read rate for the whole table
(``estimate_overdispersion``, ``estimate_read_model``) and its use of the
free-weights prior, the order of bias before background reads, the search for
the allelic bias (``fit_allele_bias``), the rules that use a fitted bias only
where it gains a margin and does not end at a limit of the search
(``BIAS_AUTO_MIN_GAIN``, ``BIAS_REJECTED_AT_LIMIT``) and the separate read
model of records called without a bias (``choose_read_model``), the fit of the
``norm`` prior by Newton's method and its fit per ploidy, the neighbour-based
``structure`` prior (``find_neighbours``), the exact search for dosages that
sum to the ploidy as a check on polyRAD's rule, and the VCF fields (GQ, GP,
PL, MPC).
The numerical agreement with both packages is checked by tests/test_call.py
against values stored in tests/fixtures/.

In a record, each allele is considered on its own against all the other
alleles together, which turns a multiallelic record into one two-outcome
problem per allele:

1. Likelihood.  For a sample of ploidy P, an allele can be present in
   0, 1, ..., P copies (its *dosage*).  The probability of the observed number
   of reads of the allele, out of all the sample's reads in the record, is
   beta-binomial.  The expected read fraction for a dosage is dosage / P,
   changed by an allelic bias if one is fitted (updog's; Gerard et al. 2018,
   Genetics 210:789-807) and by the share of reads that show the wrong
   allele (see ``expected_read_fraction``).
2. Prior.  Either every dosage is equally probable (``uniform``), or dosages
   follow Hardy-Weinberg proportions for the allele's frequency (``hwe``):
   each of the P chromosome copies carries the allele independently.  That is
   the expectation for a diploid or an autopolyploid in a randomly mating
   population.  It is not the expectation under disomic inheritance in an
   allopolyploid, or in a family.  The ``norm`` prior makes no such
   assumption: the dosage probabilities follow a normal shape whose mean and
   spread are fitted to the samples of each ploidy.  The ``structure`` prior
   keeps Hardy-Weinberg proportions but gives every sample its own allele
   frequency, taken from the samples most like it across the table's
   markers (``find_neighbours``), for samples from several groups, clines or
   families.
3. Posterior.  Prior times likelihood, normalised over the dosages.
4. Iteration (``hwe``, ``norm`` and ``structure``).  The allele frequencies,
   and for ``norm`` the prior's two parameters, are re-estimated from the
   posteriors and steps 1-3 repeat until they stop changing.  Records do not
   share any parameter, so each converges on its own (polyRAD stops all loci
   together, on their mean change).  With ``structure`` the lists of like
   samples are made once for the whole table, and each record then iterates
   on its own with them.

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

Reading map: the model (``expected_read_fraction``, ``dosage_shapes``,
``log_likelihoods_for_counts``, ``hwe_prior``, ``normal_prior``, ``fit_normal_prior``, ``dosage_posterior``,
``call_record``, ``fit_allele_bias``); from marginal posteriors to one
genotype (``correct_dosages_to_ploidy``, ``best_dosages_summing_to_ploidy``);
the estimate of the overdispersion and the background-read rate
(``estimate_overdispersion``, ``estimate_read_model``); the samples most
like each sample, for the ``structure`` prior (``find_neighbours``); input
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
import heapq
import math
import operator
import re
import sys
import time
from array import array
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
    proper_fraction,
    variant_events,
)

PRIOR_UNIFORM = "uniform"
PRIOR_HWE = "hwe"
PRIOR_NORM = "norm"
PRIOR_STRUCTURE = "structure"
PRIOR_CHOICES = (PRIOR_UNIFORM, PRIOR_HWE, PRIOR_NORM, PRIOR_STRUCTURE)
# A prior with one free probability per dosage (updog's "flex").  It is used to estimate
# the overdispersion and the background-read rate, and is not offered for calling: with
# as many parameters as dosages it follows the reads too closely to steady a call.
PRIOR_FREE = "free"

# polyRAD's defaults: AddGenotypeLikelihood(overdispersion = 9), readDArTag(contamRate = 0.001)
# and IterateHWE(tol = 1e-5).  polyRAD has no iteration limit.  ModelSettings starts from
# these two, so that a bare ModelSettings() is polyRAD's model; the command sets its own.
POLYRAD_OVERDISPERSION = 9.0
POLYRAD_CONTAMINATION_RATE = 0.001
# The command estimates the overdispersion and the background-read rate unless they are
# given, and uses these when it cannot.  The overdispersion stays at polyRAD's default
# until more real tables have been looked at; the two alfalfa tables the background rate
# was first estimated on gave 0.001 and 0.003.
FALLBACK_OVERDISPERSION = POLYRAD_OVERDISPERSION
FALLBACK_BACKGROUND_RATE = 0.002
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
    background_rate: probability that a read shows one allele although it
        came from another (see ``expected_read_fraction``).
    prior: ``uniform``, ``hwe``, ``norm`` (see ``normal_prior``) or
        ``structure`` (see ``neighbours``).
    selfing_rate: with the ``hwe`` prior, the share of individuals that come
        from self-fertilisation, from 0 (random mating) to 1 (fully inbred).
    fit_allele_bias: estimate an allelic bias for every record with two
        alleles (see ``fit_allele_bias``).  Without it both alleles are taken
        to be read equally well.
    min_bias_gain: with ``fit_allele_bias``, what a record's fitted bias must
        add to the log-likelihood of its read counts, after the prior on the
        bias, for the record to be called with it; a record that gains less
        is called without a bias.  0 uses every fitted bias.
    reject_bias_at_limit: with ``fit_allele_bias``, do not use a fitted bias
        that ends at a limit of the search (``BIAS_LOG_LIMIT``); the record
        is called without a bias.
    read_model_without_bias: with ``fit_allele_bias``, the overdispersion and
        the background rate for a record that is called without a bias: one
        whose fitted bias is not used, and one with more than two alleles,
        which gets none.  ``overdispersion`` and ``background_rate`` are then
        for the records called with their bias.  None uses those for every
        record.
    tolerance: the passes of a record stop once the change they watch is at
        or below this: the mean absolute change in the record's allele
        frequencies (``hwe``), also the largest change in any prior
        probability (``norm``), or the mean absolute change in the samples'
        own frequencies (``structure`` with neighbours).
    max_iterations: upper limit on the passes.  The ``structure`` prior has
        two stages (``hwe`` passes, then passes with the neighbours), and
        each has this tolerance and this limit.
    neighbours: with the ``structure`` prior, for each sample the indices of
        the samples most like it (from ``find_neighbours``).  A sample's
        prior is then Hardy-Weinberg at the allele frequency among those
        samples, not among all samples.  Without neighbours the ``structure``
        prior is the ``hwe`` prior.
    """

    overdispersion: float = POLYRAD_OVERDISPERSION
    contamination_rate: float = POLYRAD_CONTAMINATION_RATE
    background_rate: float = 0.0
    prior: str = PRIOR_HWE
    selfing_rate: float = 0.0
    fit_allele_bias: bool = False
    min_bias_gain: float = 0.0
    reject_bias_at_limit: bool = False
    read_model_without_bias: Optional[Tuple[float, float]] = None
    tolerance: float = DEFAULT_TOLERANCE
    max_iterations: int = DEFAULT_MAX_ITERATIONS
    neighbours: Optional[Tuple[Tuple[int, ...], ...]] = None

    def __post_init__(self) -> None:
        if self.overdispersion <= 0:
            raise ValueError("overdispersion must be positive")
        if not 0 <= self.contamination_rate < 1:
            raise ValueError("contamination_rate must be at least 0 and below 1")
        if not 0 <= self.background_rate < 0.5:
            raise ValueError("background_rate must be at least 0 and below 0.5")
        if self.prior not in PRIOR_CHOICES + (PRIOR_FREE,):
            raise ValueError(f"prior must be one of {', '.join(PRIOR_CHOICES)}")
        if not 0 <= self.selfing_rate <= 1:
            raise ValueError("selfing_rate must be between 0 and 1")
        if self.selfing_rate > 0 and self.prior != PRIOR_HWE:
            raise ValueError("selfing_rate applies to the hwe prior only")
        if self.tolerance <= 0:
            raise ValueError("tolerance must be positive")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be at least 1")
        if self.min_bias_gain < 0:
            raise ValueError("min_bias_gain must be at least 0")
        if self.fit_allele_bias and self.prior == PRIOR_STRUCTURE:
            # The bias search scores a bias by the prior shared within a ploidy; it has not
            # been worked out or tested for priors that differ from sample to sample.
            raise ValueError("an allelic bias cannot be fitted with the structure prior")

    def without_bias(self) -> "ModelSettings":
        """These settings for a record that is called without an allelic bias: those of a run that fits none."""
        overdispersion, background_rate = self.read_model_without_bias or (self.overdispersion, self.background_rate)
        return dataclasses.replace(
            self, fit_allele_bias=False, min_bias_gain=0.0, reject_bias_at_limit=False, read_model_without_bias=None,
            overdispersion=overdispersion, background_rate=background_rate,
        )


@dataclasses.dataclass
class RecordCall:
    """The model's result for one marker.

    ``log_likelihoods`` and ``posteriors`` are indexed
    ``[sample][allele][dosage]``, with dosage running from 0 to the sample's
    ploidy.  ``allele_frequencies`` are the final estimates; the posteriors
    were computed from the estimates of the pass before.  ``final_change`` is
    the last value of the change the passes stop on (see
    ``ModelSettings.tolerance``).

    ``priors`` holds the prior the posteriors were computed with, as
    ``[group][allele][dosage]``, and ``prior_groups`` each sample's group:
    its ploidy, since samples of a ploidy share a prior, or under the
    ``structure`` prior with neighbours the sample's own index.  With the
    ``norm`` prior,
    ``normal_prior_parameters`` holds the last fitted
    ``(mean, standard deviation)`` as ``[ploidy][allele]``, for the ploidy
    groups large enough to be fitted.  ``allele_biases`` holds the bias each
    allele was called with (see ``expected_read_fraction``), 1 without one;
    ``bias_fit`` describes the fit when a bias was fitted, whether or not the
    record was called with it.

    With the ``structure`` prior and neighbours, ``sample_frequencies`` holds
    the allele frequencies each sample's own prior was computed at, as
    ``[sample][allele]``.

    ``has_reads`` is False when no sample has a read at the marker: there is
    then no allele frequency to start from, and the other fields are empty.
    """

    has_reads: bool
    allele_frequencies: List[float]
    log_likelihoods: List[List[List[float]]]
    posteriors: List[List[List[float]]]
    iterations: int = 0
    final_change: Optional[float] = None
    converged: bool = True
    priors: Dict[int, List[List[float]]] = dataclasses.field(default_factory=dict)
    prior_groups: List[int] = dataclasses.field(default_factory=list)
    normal_prior_parameters: Dict[int, List[Tuple[float, float]]] = dataclasses.field(default_factory=dict)
    allele_biases: List[float] = dataclasses.field(default_factory=list)
    bias_fit: Optional["BiasFit"] = None
    sample_frequencies: List[List[float]] = dataclasses.field(default_factory=list)


# ---------------------------------------------------------------------------
# Likelihood
# ---------------------------------------------------------------------------


def log_binomial_coefficient(total: int, chosen: int) -> float:
    """Natural log of the number of ways to choose ``chosen`` items out of ``total``."""
    return math.lgamma(total + 1) - math.lgamma(chosen + 1) - math.lgamma(total - chosen + 1)


def minimum_allele_frequency(sample_count: int) -> float:
    """The smallest allele frequency the likelihood is given.

    An estimated frequency of exactly 0 or 1 would make the contamination
    term vanish for that allele.  polyRAD (``AddGenotypeLikelihood``) uses
    ``1 / nTaxa / max(ploidies)``, where the ploidies are the locus ploidies
    it is given as ``possiblePloidies``; with ``list(2)``, as in the stored
    comparisons, that is the 1 / (2 * samples) used here, whatever the
    samples' own ploidy.
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


def expected_read_fraction(
    dosage: int,
    ploidy: int,
    allele_frequency: float,
    contamination_rate: float,
    bias: float = 1.0,
    background_rate: float = 0.0,
) -> float:
    """Probability that one read of a sample shows the allele.

    A read comes from the sample itself with probability
    1 - contamination_rate, and then shows the allele in proportion to its
    dosage.  Otherwise it comes from another sample and shows the allele at
    its population frequency (polyRAD's contamination).

    ``background_rate`` is the probability that a read shows the allele
    although it came from another, and the reverse.  It turns the fraction
    ``f`` into ``f * (1 - background_rate) + (1 - f) * background_rate``,
    which is how updog models sequencing error.  Here it stands for every
    source of such reads that is left after the MADC was built: sequencing
    and amplification errors, index hopping and contamination alike.  Unlike
    polyRAD's contamination it does not depend on the allele's frequency, so
    it also accounts for stray reads of an allele no other sample carries.

    ``bias`` is updog's allelic bias (Gerard et al. 2018, Genetics
    210:789-807): how readily a copy of any other allele yields a read,
    relative to a copy of this allele.  At 1 the two are read alike; above 1
    this allele is under-represented among the reads, and the fraction
    ``f`` becomes ``f / (bias * (1 - f) + f)``.

    The bias is applied to the sample's own copies and the background reads
    are added afterwards.  updog applies its sequencing error first and the
    bias second, and fits an error rate for every marker.  With one
    background rate for the whole table, that order would let the bias of a
    marker stand in for its own rate of stray reads: samples that do not
    carry an allele would pull the bias around, and at a marker with few
    carriers they would outweigh them.  In this order such samples say
    nothing about the bias.  The two orders differ by about the background
    rate in the expected fraction, and not at all when it is 0.
    """
    fraction = (dosage / ploidy) * (1.0 - contamination_rate) + contamination_rate * allele_frequency
    if bias != 1.0:
        fraction = fraction / (bias * (1.0 - fraction) + fraction)
    if background_rate:
        fraction = fraction * (1.0 - background_rate) + (1.0 - fraction) * background_rate
    return fraction


class DosageShape(NamedTuple):
    """The beta-binomial of one dosage: its two shape parameters, and the part of the
    log-probability that depends only on them.  ``fraction`` is the expected read fraction."""

    fraction: float
    shape_allele: float
    shape_other: float
    log_normaliser: float


def dosage_shapes(
    ploidy: int, allele_frequency: float, settings: ModelSettings, bias: float = 1.0
) -> List[DosageShape]:
    """The beta-binomial of every dosage (0 to ``ploidy``) of one allele.

    This and ``log_likelihoods_for_counts`` are the likelihood of polyRAD's
    ``AddGenotypeLikelihood``, kept as logarithms.  The shape parameters are ``fraction * overdispersion`` and
    ``(1 - fraction) * overdispersion``, where ``fraction`` is the expected
    read fraction of the dosage.  They are the same for every sample of the
    ploidy, so they are worked out once and not per sample.
    ``allele_frequency`` must already be inside the open interval (see
    ``frequency_for_likelihood``).
    """
    shapes = []
    for dosage in range(ploidy + 1):
        fraction = expected_read_fraction(
            dosage, ploidy, allele_frequency, settings.contamination_rate, bias, settings.background_rate
        )
        if fraction <= 0.0 or fraction >= 1.0:
            shapes.append(DosageShape(fraction, 0.0, 0.0, 0.0))
            continue
        shape_allele = fraction * settings.overdispersion
        shape_other = (1.0 - fraction) * settings.overdispersion
        shapes.append(
            DosageShape(fraction, shape_allele, shape_other, -math.lgamma(shape_allele) - math.lgamma(shape_other))
        )
    return shapes


def count_constant(allele_depth: int, record_depth: int, overdispersion: float) -> float:
    """The part of the beta-binomial log-probability that is the same for every dosage.

    The two shape parameters add up to the overdispersion whatever the
    dosage, so these terms depend on the read counts alone.
    """
    return (
        log_binomial_coefficient(record_depth, allele_depth)
        + math.lgamma(overdispersion)
        - math.lgamma(record_depth + overdispersion)
    )


def log_likelihoods_for_counts(
    allele_depth: int, record_depth: int, constant: float, shapes: Sequence[DosageShape]
) -> List[float]:
    """Log-likelihood of each dosage for one sample's read counts.

    ``constant`` is ``count_constant`` for these counts and ``shapes`` is
    ``dosage_shapes`` for the sample's ploidy.  With no reads at the marker
    every dosage has likelihood 1 (log 0.0).
    """
    if record_depth == 0:
        return [0.0] * len(shapes)
    other_depth = record_depth - allele_depth
    log_likelihoods = []
    for shape in shapes:
        if shape.fraction <= 0.0 or shape.fraction >= 1.0:
            # Only without contamination: a dosage of 0 cannot produce a read of the
            # allele, and a dosage equal to the ploidy cannot produce any other read.
            impossible = allele_depth > 0 if shape.fraction <= 0.0 else other_depth > 0
            log_likelihoods.append(NEGATIVE_INFINITY if impossible else 0.0)
        else:
            log_likelihoods.append(
                constant
                + math.lgamma(allele_depth + shape.shape_allele)
                + math.lgamma(other_depth + shape.shape_other)
                + shape.log_normaliser
            )
    return log_likelihoods


def allele_log_likelihoods(
    allele_depth: int,
    record_depth: int,
    ploidy: int,
    allele_frequency: float,
    settings: ModelSettings,
    bias: float = 1.0,
) -> List[float]:
    """Log-likelihood of each dosage (0 to ``ploidy``) of one allele in one sample.

    ``allele_depth`` reads out of ``record_depth`` showed the allele.  The
    likelihood is the beta-binomial probability of that count (see
    ``dosage_shapes``).  This is the likelihood written for one sample; the
    calling code uses the two functions it is made of, so that samples of a
    ploidy share their work.
    """
    if record_depth == 0:
        return [0.0] * (ploidy + 1)
    return log_likelihoods_for_counts(
        allele_depth,
        record_depth,
        count_constant(allele_depth, record_depth, settings.overdispersion),
        dosage_shapes(ploidy, allele_frequency, settings, bias),
    )


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
    The values are those of polyRAD's ``.gameteProb``, from the closed form.
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
    dosage is the sum of two gametes of the same parent.  ``S`` is polyRAD's
    ``.selfmat``.
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
FITTED_PRIOR_MIN_SAMPLES = 10


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


# Range of the amount added to the curvature when a Newton step of the normal-shaped prior's
# fit overshoots.  The curvature itself is at most a few hundred (the variance of a squared
# dosage), so the largest amount leaves a step too short to change the fit.
NEWTON_MIN_DAMPING = 1e-12
NEWTON_MAX_DAMPING = 1e8


def _normal_prior_newton_fit(shares: Sequence[float], ploidy: int, mean: float, standard_deviation: float):
    """Fit the normal-shaped prior without limits on its parameters, by Newton's method.

    The quantity maximised is updog's (``obj_for_weighted_lnorm``); updog
    maximises it with L-BFGS-B, and this method is MADCATS's own.
    The prior is written as ``log probability = a * dosage + b * dosage**2 -
    constant``, with ``a = mean / sd**2`` and ``b = -1 / (2 * sd**2)``.  In
    ``a`` and ``b`` the weighted log prior has a single peak, reached where
    the prior's mean dosage and mean squared dosage equal those of
    ``shares``.  Returns ``(a, b)``, or None when the iteration does not
    settle (the shares sit on one or two dosages and the peak is at infinity).

    A Newton step that does not improve the fit is shortened by adding to the
    diagonal of the curvature (the Levenberg-Marquardt method): the more is
    added, the shorter the step and the nearer its direction to the
    gradient's.  A start whose prior is almost entirely on one dosage has
    next to no curvature, and its plain Newton step is too long by many
    orders of magnitude; halving such a step does not rescue it.  When no
    amount added improves the fit, to within rounding, the peak has been
    reached.
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
        damping = 0.0
        while damping <= NEWTON_MAX_DAMPING:
            damped_first, damped_second = variance_first + damping, variance_second + damping
            determinant = damped_first * damped_second - covariance**2
            if determinant > 1e-300:
                step_linear = (damped_second * gradient_first - covariance * gradient_second) / determinant
                step_quadratic = (damped_first * gradient_second - covariance * gradient_first) / determinant
                candidate = prior_and_objective(linear + step_linear, quadratic + step_quadratic)
                if candidate[1] >= objective:
                    break
            damping = max(NEWTON_MIN_DAMPING, 10.0 * damping)
        else:
            return linear, quadratic
        linear, quadratic = linear + step_linear, quadratic + step_quadratic
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


def log_probability(probability: float) -> float:
    """Natural log of a probability, with minus infinity for a probability of zero."""
    return math.log(probability) if probability > 0.0 else NEGATIVE_INFINITY


def log_probabilities(probabilities: Sequence[float]) -> List[float]:
    """``log_probability`` of each of several probabilities."""
    return [math.log(probability) if probability > 0.0 else NEGATIVE_INFINITY for probability in probabilities]


def log_sum_of_exponentials(log_values: Sequence[float]) -> float:
    """Log of the sum of numbers given by their logs; minus infinity when all of them are zero."""
    largest = max(log_values)
    if largest == NEGATIVE_INFINITY:
        return NEGATIVE_INFINITY
    return largest + math.log(sum(math.exp(log_value - largest) for log_value in log_values))


def dosage_posterior(
    log_likelihoods: Sequence[float], prior: Sequence[float], log_prior: Optional[Sequence[float]] = None
) -> List[float]:
    """Posterior probability of each dosage: prior times likelihood, normalised.

    When the prior gives zero probability to every dosage the reads allow,
    the product is zero everywhere and cannot be normalised; the prior is
    returned unchanged, as polyRAD does.  ``log_prior`` can be given when the
    same prior is used for many samples, to take its logs once.
    """
    if log_prior is None:
        log_prior = log_probabilities(prior)
    log_products = [log_likelihood + log_value for log_likelihood, log_value in zip(log_likelihoods, log_prior)]
    largest = max(log_products)
    if largest == NEGATIVE_INFINITY:
        return list(prior)
    # Subtract the largest term before exponentiating so that it becomes exp(0) = 1 and
    # the sum cannot underflow to zero.
    scaled = [math.exp(log_product - largest) for log_product in log_products]
    total = sum(scaled)
    return [value / total for value in scaled]


def posterior_mean_dosage(posterior: Sequence[float]) -> float:
    """Expected number of copies of the allele (polyRAD's ``GetWeightedMeanGenotypes``, per sample)."""
    return sum(dosage * probability for dosage, probability in enumerate(posterior))


# ---------------------------------------------------------------------------
# One marker
# ---------------------------------------------------------------------------


def expected_allele_share(posterior: Sequence[float], ploidy: int) -> float:
    """The share of a sample's chromosome copies expected to carry the allele.

    Shares can be compared and averaged between samples of different ploidy.
    """
    return posterior_mean_dosage(posterior) / ploidy


def depth_ratio_frequencies(depths: Sequence[Sequence[int]]) -> Optional[List[float]]:
    """Starting allele frequencies: each allele's mean share of a sample's reads.

    These are polyRAD's starting frequencies (the depth ratios averaged by
    ``.alleleFreq``).
    ``depths`` is ``[sample][allele]``.  Samples with no reads at the marker
    are left out.  Returns None when no sample has a read.
    """
    allele_count = len(depths[0])
    share_totals = [0.0] * allele_count
    samples_with_reads = 0
    for sample_depths in depths:
        record_depth = sum(sample_depths)
        if record_depth == 0:
            continue
        samples_with_reads += 1
        for allele, allele_depth in enumerate(sample_depths):
            share_totals[allele] += allele_depth / record_depth
    if samples_with_reads == 0:
        return None
    return [total / samples_with_reads for total in share_totals]


# With the structure prior, the record's allele frequency over all samples counts as this
# many extra neighbours of every sample.  It keeps a sample's frequency off 0 and 1 when
# all its neighbours agree.
NEIGHBOUR_RECORD_WEIGHT = 1.0


def _frequency_priors(
    ploidy_groups: Sequence[int], allele_frequencies: Sequence[float], settings: ModelSettings
) -> Dict[int, List[List[float]]]:
    """Priors that follow from the allele frequencies alone, as ``[ploidy][allele][dosage]``.

    Uniform for the ``uniform`` prior and Hardy-Weinberg proportions for
    every other prior: ``hwe`` uses them throughout, the fitted priors start
    from them (``_fitted_priors``), and ``structure`` uses them until the
    neighbours take over.  Samples of the same ploidy share a prior.
    """
    if settings.prior == PRIOR_UNIFORM:
        return {ploidy: [uniform_prior(ploidy)] * len(allele_frequencies) for ploidy in ploidy_groups}
    return {
        ploidy: [hwe_prior(ploidy, frequency, settings.selfing_rate) for frequency in allele_frequencies]
        for ploidy in ploidy_groups
    }


def _fitted_priors(
    posteriors: Sequence[Sequence[Sequence[float]]],
    ploidies: Sequence[int],
    allele_frequencies: Sequence[float],
    previous_parameters: Dict[int, List[Tuple[float, float]]],
    free: bool,
) -> Tuple[Dict[int, List[List[float]]], Dict[int, List[Tuple[float, float]]]]:
    """Priors fitted to the current posteriors, one per ploidy group and allele.

    Returns the priors as ``[ploidy][allele][dosage]`` and, for the
    normal-shaped prior, their ``(mean, standard deviation)`` as
    ``[ploidy][allele]``.  With ``free`` the prior of a dosage is simply its
    share of the summed posteriors, with no shape imposed.  Each ploidy group
    is fitted on its own samples, because the dosages of one ploidy say
    nothing about another.  A group with fewer than
    ``FITTED_PRIOR_MIN_SAMPLES`` samples gets the Hardy-Weinberg prior for
    the marker's allele frequencies instead, and no parameters.
    """
    priors: Dict[int, List[List[float]]] = {}
    parameters: Dict[int, List[Tuple[float, float]]] = {}
    for ploidy in set(ploidies):
        group = [
            sample_posteriors
            for sample_posteriors, sample_ploidy in zip(posteriors, ploidies)
            if sample_ploidy == ploidy
        ]
        if len(group) < FITTED_PRIOR_MIN_SAMPLES:
            priors[ploidy] = [hwe_prior(ploidy, frequency) for frequency in allele_frequencies]
            continue
        priors[ploidy] = []
        if not free:
            parameters[ploidy] = []
        for allele, frequency in enumerate(allele_frequencies):
            dosage_weights = [
                sum(sample_posteriors[allele][dosage] for sample_posteriors in group) for dosage in range(ploidy + 1)
            ]
            if free:
                priors[ploidy].append([weight / len(group) for weight in dosage_weights])
                continue
            if ploidy in previous_parameters:
                start_mean, start_sd = previous_parameters[ploidy][allele]
            else:
                # The first fit starts from the binomial's mean and spread, as updog's does
                # (flexdog_full); the floor of 0.1 on the spread is added here.
                start_mean = ploidy * frequency
                start_sd = max(0.1, math.sqrt(ploidy * frequency * (1.0 - frequency)))
            mean, standard_deviation = fit_normal_prior(dosage_weights, ploidy, start_mean, start_sd)
            parameters[ploidy].append((mean, standard_deviation))
            priors[ploidy].append(normal_prior(ploidy, mean, standard_deviation))
    return priors, parameters


def _neighbourhood_frequencies(
    posteriors: Sequence[Sequence[Sequence[float]]],
    ploidies: Sequence[int],
    record_frequencies: Sequence[float],
    neighbours: Sequence[Sequence[int]],
) -> List[List[float]]:
    """Each sample's allele frequencies among its neighbours, as ``[allele][sample]``.

    A sample's frequency of an allele is the mean share of chromosome copies
    expected to carry it among the sample's neighbours, with the record's
    frequency over all samples counted as ``NEIGHBOUR_RECORD_WEIGHT`` more
    neighbours.  A sample is not among its own neighbours, so its own reads
    do not enter its own frequency.
    """
    frequencies = []
    for allele, record_frequency in enumerate(record_frequencies):
        shares = [
            expected_allele_share(sample_posteriors[allele], ploidy)
            for sample_posteriors, ploidy in zip(posteriors, ploidies)
        ]
        frequencies.append(
            [
                (sum(map(shares.__getitem__, sample_neighbours)) + NEIGHBOUR_RECORD_WEIGHT * record_frequency)
                / (len(sample_neighbours) + NEIGHBOUR_RECORD_WEIGHT)
                for sample_neighbours in neighbours
            ]
        )
    return frequencies


def _neighbour_priors(
    sample_frequencies: Sequence[Sequence[float]], ploidies: Sequence[int]
) -> Dict[int, List[List[float]]]:
    """Hardy-Weinberg priors at each sample's own allele frequencies, as ``[sample][allele][dosage]``.

    This is the form of polyRAD's population-structure prior
    (``AddGenotypePriorProb_ByTaxa``); where the frequencies come from is
    different (``find_neighbours``).
    """
    return {
        sample: [hwe_prior(ploidy, of_allele[sample]) for of_allele in sample_frequencies]
        for sample, ploidy in enumerate(ploidies)
    }


def _log_likelihoods_of_samples(
    depths: Sequence[Sequence[int]],
    record_depths: Sequence[int],
    count_constants: Sequence[Sequence[float]],
    ploidies: Sequence[int],
    allele_frequencies: Sequence[float],
    settings: ModelSettings,
    biases: Sequence[float],
) -> List[List[List[float]]]:
    """Log-likelihood of every dosage for every sample and allele, as ``[sample][allele][dosage]``.

    ``depths`` holds the alleles being computed, which need not be all the
    alleles of the record; ``record_depths`` is each sample's total over all
    of them.  ``count_constants`` is ``count_constant`` for every sample and
    allele, which does not change between passes.  The allele frequencies
    matter only through the contamination term.
    """
    sample_count = len(depths)
    allele_count = len(allele_frequencies)
    likelihood_frequencies = [frequency_for_likelihood(frequency, sample_count) for frequency in allele_frequencies]
    # Samples of the same ploidy share the beta-binomial of each dosage.
    shapes = {
        ploidy: [
            dosage_shapes(ploidy, likelihood_frequencies[allele], settings, biases[allele])
            for allele in range(allele_count)
        ]
        for ploidy in set(ploidies)
    }
    return [
        [
            log_likelihoods_for_counts(sample_depths[allele], record_depth, constants[allele], shapes[ploidy][allele])
            for allele in range(allele_count)
        ]
        for sample_depths, record_depth, constants, ploidy in zip(depths, record_depths, count_constants, ploidies)
    ]


def _posteriors_of_samples(
    log_likelihoods: Sequence[Sequence[Sequence[float]]],
    priors: Dict[int, List[List[float]]],
    prior_groups: Sequence[int],
) -> List[List[List[float]]]:
    """Posterior of every dosage for every sample and allele, as ``[sample][allele][dosage]``.

    ``priors`` is ``[group][allele][dosage]`` and ``prior_groups`` names each
    sample's group: its ploidy when samples of a ploidy share a prior, or the
    sample itself under the ``structure`` prior.
    """
    # The samples of a group share the prior's logs.
    log_priors = {group: [log_probabilities(prior) for prior in group_priors] for group, group_priors in priors.items()}
    return [
        [
            dosage_posterior(allele_log_likelihoods, priors[group][allele], log_priors[group][allele])
            for allele, allele_log_likelihoods in enumerate(sample_log_likelihoods)
        ]
        for sample_log_likelihoods, group in zip(log_likelihoods, prior_groups)
    ]


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
            totals[allele] += expected_allele_share(posterior, ploidy)
    return [total / len(posteriors) for total in totals]


def _iterate_record(
    depths: Sequence[Sequence[int]],
    record_depths: Sequence[int],
    count_constants: Sequence[Sequence[float]],
    ploidies: Sequence[int],
    settings: ModelSettings,
    frequencies: Sequence[float],
    biases: Sequence[float],
    warm_start: Optional[RecordCall] = None,
    neighbours: Optional[Sequence[Sequence[int]]] = None,
) -> RecordCall:
    """Run the passes of ``call_record`` for the alleles in ``depths`` with fixed biases.

    The passes are polyRAD's ``IterateHWE`` for one record, with a limit on
    their number.  For the ``norm`` prior they are the two steps of an EM
    iteration, as in updog, but stopped on the change in the frequencies and
    the prior where updog stops on the change in the log-likelihood.

    Each pass computes the posteriors from the current priors, re-estimates
    the allele frequencies from them, and renews the priors: from the new
    frequencies (``hwe``), by a fit to the posteriors (``norm`` and the
    estimate's free weights), or, when ``neighbours`` is given, from each
    sample's neighbours.  The passes stop when the quantity the priors follow
    from has settled: the record's allele frequencies, for the fitted priors
    also every prior probability, and with neighbours the samples' own
    frequencies.

    ``warm_start`` is an earlier result for the same alleles.  Without
    neighbours it is usually the result at a nearby bias, and its allele
    frequencies and priors replace ``frequencies`` and the priors that follow
    from them, so that fewer passes are needed.  With neighbours it is the
    record called with shared priors, which the samples' first frequencies
    are taken from.
    """
    sample_count = len(ploidies)
    normal_parameters: Dict[int, List[Tuple[float, float]]] = {}
    sample_frequencies: List[List[float]] = []
    if neighbours is not None:
        # One prior per sample.
        prior_groups = list(range(sample_count))
        frequencies = list(warm_start.allele_frequencies)
        sample_frequencies = _neighbourhood_frequencies(warm_start.posteriors, ploidies, frequencies, neighbours)
        priors = _neighbour_priors(sample_frequencies, ploidies)
    else:
        # One prior per ploidy.
        prior_groups = list(ploidies)
        if warm_start is None:
            frequencies = list(frequencies)
            priors = _frequency_priors(sorted(set(ploidies)), frequencies, settings)
        else:
            frequencies = list(warm_start.allele_frequencies)
            priors = warm_start.priors
            normal_parameters = warm_start.normal_prior_parameters

    # Only the contamination term makes the likelihood depend on the allele frequencies.
    # Without it the likelihoods are the same in every pass and are computed once.
    likelihoods_change = settings.contamination_rate > 0
    log_likelihoods = _log_likelihoods_of_samples(
        depths, record_depths, count_constants, ploidies, frequencies, settings, biases
    )
    if settings.prior == PRIOR_UNIFORM:
        posteriors = _posteriors_of_samples(log_likelihoods, priors, prior_groups)
        return RecordCall(
            True, frequencies, log_likelihoods, posteriors, iterations=1, priors=priors, prior_groups=prior_groups,
            allele_biases=list(biases),
        )

    for iterations in range(1, settings.max_iterations + 1):
        if likelihoods_change and iterations > 1:
            log_likelihoods = _log_likelihoods_of_samples(
                depths, record_depths, count_constants, ploidies, frequencies, settings, biases
            )
        posteriors = _posteriors_of_samples(log_likelihoods, priors, prior_groups)
        new_frequencies = frequencies_from_posteriors(posteriors, ploidies)
        used_priors, used_sample_frequencies = priors, sample_frequencies
        if neighbours is not None:
            sample_frequencies = _neighbourhood_frequencies(posteriors, ploidies, new_frequencies, neighbours)
            priors = _neighbour_priors(sample_frequencies, ploidies)
            change = sum(
                abs(new - old)
                for new_allele, old_allele in zip(sample_frequencies, used_sample_frequencies)
                for new, old in zip(new_allele, old_allele)
            ) / (len(new_frequencies) * sample_count)
        else:
            change = sum(abs(new - old) for new, old in zip(new_frequencies, frequencies)) / len(frequencies)
            if settings.prior in (PRIOR_NORM, PRIOR_FREE):
                priors, normal_parameters = _fitted_priors(
                    posteriors, ploidies, new_frequencies, normal_parameters, free=settings.prior == PRIOR_FREE
                )
                largest_prior_change = max(
                    abs(new - old)
                    for ploidy in priors
                    for new_prior, old_prior in zip(priors[ploidy], used_priors[ploidy])
                    for new, old in zip(new_prior, old_prior)
                )
                change = max(change, largest_prior_change)
            else:
                priors = _frequency_priors(sorted(priors), new_frequencies, settings)
        frequencies = new_frequencies
        if change <= settings.tolerance:
            break
    return RecordCall(
        True,
        frequencies,
        log_likelihoods,
        posteriors,
        iterations=iterations,
        final_change=change,
        converged=change <= settings.tolerance,
        priors=used_priors,
        prior_groups=prior_groups,
        normal_prior_parameters=normal_parameters,
        allele_biases=list(biases),
        # As [sample][allele], the frequencies the posteriors were computed at.
        sample_frequencies=[list(of_sample) for of_sample in zip(*used_sample_frequencies)],
    )


# updog's settings for the allelic bias: the five values each fit is started from, and
# the spread of the normal prior that pulls the log of the bias towards 0 (no bias).
BIAS_LOG_STARTS = (-1.0, -0.5, 0.0, 0.5, 1.0)
BIAS_LOG_PRIOR_SD = 0.7
# The search goes no further than a bias of exp(-2) to exp(2), about 0.14 to 7.4, and
# stops when the log of the bias is known to this precision (about 1% of the bias).
BIAS_LOG_LIMIT = 2.0
BIAS_LOG_TOLERANCE = 0.01
# With --allele-bias auto a record is called with its fitted bias only when the fit raises
# the log-likelihood of its read counts, after the prior on the bias, by this much over no
# bias (``ModelSettings.min_bias_gain``).  The value is MADCATS's own, chosen on potato,
# delta smelt and polyploids pooled from smelt reads: at 10 the records of biased data are
# called with their bias and nine in ten records of unbiased data without (docs/DEVELOPMENT.md).
BIAS_AUTO_MIN_GAIN = 10.0
# Why a record was called without the bias fitted to it (``BiasFit.rejection``).  A bias at
# a limit of the search is not used with --allele-bias auto
# (``ModelSettings.reject_bias_at_limit``): where samples have few heterozygotes, as inbred
# lines do, and the prior expects them, a bias at the limit makes homozygous samples pass
# for heterozygous ones and gains far more than the margin (docs/DEVELOPMENT.md).
BIAS_REJECTED_FOR_GAIN = "gain_below_margin"
BIAS_REJECTED_AT_LIMIT = "at_search_limit"
# A warning is given when this share of a run's fitted biases end at a limit of the search.
# On five real data sets the share was 0.3 to 2.5%; on simulated inbred lines, all of them.
BIAS_AT_LIMIT_WARNING_SHARE = 0.2


class BiasFit(NamedTuple):
    """How the allelic bias of a two-allele record was fitted.

    ``log_bias`` is the natural log of the alternate allele's fitted bias, and
    ``objective`` the value that was maximised there: the log-likelihood of
    the read counts plus the log of the prior on the bias.
    ``nearest_log_start`` is the log bias of the starting value nearest the
    result, ``fits`` the number of bias values tried, and ``at_limit``
    whether the result is the largest or smallest bias the search allows.
    ``gain`` is how much the fitted bias raised the objective over no bias.
    ``rejection`` is None when the record's call was made with the bias.
    Otherwise it says why not, and the record's ``allele_biases`` are 1
    although a bias was fitted: ``BIAS_REJECTED_FOR_GAIN`` where the gain is
    less than ``ModelSettings.min_bias_gain``, ``BIAS_REJECTED_AT_LIMIT``
    where the gain is enough but the bias is at a limit of the search and
    ``ModelSettings.reject_bias_at_limit`` is set.
    """

    log_bias: float
    objective: float
    nearest_log_start: float
    fits: int
    at_limit: bool
    gain: float
    rejection: Optional[str] = None

    @property
    def used(self) -> bool:
        return self.rejection is None


def log_probability_of_counts(log_likelihoods: Sequence[float], log_prior: Sequence[float]) -> float:
    """Log-probability of one sample's read count of an allele, with the sample's genotype summed out.

    The genotype is unknown, so the probability of the count is the average
    over the dosages, weighted by their prior probabilities.  Minus infinity
    when the prior rules out every dosage the count allows.
    """
    return log_sum_of_exponentials(
        [log_likelihood + log_dosage_prior for log_likelihood, log_dosage_prior in zip(log_likelihoods, log_prior)]
    )


def marginal_log_likelihood(
    log_likelihoods: Sequence[Sequence[Sequence[float]]],
    priors: Dict[int, List[List[float]]],
    prior_groups: Sequence[int],
    allele: int = 0,
) -> float:
    """Log-likelihood of the read counts of one allele, with the genotypes summed out.

    The sum over samples of ``log_probability_of_counts``.
    ``log_likelihoods`` is ``[sample][allele][dosage]``,
    ``priors`` is ``[group][allele][dosage]`` and ``prior_groups`` gives each
    sample's group, as in ``RecordCall``.
    """
    log_priors = {group: log_probabilities(prior[allele]) for group, prior in priors.items()}
    total = 0.0
    for sample_log_likelihoods, group in zip(log_likelihoods, prior_groups):
        sample_total = log_probability_of_counts(sample_log_likelihoods[allele], log_priors[group])
        # A sample whose counts the prior rules out entirely is left out, as its posterior
        # falls back to the prior (see dosage_posterior).
        if sample_total > NEGATIVE_INFINITY:
            total += sample_total
    return total


def fit_allele_bias(
    depths: Sequence[Sequence[int]],
    record_depths: Sequence[int],
    count_constants: Sequence[Sequence[float]],
    ploidies: Sequence[int],
    settings: ModelSettings,
    frequencies: Sequence[float],
) -> RecordCall:
    """Call one allele of a two-allele record while estimating its allelic bias.

    ``depths`` holds the one allele (``[sample][0]``).  The bias (see
    ``expected_read_fraction``) is the value that makes the read counts most
    probable, after a normal prior on its log that pulls it towards 1, as in
    updog.  For any given bias the record is called as usual
    (``_iterate_record``), which also fits the allele frequency or the normal
    prior, and the log-likelihood of the counts is taken from that result.

    The bias, its prior and the starting values are updog's (``flexdog``).
    The search is MADCATS's own: updog runs its EM, which updates the bias,
    the sequencing error and the overdispersion together, from each starting
    value.

    The search has two stages.  The record is first called at each of
    updog's five starting values, each from scratch, so that a wrong first
    guess of the genotypes at one bias cannot carry over to another.  Every
    starting value that scores at least as well as its neighbours is then
    refined by a golden-section search between them, each call starting from
    the nearest one already made; an outermost starting value is first moved
    outwards in the same steps for as long as that improves the fit.  The
    result is the call with the highest value among all that were made.

    Several starting values are refined, not just the best, because their
    scores before refinement can mislead.  At a marker where all samples
    share one or two dosages, the reads fit equally well with every dosage
    shifted by one and a strong bias; only the prior on the bias tells the
    two apart, and only once each has been refined to its best bias.

    With ``settings.min_bias_gain`` the best bias is used only when it raises
    that value by as much over no bias, which is always among the calls made;
    otherwise the record's result is its call without a bias, and the fitted
    bias is only reported (``BiasFit``).  A bias fitted to a record that has
    none follows the noise in its counts, and on reads without bias this cost
    accuracy, most of all in hexaploids.  With
    ``settings.reject_bias_at_limit`` a bias at a limit of the search is not
    used either (see ``BIAS_REJECTED_AT_LIMIT``).  Both rules are MADCATS's
    own; updog always uses the bias it fits.
    """
    tried: Dict[float, Tuple[float, RecordCall]] = {}

    def objective(log_bias: float, from_scratch: bool = False) -> float:
        if log_bias not in tried:
            warm_start = None
            if tried and not from_scratch:
                warm_start = tried[min(tried, key=lambda other: abs(other - log_bias))][1]
            result = _iterate_record(
                depths, record_depths, count_constants, ploidies, settings, frequencies, [math.exp(log_bias)],
                warm_start,
            )
            log_prior = -(log_bias**2) / (2.0 * BIAS_LOG_PRIOR_SD**2)
            log_likelihood = marginal_log_likelihood(result.log_likelihoods, result.priors, result.prior_groups)
            tried[log_bias] = (log_likelihood + log_prior, result)
        return tried[log_bias][0]

    step = BIAS_LOG_STARTS[1] - BIAS_LOG_STARTS[0]
    start_scores = [objective(log_bias, from_scratch=True) for log_bias in BIAS_LOG_STARTS]
    last = len(BIAS_LOG_STARTS) - 1
    for index, centre in enumerate(BIAS_LOG_STARTS):
        neighbour_scores = [start_scores[other] for other in (index - 1, index + 1) if 0 <= other <= last]
        if start_scores[index] < max(neighbour_scores):
            continue
        if index in (0, last):
            while abs(centre) + step <= BIAS_LOG_LIMIT + 1e-9:
                further = centre + math.copysign(step, centre)
                if objective(further) <= objective(centre):
                    break
                centre = further
        _golden_section_maximum(
            objective, max(-BIAS_LOG_LIMIT, centre - step), min(BIAS_LOG_LIMIT, centre + step), BIAS_LOG_TOLERANCE
        )
    log_bias = max(tried, key=lambda value: tried[value][0])
    best_objective, result = tried[log_bias]
    gain = best_objective - tried[0.0][0]
    at_limit = abs(log_bias) >= BIAS_LOG_LIMIT
    rejection = None
    if gain < settings.min_bias_gain:
        rejection = BIAS_REJECTED_FOR_GAIN
    elif at_limit and settings.reject_bias_at_limit:
        rejection = BIAS_REJECTED_AT_LIMIT
    if rejection is not None:
        result = tried[0.0][1]
    # The starting value whose refinement gave the fitted bias: the nearest one.
    nearest_log_start = min(BIAS_LOG_STARTS, key=lambda start: abs(start - log_bias))
    result.bias_fit = BiasFit(log_bias, best_objective, nearest_log_start, len(tried), at_limit, gain, rejection)
    return result


def _with_mirrored_reference_allele(alternate: RecordCall) -> RecordCall:
    """Complete a two-allele record from the result for its alternate allele.

    With two alleles the reference allele's dosage is the ploidy minus the
    alternate allele's, so everything about it is the alternate allele's
    result read backwards: its likelihoods, prior and posterior in reverse
    dosage order, its frequency one minus the other, and its bias the
    reciprocal.
    """
    def both(values: Sequence[float]) -> List[List[float]]:
        return [list(reversed(values)), list(values)]

    return dataclasses.replace(
        alternate,
        allele_frequencies=[1.0 - alternate.allele_frequencies[0], alternate.allele_frequencies[0]],
        log_likelihoods=[both(sample[0]) for sample in alternate.log_likelihoods],
        posteriors=[both(sample[0]) for sample in alternate.posteriors],
        priors={ploidy: both(prior[0]) for ploidy, prior in alternate.priors.items()},
        normal_prior_parameters={
            ploidy: [(ploidy - parameters[0][0], parameters[0][1]), parameters[0]]
            for ploidy, parameters in alternate.normal_prior_parameters.items()
        },
        allele_biases=[1.0 / alternate.allele_biases[0], alternate.allele_biases[0]],
        sample_frequencies=[[1.0 - of_sample[0], of_sample[0]] for of_sample in alternate.sample_frequencies],
    )


def call_record(
    depths: Sequence[Sequence[int]],
    ploidies: Sequence[int],
    settings: ModelSettings = ModelSettings(),
) -> RecordCall:
    """Estimate dosage posteriors for every sample and allele of one marker.

    depths: read counts, ``[sample][allele]``.
    ploidies: one ploidy per sample.

    With the ``uniform`` prior there is one pass, and the allele frequencies
    are used only by the contamination term.  With the ``hwe`` prior the
    passes repeat until the frequencies settle (see ``ModelSettings``).  With
    the ``norm`` prior each pass also re-fits the prior to the posteriors
    (the two steps of an EM iteration), and the passes repeat until neither
    the frequencies nor any prior probability changes by more than the
    tolerance.  With the ``structure`` prior the record is first called with
    the ``hwe`` prior, and the passes are then run again with each sample's
    prior taken from its neighbours (``_neighbourhood_frequencies``).

    A record with two alleles is computed for its second (alternate) allele
    only; the first is its mirror image (``_with_mirrored_reference_allele``).
    With ``settings.fit_allele_bias`` such a record also gets an allelic bias
    (``fit_allele_bias``); records with more alleles never do.  A record that
    is called without a bias, because it has more than two alleles or because
    its fitted bias is not used, is called with ``settings.without_bias()``:
    exactly as a run that fits no bias calls it.
    """
    if len(depths) != len(ploidies):
        raise ValueError("depths and ploidies must have one entry per sample")
    # Without any neighbours the structure prior is the Hardy-Weinberg prior.
    use_neighbours = settings.prior == PRIOR_STRUCTURE and settings.neighbours is not None and any(settings.neighbours)
    if use_neighbours and len(settings.neighbours) != len(ploidies):
        raise ValueError("settings.neighbours must have one entry per sample")
    frequencies = depth_ratio_frequencies(depths)
    if frequencies is None:
        return RecordCall(has_reads=False, allele_frequencies=[], log_likelihoods=[], posteriors=[])

    record_depths = [sum(sample_depths) for sample_depths in depths]
    two_alleles = len(frequencies) == 2
    if two_alleles:
        depths = [[sample_depths[1]] for sample_depths in depths]
        frequencies = [frequencies[1]]
    elif settings.fit_allele_bias:
        settings = settings.without_bias()

    def count_constants_at(overdispersion: float) -> List[List[float]]:
        return [
            [count_constant(allele_depth, record_depth, overdispersion) if record_depth else 0.0
             for allele_depth in sample_depths]
            for sample_depths, record_depth in zip(depths, record_depths)
        ]

    count_constants = count_constants_at(settings.overdispersion)
    if two_alleles and settings.fit_allele_bias:
        result = fit_allele_bias(depths, record_depths, count_constants, ploidies, settings, frequencies)
        if not result.bias_fit.used and settings.read_model_without_bias is not None:
            fit = result.bias_fit
            settings = settings.without_bias()
            count_constants = count_constants_at(settings.overdispersion)
            result = _iterate_record(depths, record_depths, count_constants, ploidies, settings, frequencies, [1.0])
            result.bias_fit = fit
    else:
        # Without neighbours, _iterate_record gives the structure prior Hardy-Weinberg
        # proportions, which is where its second run starts from.
        result = _iterate_record(
            depths, record_depths, count_constants, ploidies, settings, frequencies, [1.0] * len(frequencies)
        )
    if use_neighbours:
        shared = result
        result = _iterate_record(
            depths, record_depths, count_constants, ploidies, settings, shared.allele_frequencies,
            shared.allele_biases, warm_start=shared, neighbours=settings.neighbours,
        )
        result.iterations += shared.iterations
        result.converged = result.converged and shared.converged
    return _with_mirrored_reference_allele(result) if two_alleles else result


# ---------------------------------------------------------------------------
# From marginal posteriors to one genotype per sample
# ---------------------------------------------------------------------------


def most_probable_dosage(posterior: Sequence[float]) -> Optional[int]:
    """The dosage with the highest posterior probability, or None when two tie for it.

    polyRAD's ``BestGenos`` takes the first of tied dosages, and compares in
    single precision, so that it also ties posteriors such as 0.4, 0.4, 0.2.
    """
    highest = max(posterior)
    if sum(1 for probability in posterior if probability == highest) > 1:
        return None
    return list(posterior).index(highest)


def correct_dosages_to_ploidy(
    dosages: Sequence[Optional[int]],
    posteriors: Sequence[Sequence[float]],
    ploidy: int,
) -> Optional[List[int]]:
    """Adjust the dosages of one sample's alleles until they add up to ``ploidy``.

    ``dosages`` holds each allele's most probable dosage and ``posteriors``
    its dosage posterior.  This is ``CorrectGenos`` of polyRAD 2.0 and later
    (``GetProbableGenotypes(multiallelic = "correct")``): while the sum
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
            # The log-probability of the whole genotype if this allele takes the step.
            genotype_log_probability = sum(
                log_probability(posteriors[other][corrected[other] + (step if other == allele else 0)])
                for other in range(len(corrected))
            )
            if genotype_log_probability > best_log_probability:
                best_allele, best_log_probability = allele, genotype_log_probability
        if best_allele is None:
            return None
        corrected[best_allele] += step
    return corrected


def best_dosages_summing_to_ploidy(posteriors: Sequence[Sequence[float]], ploidy: int) -> Optional[List[int]]:
    """The dosages that add up to ``ploidy`` with the highest product of marginal posteriors.

    It maximises the same quantity as ``correct_dosages_to_ploidy`` but over
    every combination, so its result is never less probable.  polyRAD had an
    exhaustive search for the same quantity (``BestMultiGeno``, deprecated in
    2.0); this one is a dynamic programme and is used only to count how often
    polyRAD's rule falls short of it.  The search goes
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
# Estimating the read model from the data: overdispersion and background-read rate
# ---------------------------------------------------------------------------

# Values tried by ``estimate_overdispersion``.  Steps widen as the values grow, because
# the beta-binomial changes less and less as it approaches the binomial.
OVERDISPERSION_GRID = (2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500)
# The estimate uses at most this many sample-by-marker cells, from markers taken at even
# spacing through the table, so that it stays quick and the same table always gives the
# same value.
OVERDISPERSION_MAX_CELLS = 20000
# However many samples there are, at least this many markers are used when the table has them.
OVERDISPERSION_MIN_MARKERS = 20
# With fewer sample-by-marker cells that have reads, the estimate would be mostly noise.
OVERDISPERSION_MIN_CELLS = 100
OVERDISPERSION_MAX_ROUNDS = 4
# Background-read rates tried by ``estimate_read_model``.
BACKGROUND_RATE_GRID = (0.0002, 0.0005, 0.001, 0.002, 0.003, 0.005, 0.01, 0.02)


@dataclasses.dataclass
class OverdispersionEstimate:
    """Result of ``estimate_overdispersion``.

    ``value`` is the chosen overdispersion, or None when fewer than
    ``OVERDISPERSION_MIN_CELLS`` sample-by-marker cells had reads
    (``cells_with_reads``), which is too little to estimate from.
    ``log_likelihood_by_value`` maps every value tried to the log-likelihood
    of the read counts under it, in the last round.  ``at_limit`` is ``"lowest"`` or
    ``"highest"`` when the chosen value is an end of the grid, so the best
    value may lie beyond it.  ``values_by_round`` is the value chosen in each
    round, and ``markers_used`` the number of markers the estimate is based on.
    """

    value: Optional[float]
    log_likelihood_by_value: Dict[float, float]
    markers_used: int
    values_by_round: List[float]
    at_limit: Optional[str] = None
    cells_with_reads: int = 0


def marker_overdispersion_log_likelihoods(
    depths: Sequence[Sequence[int]], ploidies: Sequence[int], settings: ModelSettings, grid: Sequence[float]
) -> Optional[List[float]]:
    """Log-likelihood of one two-allele marker's read counts at each overdispersion in ``grid``.

    depths: ``[sample][allele]`` read counts of the reference and alternate allele.

    The marker is first called with ``settings``.  The allele frequency, the
    genotype prior and the allelic bias of that call are then held fixed
    while the overdispersion takes each value of ``grid``, and the
    log-likelihood of the counts is computed with the genotypes summed out
    (as in ``marginal_log_likelihood``).  Returns None for a marker without
    reads.
    """
    result = call_record(depths, ploidies, settings)
    if not result.has_reads:
        return None
    frequency = frequency_for_likelihood(result.allele_frequencies[1], len(depths))
    bias = result.allele_biases[1]
    log_priors = {group: log_probabilities(prior[1]) for group, prior in result.priors.items()}
    # The samples with reads, each with the term of count_constant that no overdispersion changes.
    samples_with_reads = [
        (alternate_depth, reference_depth + alternate_depth,
         log_binomial_coefficient(reference_depth + alternate_depth, alternate_depth), ploidy, group)
        for (reference_depth, alternate_depth), ploidy, group in zip(depths, ploidies, result.prior_groups)
        if reference_depth + alternate_depth > 0
    ]
    totals = []
    for value in grid:
        at_value = dataclasses.replace(settings, overdispersion=value)
        shapes = {ploidy: dosage_shapes(ploidy, frequency, at_value, bias) for ploidy in set(ploidies)}
        log_gamma_of_value = math.lgamma(value)
        total = 0.0
        for alternate_depth, record_depth, log_coefficient, ploidy, group in samples_with_reads:
            # count_constant, with its terms in the same order.
            constant = log_coefficient + log_gamma_of_value - math.lgamma(record_depth + value)
            log_likelihoods = log_likelihoods_for_counts(alternate_depth, record_depth, constant, shapes[ploidy])
            sample_total = log_probability_of_counts(log_likelihoods, log_priors[group])
            if sample_total > NEGATIVE_INFINITY:
                total += sample_total
        totals.append(total)
    return totals


def _evenly_spaced(items: Sequence, limit: int) -> List:
    """At most ``limit`` of the items, taken at even spacing."""
    if len(items) <= limit:
        return list(items)
    return [items[index * len(items) // limit] for index in range(limit)]


def _markers_for_estimate(
    depths_by_marker: Sequence[Sequence[Sequence[int]]], sample_count: int
) -> Tuple[List[Sequence[Sequence[int]]], int]:
    """The markers the read-model estimate works on, and how many of their sample cells have reads.

    Markers without reads are left out, and at most enough markers to make up
    ``OVERDISPERSION_MAX_CELLS`` cells are kept, at even spacing.  Selecting
    from an earlier selection returns it unchanged.
    """
    samples_with_reads = [
        sum(1 for sample_depths in depths if sum(sample_depths) > 0) for depths in depths_by_marker
    ]
    usable = [(depths, count) for depths, count in zip(depths_by_marker, samples_with_reads) if count > 0]
    marker_limit = max(OVERDISPERSION_MIN_MARKERS, OVERDISPERSION_MAX_CELLS // max(1, sample_count))
    chosen_markers = _evenly_spaced(usable, marker_limit)
    return [depths for depths, _ in chosen_markers], sum(count for _, count in chosen_markers)


def estimate_overdispersion(
    depths_by_marker: Sequence[Sequence[Sequence[int]]],
    ploidies: Sequence[int],
    settings: ModelSettings = ModelSettings(),
    grid: Sequence[float] = OVERDISPERSION_GRID,
    start: float = POLYRAD_OVERDISPERSION,
    parallel_map=map,
) -> OverdispersionEstimate:
    """Choose the overdispersion that makes the observed read counts most probable.

    depths_by_marker: ``[marker][sample][allele]`` read counts of two alleles
        per marker (the reference-type and alternate-type rows).
    settings: the model the markers are called with while estimating; its
        overdispersion is ignored.  The command passes the free-weights prior
        whatever prior the run uses (``choose_read_model``).
    parallel_map: a function like ``map``, used to work on the markers; a
        process pool's ``map`` spreads them over its workers.

    For every value in ``grid`` the log-likelihood of the counts is added up
    over the markers (``marker_overdispersion_log_likelihoods``), and the
    value with the highest total is chosen.  Each marker's allele frequency,
    genotype prior and allelic bias come from calling it under the current
    value, so the search is repeated with the chosen value until the choice
    no longer changes, which usually takes two rounds.

    This is a maximum-likelihood estimate and MADCATS's own; updog instead
    fits an overdispersion for every marker.  polyRAD's ``TestOverdispersion``
    instead tests the read counts of samples that a rough genotype estimate
    assigns one allele copy.  That works in diploids, but in polyploids the
    rough estimate assigns the wrong dosage often enough that the counts look
    more scattered than they are (see docs/DEVELOPMENT.md).
    """
    markers, cells_with_reads = _markers_for_estimate(depths_by_marker, len(ploidies))
    if cells_with_reads < OVERDISPERSION_MIN_CELLS:
        return OverdispersionEstimate(None, {}, len(markers), [], cells_with_reads=cells_with_reads)

    current = start
    values_by_round: List[float] = []
    log_likelihood: Dict[float, float] = {}
    for _ in range(OVERDISPERSION_MAX_ROUNDS):
        per_marker = parallel_map(
            functools.partial(
                marker_overdispersion_log_likelihoods,
                ploidies=ploidies,
                settings=dataclasses.replace(settings, overdispersion=current),
                grid=grid,
            ),
            markers,
        )
        totals = [0.0] * len(grid)
        for marker_totals in per_marker:
            if marker_totals is not None:
                totals = [total + marker_total for total, marker_total in zip(totals, marker_totals)]
        log_likelihood = dict(zip(grid, totals))
        best = max(grid, key=lambda value: log_likelihood[value])
        values_by_round.append(best)
        if best == current:
            break
        current = best
    chosen = values_by_round[-1]
    at_limit = "lowest" if chosen == min(grid) else "highest" if chosen == max(grid) else None
    return OverdispersionEstimate(chosen, log_likelihood, len(markers), values_by_round, at_limit, cells_with_reads)


@dataclasses.dataclass
class ReadModelEstimate:
    """Result of ``estimate_read_model``.

    ``background_rate`` is the chosen rate and ``overdispersion_estimate``
    the overdispersion estimate made at that rate.  ``by_background_rate`` maps
    every rate tried to the overdispersion chosen at it and the
    log-likelihood of the read counts there.
    """

    background_rate: float
    overdispersion_estimate: OverdispersionEstimate
    by_background_rate: Dict[float, Tuple[Optional[float], float]]


def estimate_read_model(
    depths_by_marker: Sequence[Sequence[Sequence[int]]],
    ploidies: Sequence[int],
    settings: ModelSettings,
    background_rates: Sequence[float] = BACKGROUND_RATE_GRID,
    overdispersions: Sequence[float] = OVERDISPERSION_GRID,
    parallel_map=map,
    start_background_rate: float = FALLBACK_BACKGROUND_RATE,
    start_overdispersion: float = POLYRAD_OVERDISPERSION,
) -> ReadModelEstimate:
    """Choose the background-read rate and the overdispersion that make the read counts most probable.

    For a given background-read rate the overdispersion is chosen by
    ``estimate_overdispersion``, which also gives the log-likelihood of the
    counts at its choice.  The rate with the highest such log-likelihood is
    found by starting at the rate nearest ``start_background_rate`` and
    moving to a neighbouring rate for as long as that raises it; the search
    for the overdispersion at the first rate starts from
    ``start_overdispersion``.  Give a single value in ``background_rates`` or
    ``overdispersions`` to hold that parameter fixed.

    The two parameters can be told apart because they act on different
    samples: reads of an allele in samples that do not carry it show the
    background rate, and the spread of read fractions among samples that
    carry both alleles shows the overdispersion.  In simulated diploids and
    tetraploids with 20 or more reads per sample both were recovered; with 6
    to 25 reads the rate was still recovered, and the overdispersion was
    within a factor of two of a simulated 40 and unmeasurable at 150 (see
    docs/DEVELOPMENT.md).
    """
    rates = sorted(background_rates)
    tried: Dict[float, OverdispersionEstimate] = {}
    # Every rate is judged on the same markers, selected once.
    depths_by_marker, _ = _markers_for_estimate(depths_by_marker, len(ploidies))

    def log_likelihood_at(index: int) -> float:
        rate = rates[index]
        if rate not in tried:
            # Start each search from the overdispersion chosen at the rate tried before.
            start = start_overdispersion
            if tried:
                previous = tried[min(tried, key=lambda other: abs(other - rate))]
                start = previous.value if previous.value is not None else start
            if len(overdispersions) == 1:
                start = overdispersions[0]
            tried[rate] = estimate_overdispersion(
                depths_by_marker, ploidies, dataclasses.replace(settings, background_rate=rate),
                overdispersions, start, parallel_map,
            )
        estimate = tried[rate]
        if not estimate.log_likelihood_by_value:
            return NEGATIVE_INFINITY
        return max(estimate.log_likelihood_by_value.values())

    best = min(range(len(rates)), key=lambda index: abs(rates[index] - start_background_rate))
    log_likelihood_at(best)
    while True:
        neighbours = [index for index in (best - 1, best + 1) if 0 <= index < len(rates)]
        better = [index for index in neighbours if log_likelihood_at(index) > log_likelihood_at(best)]
        if not better:
            break
        best = max(better, key=log_likelihood_at)
    return ReadModelEstimate(
        rates[best],
        tried[rates[best]],
        {rate: (estimate.value, log_likelihood_at(rates.index(rate))) for rate, estimate in sorted(tried.items())},
    )


# ---------------------------------------------------------------------------
# Population structure: the samples most like each sample
# ---------------------------------------------------------------------------

# The numbers of neighbours tried; 0 is the plain Hardy-Weinberg prior.
NEIGHBOUR_COUNTS = (0, 5, 10, 20, 40, 80)
# Samples are compared on at most this many markers, taken at even spacing.
NEIGHBOUR_MAX_MARKERS = 1000
# Below these sizes the table is too small to tell which samples are alike.
NEIGHBOUR_MIN_MARKERS = 100
NEIGHBOUR_MIN_SAMPLES = 20
# The number of neighbours is judged on at most this many markers of each half.
NEIGHBOUR_JUDGED_MARKERS = 250
# The samples are handed to the worker processes in about this many blocks.
NEIGHBOUR_BLOCKS = 16


@dataclasses.dataclass
class NeighbourChoice:
    """Result of ``find_neighbours``.

    ``neighbours`` holds, for each sample, the indices of the samples most
    like it, nearest first.  Every sample has ``neighbours_per_sample`` of
    them, which is 0 when neighbours predict a sample no better than the
    whole table does, or when the table is too small to tell; ``reason``
    then says which.
    ``markers_used`` is the number of markers the samples were compared on.
    ``error_by_count`` maps each number of neighbours tried to the mean
    squared error with which that many neighbours predicted the samples'
    allele shares at markers they had not been chosen on.
    """

    neighbours: Tuple[Tuple[int, ...], ...]
    neighbours_per_sample: int
    markers_used: int
    error_by_count: Dict[int, float]
    reason: Optional[str] = None


def expected_alternate_shares(
    depths: Sequence[Sequence[int]], ploidies: Sequence[int], settings: ModelSettings
) -> Optional[List[float]]:
    """Each sample's expected share of chromosome copies that carry a marker's second allele.

    ``depths`` is ``[sample][allele]`` with two alleles.  The share is the
    posterior mean dosage divided by the ploidy, so that samples of different
    ploidy can be compared.  A sample without reads gets the marker's allele
    frequency.  Returns None when no sample has a read.
    """
    result = call_record(depths, ploidies, settings)
    if not result.has_reads:
        return None
    return [
        expected_allele_share(sample_posteriors[1], ploidy)
        for sample_posteriors, ploidy in zip(result.posteriors, ploidies)
    ]


def _sum_of_products(first: Sequence[float], second: Sequence[float]) -> float:
    return sum(map(operator.mul, first, second))


if hasattr(math, "sumprod"):          # Python 3.12 and later: several times faster
    _sum_of_products = math.sumprod   # noqa: F811


def _neighbour_candidates(
    samples: Sequence[int],
    shares_by_half: Sequence[Sequence[Sequence[float]]],
    judged_by_half: Sequence[Sequence[Sequence[float]]],
    judged_means: Sequence[Sequence[float]],
    neighbour_counts: Sequence[int],
) -> List[Tuple[List[int], List[float]]]:
    """For each of ``samples``: the samples most like it, and how well its nearest ones predict it.

    ``shares_by_half[half][sample]`` holds a sample's allele shares at the
    markers of one half of the table, and ``judged_by_half`` the same for a
    thinned set of each half's markers, whose means over all samples are
    ``judged_means``.  Two samples are alike when their shares differ little:
    the measure is the sum of squared differences over markers.

    Returns, per sample, its ``max(neighbour_counts)`` most similar samples over both
    halves together, nearest first, and for each entry of ``neighbour_counts``
    (ascending) the summed squared error of predicting the sample's shares at
    the judged markers of one half from that many of its nearest samples as
    found on the other half, both ways round.  The prediction is the one the
    prior uses: the neighbours' mean with the marker's mean counted as
    ``NEIGHBOUR_RECORD_WEIGHT`` more of them.
    """
    # The rows arrive as compact arrays; products of plain lists are about twice as fast.
    shares_by_half = [[list(row) for row in half] for half in shares_by_half]
    sample_count = len(shares_by_half[0])
    limit = max(neighbour_counts)
    squared_lengths = [[_sum_of_products(row, row) for row in half] for half in shares_by_half]
    results = []
    for sample in samples:
        others = [other for other in range(sample_count) if other != sample]
        differences_by_half = []
        for half, lengths in zip(shares_by_half, squared_lengths):
            own = half[sample]
            # |a - b|^2 = |a|^2 + |b|^2 - 2 a.b
            differences_by_half.append(
                [lengths[sample] + lengths[other] - 2.0 * _sum_of_products(own, half[other])
                 for other in range(sample_count)]
            )
        errors = [0.0] * len(neighbour_counts)
        for measured, judged in ((0, 1), (1, 0)):
            nearest = heapq.nsmallest(limit, others, key=differences_by_half[measured].__getitem__)
            own = judged_by_half[judged][sample]
            means = judged_means[judged]
            neighbour_totals = [0.0] * len(own)
            included = 0
            for index, count in enumerate(neighbour_counts):
                for other in nearest[included:count]:
                    neighbour_totals = list(map(operator.add, neighbour_totals, judged_by_half[judged][other]))
                included = count
                weight = count + NEIGHBOUR_RECORD_WEIGHT
                errors[index] += sum(
                    (value - (total + NEIGHBOUR_RECORD_WEIGHT * mean) / weight) ** 2
                    for value, total, mean in zip(own, neighbour_totals, means)
                )
        overall = [first + second for first, second in zip(*differences_by_half)]
        results.append((heapq.nsmallest(limit, others, key=overall.__getitem__), errors))
    return results


def find_neighbours(
    depths_by_marker: Sequence[Sequence[Sequence[int]]],
    ploidies: Sequence[int],
    settings: ModelSettings,
    parallel_map=map,
) -> NeighbourChoice:
    """Find the samples most like each sample across markers, and how many of them to use.

    This is the one step of the ``structure`` prior that looks at the whole
    table, and it is MADCATS's own.  polyRAD's population-structure prior
    (``IteratePopStruct``) gets each sample's allele frequencies from a
    principal component analysis of all markers and a regression on its
    axes; here they come from the sample's neighbours.  The prior built from
    them has polyRAD's form (``_neighbour_priors``).  ``depths_by_marker`` is ``[marker][sample][allele]`` with two
    alleles per marker (reference and alternate).

    Every marker is called with the ``hwe`` prior, which gives each sample an
    expected share of alternate copies per marker.  Samples whose shares
    differ little over the markers are alike: members of one subpopulation,
    of one family, or close on a cline.  The sum of squared differences used
    here ranks pairs of samples as a genomic relationship matrix would.

    How many neighbours to use is chosen from ``NEIGHBOUR_COUNTS`` by
    prediction.  The markers are split in two halves; the neighbours found on
    one half predict each sample's shares at the other half, and the number
    with the smallest error wins.  With too many neighbours a small group or
    family is diluted, and with too few the frequency is noisy.  When no
    number predicts better than the table's own frequencies (0 neighbours)
    there is no structure worth modelling and the prior stays Hardy-Weinberg.
    The halves matter: judged on the markers they were found on, neighbours
    always look better than they are, because each of those markers helped
    pick them.

    On simulated groups, clines and families and on thinned potato DArTag
    counts, the prior built from these lists called as many genotypes
    correctly as polyRAD's population-structure prior or more (see
    docs/DEVELOPMENT.md).

    ``parallel_map`` is a replacement for ``map``, for calling the markers and
    comparing the samples in worker processes.
    """
    sample_count = len(ploidies)
    no_neighbours = tuple(() for _ in ploidies)
    if sample_count < NEIGHBOUR_MIN_SAMPLES:
        return NeighbourChoice(no_neighbours, 0, 0, {}, f"fewer than {NEIGHBOUR_MIN_SAMPLES} samples")
    first_stage = dataclasses.replace(settings, prior=PRIOR_HWE, fit_allele_bias=False, neighbours=None)
    shares_by_marker = [
        shares
        for shares in parallel_map(
            functools.partial(expected_alternate_shares, ploidies=ploidies, settings=first_stage),
            _evenly_spaced(depths_by_marker, NEIGHBOUR_MAX_MARKERS),
        )
        if shares is not None
    ]
    if len(shares_by_marker) < NEIGHBOUR_MIN_MARKERS:
        return NeighbourChoice(
            no_neighbours, 0, len(shares_by_marker), {}, f"fewer than {NEIGHBOUR_MIN_MARKERS} markers with reads"
        )

    def rows_by_sample(markers: Sequence[Sequence[float]]) -> List[array]:
        """Turn ``[marker][sample]`` into one compact row of shares per sample."""
        return [array("d", (marker_shares[sample] for marker_shares in markers)) for sample in range(sample_count)]

    halves = (shares_by_marker[0::2], shares_by_marker[1::2])
    judged_halves = [_evenly_spaced(half, NEIGHBOUR_JUDGED_MARKERS) for half in halves]
    neighbour_counts = [count for count in NEIGHBOUR_COUNTS if count < sample_count]
    block_size = -(-sample_count // NEIGHBOUR_BLOCKS)
    blocks = [range(start, min(start + block_size, sample_count)) for start in range(0, sample_count, block_size)]
    candidates = [
        candidate
        for block_candidates in parallel_map(
            functools.partial(
                _neighbour_candidates,
                shares_by_half=[rows_by_sample(half) for half in halves],
                judged_by_half=[rows_by_sample(half) for half in judged_halves],
                judged_means=[[sum(marker_shares) / sample_count for marker_shares in half] for half in judged_halves],
                neighbour_counts=neighbour_counts,
            ),
            blocks,
        )
        for candidate in block_candidates
    ]
    judged_cells = sample_count * sum(len(half) for half in judged_halves)
    error_by_count = {
        count: sum(errors[index] for _, errors in candidates) / judged_cells
        for index, count in enumerate(neighbour_counts)
    }
    best = min(neighbour_counts, key=lambda count: (error_by_count[count], count))
    reason = None if best else "neighbours predict the samples' allele shares no better than the whole table does"
    return NeighbourChoice(
        tuple(tuple(nearest[:best]) for nearest, _ in candidates), best, len(shares_by_marker), error_by_count, reason
    )


# ---------------------------------------------------------------------------
# Input: the MADC count table and the ploidy table
# ---------------------------------------------------------------------------

VCF_TARGET = "target"
VCF_ALL = "all"
VCF_HAPLOTYPE = "haplotype"
VCF_TYPES = (VCF_TARGET, VCF_ALL, VCF_HAPLOTYPE)
VCF_FILE_SUFFIX = {VCF_TARGET: ".targets.vcf", VCF_ALL: ".all.vcf", VCF_HAPLOTYPE: ".haplotypes.vcf"}

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


def _without_negative_cells(rows: Sequence[MadcRow]) -> Tuple[List[List[int]], int]:
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
    background_rate_estimated: bool = False


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
    """Phred-scaled probability that the called dosage is wrong, capped at 99.

    The VCF 4.3 definition of GQ, applied to the dosage posterior.  polyRAD's
    VCF export writes no such field.
    """
    chance_wrong = sum(probability for dosage, probability in enumerate(posterior) if dosage != called_dosage)
    if chance_wrong <= 0.0:
        return MAX_GENOTYPE_QUALITY
    return min(MAX_GENOTYPE_QUALITY, int(round(-10.0 * math.log10(chance_wrong))))


def phred_scaled_likelihoods(log_likelihoods: Sequence[float]) -> List[int]:
    """VCF PL values (VCF 4.3): each dosage's likelihood relative to the most likely one, Phred-scaled."""
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
) -> Tuple[str, Counter, Optional["BiasRow"]]:
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
    share of its reads was excluded.  Returns the line, the tallies for the
    QC report, and a ``BiasRow`` when an allelic bias was fitted (None
    otherwise).
    """
    result = call_record(depths, ploidies, options.settings)
    tallies: Counter = Counter(records=1)
    if result.has_reads and not result.converged:
        tallies["records_not_converged"] += 1
    used_reads = sum(sum(sample_depths) for sample_depths in depths)
    excluded_reads = sum(excluded_depths) if excluded_depths is not None else 0
    tallies["reads_used"] += used_reads
    tallies["reads_excluded"] += excluded_reads
    if excluded_reads and excluded_reads >= NOTABLE_EXCLUDED_READ_SHARE * (used_reads + excluded_reads):
        tallies["records_with_many_excluded_reads"] += 1

    # What follows GT, AD, DP and XD in a sample's column, and the confidence a genotype needs.
    if two_allele_model:
        value_keys = ["GQ", "DS", "GP", "PL"]
        confidence_threshold = options.min_genotype_quality
    else:
        value_keys = ["DS", "MPC"]
        confidence_threshold = options.min_marginal_posterior
    sample_columns = []
    samples_called = 0
    for sample, (sample_depths, ploidy) in enumerate(zip(depths, ploidies)):
        depth = sum(sample_depths)
        excluded = excluded_depths[sample] if excluded_depths is not None else 0
        dosages = None
        if not result.has_reads or depth == 0:
            missing_reason = "no_reads"
            values = ["."] * len(value_keys)
        else:
            if two_allele_model:
                dosages, confidence, values = _two_allele_genotype(
                    result.posteriors[sample], result.log_likelihoods[sample], ploidy
                )
            else:
                dosages, confidence, values = _multiallelic_genotype(result.posteriors[sample], ploidy, tallies)
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
    if result.has_reads and len(alleles) > 1:
        info.append("AF=" + ",".join(_decimal(frequency) for frequency in result.allele_frequencies[1:]))
    bias_row = None
    if result.bias_fit is not None:
        # The fitted bias, whether or not the genotypes were called with it.  Written for the
        # reference allele, which is what updog reports when it is given reference read
        # counts; the search itself works on the alternate allele.
        fit = result.bias_fit
        reference_bias = math.exp(-fit.log_bias)
        info.append("BIAS=" + _decimal(reference_bias, 3))
        if not fit.used:
            info.append("BIASUNUSED")
        bias_row = BiasRow(
            position, record_id, reference_bias, fit.objective, fit.gain, fit.rejection,
            math.exp(-fit.nearest_log_start), fit.fits, result.converged, fit.at_limit,
        )
    format_keys = ["GT", "AD", "DP"] + (["XD"] if excluded_depths is not None else []) + value_keys
    line = "\t".join(
        [marker, str(position), record_id, alleles[0], ",".join(alleles[1:]) or ".", ".", ".", ";".join(info),
         ":".join(format_keys)]
        + sample_columns
    )
    return line, tallies, bias_row


class BiasRow(NamedTuple):
    """One line of the bias table: a record's fitted allelic bias and how the fit went.

    ``bias`` and ``best_start`` are for the reference allele, as in the VCF's
    BIAS field.  ``objective`` is the log-likelihood of the record's read
    counts at the fitted bias plus the log of the prior on the bias, ``gain``
    what the fitted bias added to it over no bias, and ``rejection`` why the
    genotypes were called without the bias, or None when they were called
    with it (see ``BiasFit``).
    """

    position: int
    record_id: str
    bias: float
    objective: float
    gain: float
    rejection: Optional[str]
    best_start: float
    fits: int
    converged: bool
    at_limit: bool


@dataclasses.dataclass
class MarkerOutput:
    """What ``process_marker`` produced for one marker.

    ``lines`` holds the VCF lines per VCF type, ``tallies`` the QC counts
    keyed by ``(vcf type or "markers", name)``, and ``haplotype_alleles`` one
    ``(allele index, MADC allele IDs, sequence, total reads)`` entry per
    allele of the haplotype record.  ``bias_rows`` has one entry per record
    with a fitted allelic bias.  ``tag_length`` is 0 for a marker without a
    Ref tag, which gets no record.
    """

    clone_id: str
    tag_length: int = 0
    lines: Dict[str, List[str]] = dataclasses.field(default_factory=dict)
    tallies: Counter = dataclasses.field(default_factory=Counter)
    haplotype_alleles: List[Tuple[int, str, str, int]] = dataclasses.field(default_factory=list)
    bias_rows: List[BiasRow] = dataclasses.field(default_factory=list)


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
    counts_by_row, samples_cleared = _without_negative_cells(rows)
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
            line, tallies, bias_row = format_record(
                task.clone_id, target.position, task.clone_id, [target.ref_allele, target.alt_allele],
                depths, other_depths, ploidies, options, two_allele_model=True, is_target=True,
            )
            for vcf_type in site_types:
                add(vcf_type, line, tallies)
            if bias_row is not None:
                output.bias_rows.append(bias_row)

    if VCF_ALL in options.vcf_types:
        for line, tallies, bias_row in _off_target_records(
            task.clone_id, rows, counts_by_row, ref_tag, alt_tag, target, ploidies, options
        ):
            add(VCF_ALL, line, tallies)
            if bias_row is not None:
                output.bias_rows.append(bias_row)
        output.bias_rows.sort(key=lambda row: row.position)
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
    """Yield ``format_record``'s result for each substitution site away from the target.

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
        yield format_record(
            clone_id, position, f"{clone_id}_{position}", alleles,
            [[depth_by_base[base][sample] for base in alleles] for sample in range(sample_count)],
            _summed_counts(counts_by_row, rows_not_showing, sample_count),
            ploidies, options, two_allele_model=len(alleles) == 2,
        )


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
    # Haplotype records are called without an allelic bias, including those that happen to
    # have two alleles, so that every record of the file is called the same way.
    without_bias = dataclasses.replace(options, settings=options.settings.without_bias())
    line, tallies, _ = format_record(
        clone_id, 1, clone_id, alleles,
        [[depth_by_sequence[sequence][sample] for sequence in alleles] for sample in range(sample_count)],
        None, ploidies, without_bias, two_allele_model=False,
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


def model_description(
    settings: ModelSettings, overdispersion_estimated: bool, background_rate_estimated: bool
) -> str:
    """One sentence stating the genotype model and its inheritance assumption, for the VCF header."""
    source = ", estimated from the read counts" if overdispersion_estimated else ""
    background_source = ", estimated from the read counts" if background_rate_estimated else ""
    if settings.prior == PRIOR_HWE:
        mating = (
            f"a selfing rate of {settings.selfing_rate:g}" if settings.selfing_rate > 0 else "random mating"
        )
        prior = (
            "Hardy-Weinberg prior with one allele frequency per record estimated over all samples together; "
            f"it assumes {mating} and that each chromosome copy is drawn independently "
            "(polysomic inheritance in polyploids)"
        )
    elif settings.prior == PRIOR_STRUCTURE:
        neighbour_count = len(settings.neighbours[0]) if settings.neighbours else 0
        if neighbour_count:
            prior = (
                "Hardy-Weinberg prior at each sample's own allele frequency, which is the mean over the "
                f"{neighbour_count} samples most like it across the table's markers and the record's frequency "
                "counted as one more sample; it assumes random mating among like samples and that each chromosome "
                "copy is drawn independently (polysomic inheritance in polyploids)"
            )
        else:
            prior = (
                "Hardy-Weinberg prior with one allele frequency per record estimated over all samples together "
                "(the structure prior was asked for and found no like samples to use); it assumes random mating "
                "and that each chromosome copy is drawn independently (polysomic inheritance in polyploids)"
            )
    elif settings.prior == PRIOR_NORM:
        prior = (
            "normal-shaped prior over dosages (updog's norm; Gerard and Ferrao 2020), with a mean and standard "
            f"deviation fitted per record and ploidy; a ploidy with fewer than {FITTED_PRIOR_MIN_SAMPLES} samples "
            "uses Hardy-Weinberg proportions"
        )
    else:
        prior = "uniform prior over dosages"
    if settings.fit_allele_bias:
        prior += (
            "; allelic bias fitted for each record with two alleles (updog's bias; Gerard et al. 2018), with a "
            f"normal prior of standard deviation {BIAS_LOG_PRIOR_SD:g} on its log"
        )
        if settings.min_bias_gain > 0:
            prior += f", and used where it raises the log-likelihood by {settings.min_bias_gain:g} or more"
        if settings.reject_bias_at_limit:
            prior += " and does not end at a limit of the search"
    without_bias = ""
    if settings.read_model_without_bias is not None:
        overdispersion, background_rate = settings.read_model_without_bias
        without_bias = (
            f"; for records called without a bias, overdispersion {overdispersion:g} and background-read rate "
            f"{background_rate:g}"
        )
    return (
        f"beta-binomial read-count likelihood (overdispersion {settings.overdispersion:g}{source}, contamination rate "
        f"{settings.contamination_rate:g}, background-read rate {settings.background_rate:g}{background_source}"
        f"{without_bias}) adapted from polyRAD; {prior}"
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
        "##madcatsModel="
        + model_description(options.settings, options.overdispersion_estimated, options.background_rate_estimated),
        f"##madcatsFilters=min_dp={options.min_depth}, min_gq={options.min_genotype_quality}, "
        f"min_marginal_posterior={options.min_marginal_posterior:g}, "
        f"max_other_fraction={other_fraction}"
        "; a genotype that fails is written as missing with its AD and DP kept",
    ]
    lines += [f"##contig=<ID={marker},length={length}>" for marker, length in contigs]
    if vcf_type != VCF_HAPLOTYPE:
        lines.append('##INFO=<ID=TARGET,Number=0,Type=Flag,Description="The variant the marker was designed for">')
    if vcf_type != VCF_HAPLOTYPE and options.settings.fit_allele_bias:
        lines.append(
            '##INFO=<ID=BIAS,Number=1,Type=Float,Description="Fitted allelic bias of the reference allele, as updog '
            'reports it for reference read counts: below 1 the reference allele is over-represented among the reads, '
            'above 1 under-represented (two-allele records)">'
        )
        if options.settings.min_bias_gain > 0 or options.settings.reject_bias_at_limit:
            lines.append(
                '##INFO=<ID=BIASUNUSED,Number=0,Type=Flag,Description="The genotypes of this record were called '
                'without its fitted bias, which raised the log-likelihood of the read counts by less than '
                f'{options.settings.min_bias_gain:g} or ended at a limit of the search (--allele-bias auto)">'
            )
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


AUTO = "auto"


def _number_or_auto(parse):
    """An argparse type that accepts ``auto`` or whatever ``parse`` accepts."""
    def parse_value(value: str):
        return AUTO if value.strip().lower() == AUTO else parse(value)

    return parse_value


def target_depths_by_marker(tasks: Sequence[MarkerTask], sample_count: int) -> List[List[List[int]]]:
    """Reference-type and alternate-type read counts of every marker that has an alternate row.

    Returns ``[marker][sample][allele]``.  These are the counts the target
    records are called from; the read-model estimate and the search for like
    samples use them because every marker has them, whatever VCF types were
    asked for.
    """
    depths_by_marker = []
    for task in tasks:
        rows = [row for row in task.rows if row.sequence]
        if not any(row.category in ALTERNATE_CATEGORIES for row in rows):
            continue
        counts_by_row, _ = _without_negative_cells(rows)
        depths_by_marker.append(target_allele_depths(rows, counts_by_row, sample_count)[0])
    return depths_by_marker


@dataclasses.dataclass
class ReadModelChoice:
    """The overdispersion and background-read rate a run uses, and where each came from.

    ``overdispersion_estimated`` and ``background_rate_estimated`` are False
    when the value was given on the command line, and also when the estimate
    could not be used and the fallback value took its place.  The two reports
    go into the QC JSON.

    ``without_bias`` is set when the settings can leave a record without the
    bias fitted to it and the read model is estimated: the read model for
    such records, estimated without a bias.  The values of this object are
    then for the records called with their bias.
    """

    overdispersion: float
    background_rate: float
    overdispersion_estimated: bool
    background_rate_estimated: bool
    overdispersion_report: Dict[str, object]
    background_report: Dict[str, object]
    without_bias: Optional["ReadModelChoice"] = None


def _read_model_from_estimate(
    estimate: ReadModelEstimate, requested_overdispersion, requested_background_rate, used_for: Optional[str] = ""
) -> ReadModelChoice:
    """The read model a run takes from an estimate, with its reports for the QC JSON.

    A value given on the command line is kept.  A value that could not be
    estimated is replaced by its fallback: both when the table is too small,
    and the overdispersion also when its estimate is the lowest value tried.
    What is used is said on stderr, followed by ``used_for`` when the run has
    two read models, and nothing is said when ``used_for`` is None.
    """
    def say(message: str) -> None:
        if used_for is not None:
            print(message + (f", {used_for}" if used_for else ""), file=sys.stderr)

    overdispersion_report: Dict[str, object] = {"requested": requested_overdispersion}
    background_report: Dict[str, object] = {"requested": requested_background_rate}
    overdispersion_estimate = estimate.overdispersion_estimate
    too_small = overdispersion_estimate.value is None

    background_rate = requested_background_rate
    background_rate_estimated = requested_background_rate == AUTO and not too_small
    if requested_background_rate == AUTO:
        background_report["estimate"] = estimate.background_rate if background_rate_estimated else None
        background_report["by_rate"] = (
            {
                format(rate, "g"): {"overdispersion": value, "log_likelihood": round(total, 3)}
                for rate, (value, total) in estimate.by_background_rate.items()
            }
            if background_rate_estimated
            else {}
        )
        background_rate = estimate.background_rate if background_rate_estimated else FALLBACK_BACKGROUND_RATE
        at_end = background_rate_estimated and background_rate in (
            min(BACKGROUND_RATE_GRID), max(BACKGROUND_RATE_GRID)
        )
        say(
            f"{'Estimated' if background_rate_estimated else 'Using'} background-read rate: {background_rate:g}"
            + (" (the end of the range tried)" if at_end else "")
        )
    background_report["used"] = background_rate

    overdispersion = requested_overdispersion
    overdispersion_estimated = False
    if requested_overdispersion == AUTO:
        overdispersion_report.update(
            estimate=overdispersion_estimate.value,
            markers_used=overdispersion_estimate.markers_used,
            cells_with_reads=overdispersion_estimate.cells_with_reads,
            values_by_round=overdispersion_estimate.values_by_round,
            at_limit=overdispersion_estimate.at_limit,
            log_likelihood_by_value={
                format(value, "g"): round(total, 3)
                for value, total in overdispersion_estimate.log_likelihood_by_value.items()
            },
        )
        if too_small:
            overdispersion = FALLBACK_OVERDISPERSION
            say(f"Using overdispersion: {overdispersion:g}")
        elif overdispersion_estimate.at_limit == "lowest":
            overdispersion = FALLBACK_OVERDISPERSION
            say(
                f"warning: the overdispersion estimate is the lowest value tried ({overdispersion_estimate.value:g}), "
                "as it is when no sample is heterozygous (inbred lines) and the reads show no scatter to measure. "
                f"Using {overdispersion:g}; set --overdispersion to a number to choose it yourself."
            )
        else:
            overdispersion = float(overdispersion_estimate.value)
            overdispersion_estimated = True
            limit_note = ""
            if overdispersion_estimate.at_limit == "highest":
                limit_note = (
                    "; this is the highest value tried, so the counts scatter no more than random sampling of reads"
                )
            say(
                f"Estimated overdispersion: {overdispersion:g} "
                f"(from {overdispersion_estimate.markers_used:,} markers{limit_note})"
            )
    overdispersion_report["used"] = overdispersion
    return ReadModelChoice(
        overdispersion, background_rate, overdispersion_estimated, background_rate_estimated, overdispersion_report,
        background_report,
    )


def choose_read_model(
    depths_by_marker: Sequence[Sequence[Sequence[int]]],
    ploidies: Sequence[int],
    settings: ModelSettings,
    requested_overdispersion,
    requested_background_rate,
    parallel_map=map,
) -> ReadModelChoice:
    """Settle the overdispersion and the background-read rate of a run.

    ``requested_overdispersion`` and ``requested_background_rate`` are
    numbers or ``auto``.  Whatever is ``auto`` is estimated from each
    marker's reference-type and alternate-type read counts,
    ``depths_by_marker`` (``estimate_read_model``).

    The estimate is made with the free-weights prior whatever prior the run
    uses, because the estimate is only as good as the genotype proportions
    it assumes: with Hardy-Weinberg proportions it was badly off for
    samples from two subpopulations, and with those or the ``norm`` prior,
    neither of which can describe an excess of homozygotes, it failed for
    inbred samples.  An overdispersion that cannot be estimated falls back
    to ``FALLBACK_OVERDISPERSION``: when the table is too small, and when the
    best value is the lowest one tried, which is what a table without
    heterozygous samples and without stray reads gives, since it has no
    scatter to measure.

    With an allelic bias there are two estimates.  The first is made without
    a bias, as for a run that fits none.  The second fits every record's bias
    and keeps it, and its search starts from the first: started from
    polyRAD's value it stopped at an overdispersion of 10 on a simulated
    hexaploid table whose counts are far more probable at 20, with or without
    a bias, because a bias and dosage weights fitted at one overdispersion
    make that value look best.  The second estimate is the read model of the
    records called with their bias.  Where the settings can leave a record
    without its bias (``min_bias_gain``, ``reject_bias_at_limit``), the
    first is the read model of the records called without one, so that such
    a record is called exactly as a run without a bias calls it
    (docs/DEVELOPMENT.md).
    """
    if AUTO not in (requested_overdispersion, requested_background_rate):
        return ReadModelChoice(
            requested_overdispersion, requested_background_rate, False, False,
            {"requested": requested_overdispersion, "used": requested_overdispersion},
            {"requested": requested_background_rate, "used": requested_background_rate},
        )

    background_rates = BACKGROUND_RATE_GRID if requested_background_rate == AUTO else (requested_background_rate,)
    overdispersions = OVERDISPERSION_GRID if requested_overdispersion == AUTO else (requested_overdispersion,)
    estimate_settings = dataclasses.replace(
        settings, prior=PRIOR_FREE, selfing_rate=0.0, min_bias_gain=0.0, reject_bias_at_limit=False,
        read_model_without_bias=None,
    )
    estimate_without_bias = estimate_read_model(
        depths_by_marker, ploidies, dataclasses.replace(estimate_settings, fit_allele_bias=False), background_rates,
        overdispersions, parallel_map,
    )
    baseline = estimate_without_bias.overdispersion_estimate
    if baseline.value is None:
        print(
            f"note: too few read counts to estimate the read model ({baseline.cells_with_reads:,} "
            f"sample-by-marker cells with reads; {OVERDISPERSION_MIN_CELLS} needed).",
            file=sys.stderr,
        )
    if not settings.fit_allele_bias:
        return _read_model_from_estimate(estimate_without_bias, requested_overdispersion, requested_background_rate)

    some_records_without_bias = settings.min_bias_gain > 0 or settings.reject_bias_at_limit
    without_bias = _read_model_from_estimate(
        estimate_without_bias, requested_overdispersion, requested_background_rate,
        used_for="for records called without a bias" if some_records_without_bias else None,
    )
    starts = {}
    if baseline.value is not None:
        starts = dict(
            start_background_rate=estimate_without_bias.background_rate, start_overdispersion=baseline.value
        )
    with_bias = _read_model_from_estimate(
        estimate_read_model(
            depths_by_marker, ploidies, estimate_settings, background_rates, overdispersions, parallel_map, **starts
        ),
        requested_overdispersion, requested_background_rate,
        used_for="for records called with their bias" if some_records_without_bias else "",
    )
    # The estimate without a bias is reported beside the one with it; its "used" is null when
    # no record is called with it.
    for report, baseline_report in (
        (with_bias.overdispersion_report, without_bias.overdispersion_report),
        (with_bias.background_report, without_bias.background_report),
    ):
        if report["requested"] == AUTO:
            report["without_bias"] = {
                key: value if key != "used" or some_records_without_bias else None
                for key, value in baseline_report.items() if key != "requested"
            }
    if some_records_without_bias:
        with_bias.without_bias = without_bias
    return with_bias


def _background_rate(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed < 0.5:
        raise argparse.ArgumentTypeError("must be at least 0 and below 0.5")
    return parsed


def write_vcfs(
    prefix: str,
    vcf_types: Sequence[str],
    sample_names: Sequence[str],
    contigs: Sequence[Tuple[str, int]],
    options: CallOptions,
    outputs: Sequence[MarkerOutput],
) -> Dict[str, Path]:
    """Write one VCF per requested type; return each type's path."""
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
    return written


def write_haplotype_allele_table(prefix: str, outputs: Sequence[MarkerOutput]) -> Path:
    """Write the table that ties each allele number of the haplotype VCF to the MADC rows it was counted from."""
    alleles_path = Path(prefix + ".haplotypes.alleles.csv")
    with alleles_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["CloneID", "AlleleIndex", "AlleleIDs", "AlleleSequence", "TotalReads"])
        for output in outputs:
            for allele_index, allele_ids, sequence, total_reads in output.haplotype_alleles:
                writer.writerow([output.clone_id, allele_index, allele_ids, sequence, total_reads])
    return alleles_path


def write_bias_table(prefix: str, outputs: Sequence[MarkerOutput]) -> Tuple[Optional[Path], Dict[str, object]]:
    """Write the fitted allelic biases, one row per record; return the path and a summary for the QC JSON.

    Nothing is written, and the summary is empty, when no bias was fitted.
    The summary counts the records called without their fitted bias by
    reason, and names those whose bias was rejected for ending at a limit of
    the search.
    """
    bias_rows = [(output.clone_id, row) for output in outputs for row in output.bias_rows]
    if not bias_rows:
        return None, {}
    bias_path = Path(prefix + ".bias.csv")
    with bias_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["CloneID", "Position", "RecordID", "Bias", "LogLikelihood", "Gain", "Used", "Rejected", "BestStart",
             "Fits", "Converged", "AtLimit"]
        )
        for clone_id, row in bias_rows:
            writer.writerow(
                [clone_id, row.position, row.record_id, _decimal(row.bias), _decimal(row.objective, 3),
                 _decimal(row.gain, 3), int(row.rejection is None), row.rejection or "", _decimal(row.best_start),
                 row.fits, int(row.converged), int(row.at_limit)]
            )
    biases = sorted(row.bias for _, row in bias_rows)
    summary = {
        "records_fitted": len(biases),
        "records_called_with_the_bias": sum(1 for _, row in bias_rows if row.rejection is None),
        "records_called_without_it": {
            reason: sum(1 for _, row in bias_rows if row.rejection == reason)
            for reason in (BIAS_REJECTED_FOR_GAIN, BIAS_REJECTED_AT_LIMIT)
        },
        "rejected_at_search_limit": [
            f"{clone_id}:{row.position}" for clone_id, row in bias_rows if row.rejection == BIAS_REJECTED_AT_LIMIT
        ],
        "median": round(biases[len(biases) // 2], 4),
        "outside_0.8_to_1.25": sum(1 for bias in biases if not 0.8 <= bias <= 1.25),
        "outside_0.67_to_1.5": sum(1 for bias in biases if not 2 / 3 <= bias <= 1.5),
        "at_search_limit": sum(1 for _, row in bias_rows if row.at_limit),
        "not_converged": sum(1 for _, row in bias_rows if not row.converged),
    }
    return bias_path, summary


def _warn_about_biases_at_limit(bias_summary: Dict[str, object]) -> None:
    """Say on stderr when fitted biases ended at a limit of the search (``write_bias_table``'s summary)."""
    fitted = bias_summary.get("records_fitted", 0)
    at_limit = bias_summary.get("at_search_limit", 0)
    if not at_limit:
        return
    rejected = bias_summary["records_called_without_it"][BIAS_REJECTED_AT_LIMIT]
    not_used = (
        f" The {rejected:,} that gained enough to be used were not used: those records were called without a "
        "bias, are marked in the bias table and are listed in the QC JSON."
        if rejected else ""
    )
    if at_limit >= BIAS_AT_LIMIT_WARNING_SHARE * fitted:
        print(
            f"warning: {at_limit:,} of {fitted:,} fitted biases ended at a limit of the search.{not_used} Samples "
            "with few heterozygotes give this, as inbred lines do, and a bias has nothing to correct in them: call "
            "them with --allele-bias off, and with --selfing-rate if they are inbred",
            file=sys.stderr,
        )
    elif rejected:
        print(
            f"warning: the fitted bias of {rejected:,} of {fitted:,} records ended at a limit of the search and "
            "was not used; those records were called without a bias. They are marked in the bias table and "
            "listed in the QC JSON; --allele-bias fit uses such a bias",
            file=sys.stderr,
        )


def run_call(args: argparse.Namespace) -> int:
    """Call genotypes from a MADC and write the requested VCFs (see the module docstring)."""
    madc = Path(args.madc)
    if not madc.is_file():
        raise MadcError(f"MADC not found: {madc}")
    started = time.time()
    vcf_types = tuple(vcf_type for vcf_type in VCF_TYPES if vcf_type in args.vcf)
    if args.selfing_rate > 0 and args.prior != PRIOR_HWE:
        raise MadcError("--selfing-rate applies to --prior hwe only")
    if args.allele_bias != "off" and args.prior == PRIOR_STRUCTURE:
        raise MadcError(f"--allele-bias {args.allele_bias} cannot be combined with --prior structure")

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
        ploidy for ploidy, count in Counter(ploidies).items() if count < FITTED_PRIOR_MIN_SAMPLES
    )
    if args.prior == PRIOR_NORM and small_ploidy_groups:
        print(
            f"note: fewer than {FITTED_PRIOR_MIN_SAMPLES} samples have ploidy "
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

    model = ModelSettings(
        contamination_rate=args.contamination_rate,
        prior=args.prior,
        selfing_rate=args.selfing_rate,
        fit_allele_bias=args.allele_bias != "off",
        min_bias_gain=BIAS_AUTO_MIN_GAIN if args.allele_bias == "auto" else 0.0,
        reject_bias_at_limit=args.allele_bias == "auto",
        tolerance=args.tol,
        max_iterations=args.max_iterations,
    )
    # One pool of worker processes serves the read-model estimate, the search for like
    # samples and the calling.
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) if args.jobs > 1 else None

    def parallel_map(function, items):
        if executor is None or len(items) < 2:
            return list(map(function, items))
        return list(executor.map(function, items, chunksize=max(1, len(items) // (args.jobs * 8))))

    try:
        # The read-model estimate and the search for like samples both work on these counts.
        needs_target_depths = AUTO in (args.overdispersion, args.background_rate) or args.prior == PRIOR_STRUCTURE
        target_depths = target_depths_by_marker(tasks, len(ploidies)) if needs_target_depths else []
        read_model = choose_read_model(
            target_depths, ploidies, model, args.overdispersion, args.background_rate, parallel_map
        )
        settings = dataclasses.replace(
            model, overdispersion=read_model.overdispersion, background_rate=read_model.background_rate,
            read_model_without_bias=(
                (read_model.without_bias.overdispersion, read_model.without_bias.background_rate)
                if read_model.without_bias is not None else None
            ),
        )
        structure_report: Dict[str, object] = {}
        if args.prior == PRIOR_STRUCTURE:
            choice = find_neighbours(target_depths, ploidies, settings, parallel_map)
            settings = dataclasses.replace(settings, neighbours=choice.neighbours)
            structure_report = {
                "neighbours_per_sample": choice.neighbours_per_sample,
                "markers_compared": choice.markers_used,
                "prediction_error_by_neighbours": {
                    str(count): round(error, 6) for count, error in choice.error_by_count.items()
                },
                "reason_none_used": choice.reason,
            }
            if choice.neighbours_per_sample:
                print(
                    f"Structure prior: each sample's allele frequencies are taken from the "
                    f"{choice.neighbours_per_sample} samples most like it "
                    f"(compared on {choice.markers_used:,} markers)",
                    file=sys.stderr,
                )
            else:
                print(
                    f"note: --prior structure uses Hardy-Weinberg proportions for the whole table: {choice.reason}",
                    file=sys.stderr,
                )
        del target_depths
        options = CallOptions(
            settings=settings,
            vcf_types=vcf_types,
            min_depth=args.min_dp,
            min_genotype_quality=args.min_gq,
            min_marginal_posterior=args.min_marginal_posterior,
            max_other_fraction=args.max_other_fraction,
            gap_band=args.max_indel_length,
            overdispersion_estimated=read_model.overdispersion_estimated,
            background_rate_estimated=read_model.background_rate_estimated,
        )
        outputs = parallel_map(functools.partial(process_marker, ploidies=ploidies, options=options), tasks)
    finally:
        if executor is not None:
            executor.shutdown()

    contigs = [(output.clone_id, output.tag_length) for output in outputs if output.tag_length]
    tallies: Counter = Counter()
    for output in outputs:
        tallies.update(output.tallies)

    prefix = str(args.output_prefix)
    Path(prefix).parent.mkdir(parents=True, exist_ok=True)
    written = write_vcfs(prefix, vcf_types, sample_names, contigs, options, outputs)
    alleles_path = write_haplotype_allele_table(prefix, outputs) if VCF_HAPLOTYPE in vcf_types else None
    bias_path, bias_summary = write_bias_table(prefix, outputs)

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
            "allele_bias": args.allele_bias,
            "overdispersion": args.overdispersion,
            "contamination_rate": args.contamination_rate,
            "background_rate": args.background_rate,
            "tol": args.tol,
            "max_iterations": args.max_iterations,
            "min_dp": args.min_dp,
            "min_gq": args.min_gq,
            "min_marginal_posterior": args.min_marginal_posterior,
            "max_other_fraction": args.max_other_fraction,
            "max_indel_length": args.max_indel_length,
            "jobs": args.jobs,
        },
        "model": model_description(
            options.settings, options.overdispersion_estimated, options.background_rate_estimated
        ),
        "overdispersion": read_model.overdispersion_report,
        "background_rate": read_model.background_report,
        "structure_prior": structure_report,
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
        "allele_bias": bias_summary,
        "outputs": {
            "haplotype_alleles": str(alleles_path) if alleles_path else None,
            "bias_table": str(bias_path) if bias_path else None,
        },
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
    _warn_about_biases_at_limit(bias_summary)
    print(f"Wrote QC to {qc_path}", file=sys.stderr)
    return 0


def add_call_command(subparsers) -> None:
    calling = subparsers.add_parser(
        "call",
        help="call genotypes from a MADC and write target, all-variant or haplotype VCFs",
        description=(
            "Estimate each sample's genotype from the read counts in a MADC and write one VCF for each record "
            "type chosen with --vcf. Ploidy can differ between samples. The model is adapted from polyRAD "
            "(beta-binomial likelihood, Hardy-Weinberg prior) and updog (allelic bias, read-error term, norm "
            "prior). Coordinates are marker-local: CHROM is the marker and POS a position on its Ref tag."
        ),
    )

    files = calling.add_argument_group("input and output")
    files.add_argument(
        "--madc",
        required=True,
        help=(
            "MADC to call from, plain or gzipped: one written by generate or filter-paralogs, a vendor DArT "
            "MADC (with or without its starred preamble) or a HapApp fixed-allele-ID MADC"
        ),
    )
    files.add_argument(
        "--output-prefix",
        required=True,
        metavar="PREFIX",
        help=(
            "path prefix of the outputs: PREFIX.targets.vcf, PREFIX.all.vcf and PREFIX.haplotypes.vcf (with "
            "PREFIX.haplotypes.alleles.csv) as chosen by --vcf, PREFIX.bias.csv when an allelic bias is "
            "fitted, and PREFIX.qc.json"
        ),
    )
    files.add_argument(
        "--vcf",
        nargs="+",
        choices=VCF_TYPES,
        default=[VCF_TARGET],
        metavar="TYPE",
        help=(
            "record types to write, one VCF each; one or more of target (the variant each marker was designed "
            "for), all (the targets plus every other base substitution in the tags) and haplotype (whole tag "
            f"sequences as alleles) (default: {VCF_TARGET})"
        ),
    )
    files.add_argument("--qc-output", help="QC JSON (default: PREFIX.qc.json)")

    panel = calling.add_argument_group("samples and panel")
    panel.add_argument(
        "--ploidy", type=positive_int, default=DEFAULT_PLOIDY,
        help=f"ploidy of every sample not listed in --ploidy-file (default: {DEFAULT_PLOIDY})",
    )
    panel.add_argument(
        "--ploidy-file",
        help="two-column table of sample name and ploidy, tab- or comma-separated, with or without a header line",
    )
    panel.add_argument(
        "--reference-tags",
        help=(
            "reference-tag CSV of the panel (from design or infer-design); its tags and "
            "VariantPositionInTagSequence column identify each marker's target. Without it the MADC's Ref and "
            "Alt rows are used, and a marker whose two tags differ in more than one place gets no target record"
        ),
    )
    panel.add_argument(
        "--max-indel-length", type=nonnegative_int, default=GAP_BAND_DEFAULT,
        help=(
            "widest insertion or deletion considered when a row is compared with the Ref tag; use the value "
            f"given to generate (default: {GAP_BAND_DEFAULT})"
        ),
    )

    model = calling.add_argument_group("model")
    model.add_argument(
        "--prior",
        choices=PRIOR_CHOICES,
        default=PRIOR_HWE,
        help=(
            "genotype prior: uniform (every dosage equally probable), hwe (Hardy-Weinberg proportions at "
            "allele frequencies estimated from all samples; assumes random mating and polysomic inheritance), "
            "norm (a normal shape over the dosages, fitted per record and ploidy as in updog; does not assume "
            "random mating) or structure (Hardy-Weinberg proportions at each sample's own allele frequencies, "
            "estimated from the samples most like it; for subpopulations, clines and families) "
            f"(default: {PRIOR_HWE})"
        ),
    )
    model.add_argument(
        "--selfing-rate", type=_fraction, default=0.0,
        help=(
            "with --prior hwe, the share of individuals that come from self-fertilisation, from 0 (random "
            "mating) to 1 (fully inbred lines); it moves prior probability from the intermediate dosages to "
            "the two homozygous ones (default: 0)"
        ),
    )
    model.add_argument(
        "--allele-bias",
        choices=("off", "fit", "auto"),
        default="off",
        help=(
            "allelic bias, where one allele yields more reads per copy than the other: off (none assumed), fit "
            "(estimated for every record with two alleles, as in updog, and used) or auto (estimated the same "
            f"way, but used only for records where it raises the log-likelihood by {BIAS_AUTO_MIN_GAIN:g} or "
            "more and does not end at a limit of the search; the rest are called without it). fit and auto "
            "take several times as long and cannot be combined with --prior structure (default: off)"
        ),
    )
    model.add_argument(
        "--overdispersion", type=_number_or_auto(_positive_float), default=AUTO,
        help=(
            "beta-binomial overdispersion parameter, as in polyRAD: lower values allow read fractions to stray "
            "further from the expected fraction. 'auto' estimates it from the read counts by maximum "
            f"likelihood, and uses {FALLBACK_OVERDISPERSION:g} when it cannot; or give a number "
            f"(default: {AUTO})"
        ),
    )
    model.add_argument(
        "--background-rate", type=_number_or_auto(_background_rate), default=AUTO,
        help=(
            "probability that a read shows one allele although it came from another: sequencing and "
            "amplification errors, index hopping and contamination together. 'auto' estimates it from the read "
            f"counts along with the overdispersion, and uses {FALLBACK_BACKGROUND_RATE:g} when it cannot; or give "
            f"a number (default: {AUTO})"
        ),
    )
    model.add_argument(
        "--contamination-rate", type=proper_fraction, default=0.0,
        help=(
            "polyRAD's contamination term: the expected share of a sample's reads that come from other "
            "samples. --background-rate already covers such reads, so this is 0 unless set; polyRAD's own "
            f"default is {POLYRAD_CONTAMINATION_RATE:g} (default: 0)"
        ),
    )
    model.add_argument(
        "--tol", type=_positive_float, default=DEFAULT_TOLERANCE,
        help=(
            "with --prior hwe, norm or structure, stop iterating a record once its allele frequencies (and, "
            f"with norm, its prior probabilities) change by no more than this (default: {DEFAULT_TOLERANCE:g})"
        ),
    )
    model.add_argument(
        "--max-iterations", type=positive_int, default=DEFAULT_MAX_ITERATIONS,
        help=(
            "with --prior hwe, norm or structure, the most iterations a record may take; with structure, in "
            f"each of its two stages (default: {DEFAULT_MAX_ITERATIONS})"
        ),
    )

    filters = calling.add_argument_group("filters (a genotype that fails one is written as missing)")
    filters.add_argument(
        "--min-dp", type=nonnegative_int, default=DEFAULT_MIN_DEPTH,
        help=f"fewest reads a sample needs in a record (default: {DEFAULT_MIN_DEPTH})",
    )
    filters.add_argument(
        "--min-gq", type=nonnegative_int, default=DEFAULT_MIN_GENOTYPE_QUALITY,
        help=f"two-allele records: lowest genotype quality, GQ (default: {DEFAULT_MIN_GENOTYPE_QUALITY})",
    )
    filters.add_argument(
        "--min-marginal-posterior", type=_fraction, default=DEFAULT_MIN_MARGINAL_POSTERIOR,
        help=(
            "haplotype records and other records with more than two alleles: lowest MPC, the smallest "
            f"marginal posterior among the called dosages (default: {DEFAULT_MIN_MARGINAL_POSTERIOR:g})"
        ),
    )
    filters.add_argument(
        "--max-other-fraction", type=_fraction, default=None,
        help=(
            "target and all records: largest share of the sample's reads that may be left out of the call, "
            "XD (default: no limit)"
        ),
    )

    performance = calling.add_argument_group("performance")
    performance.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(),
        help="worker processes; markers are called independently (default: up to 4)",
    )
    calling.set_defaults(func=run_call)
