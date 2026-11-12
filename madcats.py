#!/usr/bin/env python3
"""MADCATS — Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing.

A platform-independent tool for recovering design-specified alleles and
discovering additional alleles from targeted sequencing reads, with MADC
and VCF-ready output.  The name also reflects the tool's origins as a
FASTQ-to-MADC workflow, where MADC denotes the Missing Allele Discovery
Counts format.

The implementation intentionally uses only the Python standard library.  The
one exception is optional: ``filter-paralogs``, in its own module
``madcats_paralogs.py``, imports the ``mappy`` aligner on first use
(``pip install madcats[genome]``).

Seven subcommands cover the read models the pipeline understands, an eighth
cleans a finished table, and a ninth calls genotypes from it:

- ``design`` reconstructs each marker's canonical Ref/Alt tags from vendor
  primer + probe-design files.  This is the primer-started, single-end model
  (DArTag-style): every read begins with a synthesized primer, so the assay
  itself fixes the tag and the design files describe it exactly.
- ``infer-design`` recovers the same panel definition from the probe design
  plus FASTQ evidence alone, for when the primer file is unavailable.  The
  tag must be a tag-length window of the expanded target that contains the
  variant, so inference selects among a closed candidate set — exact
  read-prefix votes first, a flagged fuzzy tier for diverged loci.
- ``merge-pairs`` collapses overlapping paired-end reads into single
  fragments before normalization (quality-picked overlap bases; engulfed
  read-through geometry discards adapter tails by construction).  Paired
  data whose mates cannot overlap (dual-anchored capture designs) is
  instead pooled — both mates fed as independent reads — an explicit
  per-run decision made from measured geometry, never auto-detected.
- ``normalize-reads`` adapts platforms whose reads are neither oriented nor
  primer-started (AgriSeq-style bidirectional amplicons): each read is
  located against the design targets by k-mer, oriented forward, and a
  fixed-length window containing the variant is excised into pipeline-ready
  FASTQs plus a matching window panel.  Unlike a primer-defined tag, the
  window is an analytical choice, not a physical one — so later runs of the
  same panel must freeze it with ``--reuse-tags``, keeping microhaplotype
  sequences comparable across runs (an identical sequence must mean
  identical biology, never a moved window).  Markers may share a physical
  locus (``LocusGroup``): tiles share one k-mer territory, reads route
  among member windows (``--group-routing``), and member windows are
  forced to excise distinct sequences.
- ``window-scan`` sweeps candidate tag lengths on a sample subset and
  reports the three selection curves — retention, cleanliness, and
  microhaplotype information — with a recommended length, so the
  window-size choice is made from the project's own reads instead of an
  external truth comparison.
- ``generate`` turns per-sample FASTQs plus a panel from any of the three
  panel routes above into the MADC table, via the stages listed below.
- ``validate`` compares a generated MADC against a reference MADC row-by-row
  and cell-by-cell, reporting precision/recall by allele category.
- ``filter-paralogs`` (optional; needs ``mappy``; in ``madcats_paralogs.py``)
  cleans a MADC from ``generate``, a vendor delivery or HapApp with a
  reference genome.  Each marker is placed by aligning its Ref tag, every
  haplotype is aligned, and rows whose best genome hit lies at another locus
  (paralogs) or ties between loci (ambiguous) are dropped or flagged.  It is
  meant for panels where co-amplified paralogous copies inflate
  heterozygosity when counts are summed per SNP; canonical rows are never
  dropped.
- ``call`` (in ``madcats_call.py``) estimates each sample's genotype from the
  read counts of a MADC and writes VCFs of the designed targets, of every
  base substitution in the tags, or of whole tag sequences as alleles.  The
  model is adapted from polyRAD and ploidy is set per sample.

Pipeline stages (the ``generate`` command; see ``run_generate``):

1. Panel definition.  Reconstruct each marker's canonical Ref/Alt tags either
   from the DArTag primer + probe-design files (``build_panel``) or from a
   platform-neutral reference-tag CSV (``load_reference_tags``).  The tags are
   the ground truth every later stage measures reads against, so this stage
   validates aggressively and fails loudly on design inconsistencies.
2. Per-sample read processing (``process_sample``).  Assign each read to one
   locus by its left-primer prefix (``PrimerIndex`` trie), keep the fixed-
   length tag, and accept it within the substitution-only identity envelope
   of the locus's Ref/Alt tags — or, when every candidate rejects it, under
   one gapped reading within the same edit budget, so undesigned indel
   haplotypes are counted rather than dropped (``select_marker``).  Tags
   carrying an N are credited to the one N-free sequence they fit or
   dropped; they never become rows.  Counting exact observed sequences —
   never consensus-corrected ones — is what makes the output auditable
   against the raw FASTQs.  Results are cached per sample
   keyed by panel + read-processing parameters, because this stage is the
   only expensive one and downstream thresholds do not invalidate it.
3. Partial-read salvage (``apply_partial_reads``).  On by default; credits a
   primer-assigned short read to the one already-learned full tag it is a
   prefix of, and only when it is longer than every prefix two of the locus's
   alleles share, so truncated reads can support but never invent alleles or
   move an allele ratio.
4. Evidence collection (``collect_evidence`` under an ``EvidenceRule``).
   A novel sequence earns an output row through reproducibility across
   samples, single-sample depth, or optional whole-tag quality.  The same
   rule drives stages 5 and 6, so protection, discovery, and emission can
   never disagree about which sequences are real.
5. Denoising (``finalize_counts``).  UNOISE-style abundance-skew absorption
   of likely sequencing errors into their parent sequences within each
   sample, then a pooled pass that folds an evidence-failing sequence into
   the protected parent the project-wide skew names — only when both carry
   the same target allele, so no move can change a sample's dosage.
   Evidence-passing and canonical sequences are protected from absorption;
   absorbed reads move to their parent row rather than disappearing, and
   the moves are reported in QC.
6. Emission (``write_madc``).  Every panel locus always emits its canonical
   Ref and Alt rows — even at zero reads — because downstream dosage tools
   need the biallelic axis present unconditionally.  Novel rows must pass
   the evidence rule; excluded reads are tallied, never silently dropped.
7. QC/provenance.  A JSON sidecar records inputs, parameters, panel
   fingerprints, per-sample QC, and a full read-accounting ledger from raw
   records through emitted rows, so any run can be audited and reproduced.

Reading map.  The file is laid out in the order a read travels; each
section opens with a banner comment:

0. Constants, errors and data types (``Locus``, ``Sample``, ``SampleResult``).
1. Sequence and file utilities (IUPAC matching, mismatch counts, CSV/FASTQ
   helpers).
2. Alignment engine: the Myers prescreen, the banded DP, and
   ``accepted_gapped_reading``/``row_reading`` — the one gapped reading a
   row may be given.
3. Panel definition: design-file reconstruction, reference-tag tables, and
   primer-free inference.
4. Platform intake: pair merging, window location and freezing, tiling and
   routing, window-scan.
5. Per-sample read assignment: ``PrimerIndex``, ``select_marker``,
   ``process_sample`` and its cache.
6. Evidence, denoising and partial-read salvage: ``EvidenceRule``,
   ``collect_evidence``, ``finalize_counts``, ``apply_partial_reads``.
7. Classification and annotation: ``classify_sequence``, ``variant_events``,
   ``marker_type_label``.
8. MADC input/output and QC: ``read_madc``, ``write_madc``, ``validate_madc``.
   The genome paralog filter is in ``madcats_paralogs.py``.
9. Commands and command line: one ``run_*`` per subcommand, ``build_parser``.

Vocabulary.  *tag*: the fixed-length sequence a read is measured as (intake
code also says *window*, the position of that sequence on a longer read);
*marker* and *locus*: the same panel entry (``Locus`` is the record);
*canonical*: the design's Ref and Alt tags; *category*: Ref, Alt, RefMatch,
AltMatch or Other, read at the target base.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import dataclasses
import datetime as dt
import functools
import gzip
import hashlib
import itertools
import json
import math
import os
import random
import re
import shlex
import statistics
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Mapping, MutableMapping, NamedTuple, Optional, Sequence, Set, Tuple


# ============================================================================
# 0. Constants, errors and data types
# ============================================================================


__version__ = "0.13.0"

# Increment only when the serialized per-sample cache or raw-read processing
# changes. Output-row filtering does not invalidate reusable FASTQ counts.
# 0.4.0: fuzzy primer matches resolve across all matched depths by tag
# identity instead of preferring the deepest primer, so caches built with
# --primer-mismatches 1 (or prefix-nested panels) must be regenerated.
CACHE_FORMAT_VERSION = "0.5.0"
TAG_LENGTH_DEFAULT = 81
QUALITY_THRESHOLDS = (10, 20, 30, 40)

IUPAC_BASES = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),
    "Y": frozenset("CT"),
    "S": frozenset("CG"),
    "W": frozenset("AT"),
    "K": frozenset("GT"),
    "M": frozenset("AC"),
    "B": frozenset("CGT"),
    "D": frozenset("AGT"),
    "H": frozenset("ACT"),
    "V": frozenset("ACG"),
    "N": frozenset("ACGT"),
}
DNA_COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
_CONCRETE = frozenset("ACGT")
ALLELE_TYPE_ORDER = {"Ref": 0, "Alt": 1, "RefMatch": 2, "AltMatch": 3, "Other": 4}
ALLELE_GROUP = {"Ref": "ref", "RefMatch": "ref", "Alt": "alt", "AltMatch": "alt", "Other": "other"}

MADC_METADATA_COLUMNS = [
    "AlleleID",
    "CloneID",
    "AlleleSequence",
    "ClusterConsensusSequence",
    "CallRate",
    "OneRatioRef",
    "OneRatioSnp",
    "FreqHomRef",
    "FreqHomSnp",
    "FreqHets",
    "PICRef",
    "PICSnp",
    "AvgPIC",
    "AvgCountRef",
    "AvgCountSnp",
    "RatioAvgCountRefAvgCountSnp",
]

MARKER_TYPE_COLUMN = "AvgPIC"   # the vendor statistics column the MarkerType label replaces
VARIANTS_COLUMN = "PICSnp"      # the vendor statistics column the Variants receipt replaces
ANNOTATION_COLUMNS = {VARIANTS_COLUMN: "Variants", MARKER_TYPE_COLUMN: "MarkerType"}
VARIANTS_SIDECAR_COLUMNS = ["CloneID", "AlleleID", "AlleleSequence", "Position", "Ref", "Alt", "Type"]

GAP_BAND_DEFAULT = 10        # --max-indel-length: the widest gap a gapped reading may use
GAP_MAX_SUBSTITUTIONS = 6    # composition guard: an indel haplotype is its canonical plus gaps and FEW substitutions
GAP_MAX_EVENTS = 6           # ... in few separate gap events; a co-amplified paralog needs many substitutions
GAP_VERIFIED_SUFFIX = 12     # bases that must follow the last gap before a gapped reading may re-read the target
PARTIAL_TAIL_DEFAULT = 9     # --min-partial-length defaults to --tag-length minus this: at most 9 unobserved tail bases
LEADING_ANCHOR_SEARCH = 40   # --leading-anchor: the anchor must end within this many bases beyond its own length
DENOISE_LOOKUP_THRESHOLD = 2 # denoiser: look up single-substitution neighbours once centroids exceed this many per tag base
GAP_OPEN_SCALE = 1 << 10     # annotation alignments weigh one edit this much more than one gap opening


class MadcError(RuntimeError):
    """A user-facing pipeline error."""


@dataclasses.dataclass(frozen=True)
class Locus:
    index: int
    marker: str
    left_primer: str
    right_primer: str
    ref_base: str
    alt_base: str
    variant_pos: int  # zero based in the tag
    ref_tag: str
    alt_tag: str
    # Markers sharing a locus_group are windows (tiles) of one physical
    # locus: normalize-reads locates reads to the group and routes them
    # among the member windows (--group-routing: one tile per read, or
    # every fully covered tile).  Empty = the marker is its own group
    # (the universal pre-0.9 case).
    locus_group: str = ""

    @property
    def group(self) -> str:
        return self.locus_group or self.marker


@dataclasses.dataclass(frozen=True)
class Sample:
    sample_id: str
    display_name: str
    fastq: str
    metadata: Mapping[str, str]


@dataclasses.dataclass
class SampleResult:
    sample: Sample
    counts: Dict[Tuple[int, str], int]
    quality_support: Dict[Tuple[int, str], Tuple[int, ...]]
    partial_counts: Dict[Tuple[int, str], int]
    partial_quality_support: Dict[Tuple[int, str], Tuple[int, ...]]
    qc: Dict[str, int]
    cache_hit: bool = False


class EvidenceSupport(NamedTuple):
    """One sequence's pooled evidence, in the order ``EvidenceRule.evaluate`` takes it."""

    global_count: int        # reads across the project
    clean_reads: int         # reads whose every tag base met --discovery-min-quality
    supporting_samples: int  # samples with at least --reproducible-min-reads
    max_sample_count: int    # the deepest single sample


class GappedReading(NamedTuple):
    """An accepted gapped reading of a row against one canonical tag.

    ``ops`` is the fewest-openings alignment (see gapped_alignment); the
    counts describe it; ``last_gap_end`` is the read index just after the
    last gap; ``ref_to_read`` maps aligned reference indices to read indices
    (reference bases inside a deletion, or past the aligned end, are absent).
    """

    edits: int
    substitutions: int
    gap_bases: int
    gap_events: int
    ops: Tuple[Tuple[str, int, int], ...]
    end_read: int
    end_ref: int
    last_gap_end: int
    ref_to_read: Dict[int, int]


@dataclasses.dataclass(frozen=True)
class AssignmentSettings:
    """Everything the per-sample read pass depends on — and therefore the cache key.

    ``cache_arguments`` is the single source of the fingerprint: a setting
    that changes what a read pass produces must be a field here, so it can
    never be passed to ``process_sample`` without also invalidating the cache.
    """

    tag_length: int
    max_offtarget_mismatches: int   # floor((1 - min_identity) * tag_length)
    primer_mismatches: int
    partial_min_length: int
    gap_band: int                   # --max-indel-length, or 0 with --gapped-rescue off
    leading_anchor: str = ""        # --leading-anchor: constant sequence preceding the amplicon in raw deliveries

    def cache_arguments(self) -> Dict[str, object]:
        return {
            "tag_length": self.tag_length,
            "max_offtarget_mismatches": self.max_offtarget_mismatches,
            "gap_band": self.gap_band,
            "primer_mismatches": self.primer_mismatches,
            "partial_min_length": self.partial_min_length,
            "leading_anchor": self.leading_anchor,
            "quality_thresholds": QUALITY_THRESHOLDS,
        }


@dataclasses.dataclass(frozen=True)
class DenoiseSettings:
    """The abundance-skew denoiser's knobs (``finalize_counts``)."""

    sample_distance: int = 1   # per-sample closest-parent absorption; 0 disables
    global_distance: int = 1   # pooled same-allele fold; 0 disables
    alpha: float = 2.0         # UNOISE-style skew exponent


@dataclasses.dataclass(frozen=True)
class OutputOptions:
    """How ``write_madc`` lays the table out; nothing here changes a count."""

    number_match_ids: bool = False
    preamble: Optional[Sequence[Sequence[str]]] = None
    annotations: bool = True
    variants_path: Optional[Path] = None


# ============================================================================
# 1. Sequence and file utilities
# ============================================================================


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(DNA_COMPLEMENT)[::-1]


def bases_compatible(a: str, b: str) -> bool:
    """True when two bases could be the same nucleotide under IUPAC codes."""
    return bool(IUPAC_BASES.get(a, frozenset(a)) & IUPAC_BASES.get(b, frozenset(b)))


def compatible_mismatches(a: str, b: str, stop_after: Optional[int] = None) -> int:
    """Count substitutions, treating IUPAC codes as compatible base sets.

    This is the design-side comparison (tag against tag); observed reads are
    compared with ``tag_mismatches``, where an observed N is a mismatch.
    """
    if len(a) != len(b):
        return max(len(a), len(b))
    if a == b:
        return 0
    mismatches = 0
    for left, right in zip(a, b):
        if not (IUPAC_BASES.get(left, frozenset(left)) & IUPAC_BASES.get(right, frozenset(right))):
            mismatches += 1
            if stop_after is not None and mismatches > stop_after:
                break
    return mismatches


def hamming_distance(a: str, b: str, stop_after: Optional[int] = None) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    mismatches = 0
    for left, right in zip(a, b):
        if left != right:
            mismatches += 1
            if stop_after is not None and mismatches > stop_after:
                break
    return mismatches


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, newline="") if path.name.lower().endswith(".gz") else path.open(mode, newline="")


def normalized_row(row: Mapping[str, str]) -> Dict[str, str]:
    return {(key or "").strip(): (value or "").strip() for key, value in row.items()}


def read_dict_rows(path: Path, delimiter: Optional[str] = None) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        if delimiter is None:
            sample = handle.read(8192)
            handle.seek(0)
            try:
                delimiter = csv.Sniffer().sniff(sample, delimiters=",\t").delimiter
            except csv.Error:
                delimiter = "\t" if "\t" in sample.partition("\n")[0] else ","
        return [normalized_row(row) for row in csv.DictReader(handle, delimiter=delimiter)]


def require_columns(rows: Sequence[Mapping[str, str]], columns: Sequence[str], path: Path) -> None:
    if not rows:
        raise MadcError(f"Input table is empty: {path}")
    missing = [column for column in columns if column not in rows[0]]
    if missing:
        raise MadcError(f"{path} is missing required column(s): {', '.join(missing)}")


def max_tandem_run(segment: str) -> int:
    """Longest tandem-repeat stretch in bp for unit sizes 1-3 (>= 2 full units).

    Homopolymers are the unit-1 case; DArT-adjacent failure analysis showed
    dinucleotide repeats such as (AC)6 are equally hostile to flow-based
    chemistries, and a variant allele that extends an adjacent repeat suffers
    allele-biased read dropout inside any window that spans the full repeat.
    """
    best = 1
    for period in (1, 2, 3):
        run = period
        for i in range(period, len(segment)):
            if segment[i] == segment[i - period]:
                run += 1
                if run >= 2 * period:
                    best = max(best, run)
            else:
                run = period
    return best


def tag_mismatches(tag: str, template: str, stop_after: Optional[int] = None) -> int:
    """Substitutions between an OBSERVED tag and a design template.

    Design-side IUPAC codes are compatible sets; the observed base must be a
    concrete member of the set.  An observed N (or any ambiguity code) is
    missing information and counts as a mismatch — ``compatible_mismatches``,
    the design-side comparison, treats it as compatible with everything,
    which let N-bearing reads sail through the identity gate.
    """
    if len(tag) != len(template):
        return max(len(tag), len(template))
    mismatches = 0
    for observed, expected in zip(tag, template):
        if observed != expected and (
            observed not in _CONCRETE or observed not in IUPAC_BASES.get(expected, frozenset(expected))
        ):
            mismatches += 1
            if stop_after is not None and mismatches > stop_after:
                break
    return mismatches


def _design_matches(observed: str, design: str) -> bool:
    """A concrete observed base matches a design base or lies in its IUPAC set; N never matches."""
    return observed == design or (observed in _CONCRETE and observed in IUPAC_BASES.get(design, frozenset(design)))


def strip_fastq_suffix(name: str) -> str:
    lower = name.lower()
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def natural_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def mean_nonzero(values: Sequence[int]) -> Optional[float]:
    nonzero = [value for value in values if value > 0]
    return statistics.mean(nonzero) if nonzero else None


def format_number(value: Optional[float]) -> str:
    if value is None:
        return "0"
    if abs(value - round(value)) < 1e-12:
        return str(int(round(value)))
    return f"{value:.6f}".rstrip("0").rstrip(".")


def pearson_from_sums(n: int, sx: float, sy: float, sxx: float, syy: float, sxy: float) -> Optional[float]:
    if n == 0:
        return None
    numerator = n * sxy - sx * sy
    denominator = math.sqrt(max(0.0, n * sxx - sx * sx) * max(0.0, n * syy - sy * sy))
    return numerator / denominator if denominator else None


# ============================================================================
# 2. Alignment engine (prescreen, banded DP, the accepted reading)
# ============================================================================


def _myers_prefix_distances(pattern: str, text: str, design_is_pattern: bool) -> List[int]:
    """Global edit distance of ``pattern`` against every prefix of ``text``.

    Hyyro's bit-vector formulation of Myers' algorithm, anchored at the
    start (the first DP row counts text characters).  Python integers are
    the bit-vectors, so each text character costs a handful of C-speed word
    operations regardless of pattern length.  ``result[j]`` is the edit
    distance between ``pattern`` and ``text[:j]``.  One side is the design
    template and may carry IUPAC codes; a match is the same relation
    ``tag_mismatches`` uses (the design code itself, or a concrete observed
    base inside its set; an observed N matches only a design N), so the
    bound agrees with the alignment and with canonical tags that share a
    code.
    """
    m = len(pattern)
    if m == 0:
        return list(range(len(text) + 1))
    peq: Dict[str, int] = {}
    if design_is_pattern:
        for i, design in enumerate(pattern):
            for observed in IUPAC_BASES.get(design, frozenset(design)) | {design}:
                peq[observed] = peq.get(observed, 0) | (1 << i)
    else:
        for i, observed in enumerate(pattern):
            peq[observed] = peq.get(observed, 0) | (1 << i)
    eq_cache: Dict[str, int] = {}

    def eq_for(design: str) -> int:
        # text is the design: every read position whose base lies in the code's set
        got = eq_cache.get(design)
        if got is None:
            got = 0
            for observed in IUPAC_BASES.get(design, frozenset(design)) | {design}:
                got |= peq.get(observed, 0)
            eq_cache[design] = got
        return got

    mask = (1 << m) - 1
    high = 1 << (m - 1)
    pv, mv, score = mask, 0, m
    out = [m]
    for base in text:
        eq = peq.get(base, 0) if design_is_pattern else eq_for(base)
        xv = eq | mv
        xh = (((eq & pv) + pv) ^ pv) | eq
        ph = mv | (~(xh | pv) & mask)
        mh = pv & xh
        if ph & high:
            score += 1
        elif mh & high:
            score -= 1
        ph = ((ph << 1) | 1) & mask
        mh = (mh << 1) & mask
        pv = mh | (~(xv | ph) & mask)
        mv = ph & xv
        out.append(score)
    return out


def gapped_prescreen(read: str, ref: str, band: int) -> int:
    """Lower bound on the start-anchored, end-free gapped edit distance.

    The end may overhang by up to ``band`` bases on either side (a read
    tail past the reference is a deletion's signature; a reference tail
    past the read an insertion's).  Unbanded, so never above the banded
    distance: a value over the budget proves the banded alignment is too.
    """
    n, m = len(read), len(ref)
    read_consumed = _myers_prefix_distances(read, ref, design_is_pattern=False)
    ref_consumed = _myers_prefix_distances(ref, read, design_is_pattern=True)
    return min(min(read_consumed[max(0, m - band):]), min(ref_consumed[max(0, n - band):]))


def _gapped_dp(read: str, ref: str, band: int):
    """Banded DP behind gapped_stats and gapped_alignment.

    Returns ``(rows, ref_sets, edits, end_read, end_ref)``: the filled table,
    the design-side IUPAC sets, and the cheapest end cell.  Ties go to the
    shorter aligned read, then the shorter aligned reference, so every
    traceback over the same table follows the same path.
    """
    n, m = len(read), len(ref)
    INF = 1 << 30
    # Design-side IUPAC codes are compatible sets (the code itself included, so
    # canonical tags sharing a code agree); an observed N matches only a design N.
    ref_sets = [IUPAC_BASES.get(design, frozenset(design)) | {design} for design in ref]
    rows = [[INF] * (m + 1) for _ in range(n + 1)]
    for j in range(min(m, band) + 1):
        rows[0][j] = j
    for i in range(1, n + 1):
        lo, hi = max(0, i - band), min(m, i + band)
        current, previous = rows[i], rows[i - 1]
        if lo == 0:
            current[0] = i
        base = read[i - 1]
        for j in range(max(1, lo), hi + 1):
            cost = previous[j - 1] + (base not in ref_sets[j - 1])
            up = previous[j] + 1
            if up < cost:
                cost = up
            left = current[j - 1] + 1
            if left < cost:
                cost = left
            current[j] = cost
    ends = [(rows[n][j], n, j) for j in range(max(0, m - band), m + 1)]
    ends += [(rows[i][m], i, m) for i in range(max(0, n - band), n + 1)]
    edits, i, j = min(ends)
    return rows, ref_sets, edits, i, j


def gapped_stats(read: str, ref: str, band: int):
    """Best start-anchored, end-free banded alignment of ``read`` to ``ref``.

    Returns ``(edits, substitutions, gap_bases, gap_events, last_gap_end,
    ref_to_read)``: ``last_gap_end`` is the read index just after the last
    gap (so ``len(read) - last_gap_end`` bases verify it) and
    ``ref_to_read`` maps every aligned reference index to its read index;
    reference bases inside a deletion, and those past the aligned end, are
    absent from the map.
    """
    rows, ref_sets, edits, i, j = _gapped_dp(read, ref, band)
    ref_to_read: Dict[int, int] = {}
    substitutions = gap_bases = gap_events = 0
    last_gap_end = -1
    in_gap = False
    while i > 0 and j > 0:
        mismatch = read[i - 1] not in ref_sets[j - 1]
        if rows[i][j] == rows[i - 1][j - 1] + mismatch:
            ref_to_read[j - 1] = i - 1
            substitutions += mismatch
            in_gap = False
            i -= 1
            j -= 1
        else:
            if not in_gap:
                gap_events += 1
                in_gap = True
                last_gap_end = max(last_gap_end, i)
            gap_bases += 1
            if rows[i][j] == rows[i - 1][j] + 1:
                i -= 1
            else:
                j -= 1
    return edits, substitutions, gap_bases, gap_events, last_gap_end, ref_to_read


@functools.lru_cache(maxsize=4096)
def gapped_alignment(read: str, ref: str, band: int):
    """The alignment the annotations describe: ``(ops, end_read, end_ref, edits)``.

    Same start-anchored, end-free band and the same minimum edit count as
    gapped_stats, but among equal-edit explanations it takes the one with
    the fewest gap openings (a three-state table costing
    ``edits * GAP_OPEN_SCALE + openings``), so an inserted block whose bases
    happen to match the reference is one event rather than fragments, as an
    affine-gap aligner would describe it.  ``ops`` runs left to right; each
    is ``("M", i, j)`` (read base ``i`` aligned to reference base ``j``,
    matching or substituted), ``("I", i, j)`` (read base ``i`` inserted
    before reference index ``j``) or ``("D", i, j)`` (reference base ``j``
    deleted before read index ``i``).  ``end_read`` and ``end_ref`` are the
    aligned lengths: read bases from ``end_read`` on overhang the reference
    (a deletion pulled downstream sequence into the window) and reference
    bases from ``end_ref`` on were pushed out of it by an insertion.  Cached
    because a row is described more than once per run.
    """
    return _three_state_alignment(read, ref, band, open_penalty=1, sub_penalty=0)


@functools.lru_cache(maxsize=4096)
def fewest_substitution_alignment(read: str, ref: str, band: int):
    """The equal-edit explanation that needs the fewest substitutions.

    Same band, same minimum edit count as gapped_alignment; the tie-break
    prefers gaps over substitutions, so this is the explanation the
    composition guard judges: an indel haplotype is its canonical plus gaps
    and few substitutions, and that question must be asked of the
    explanation that spends the fewest substitutions, not of the one that
    opens the fewest gaps (which trades a gap opening for a substitution
    whenever it can).
    """
    return _three_state_alignment(read, ref, band, open_penalty=0, sub_penalty=1)


def _three_state_alignment(read: str, ref: str, band: int, open_penalty: int, sub_penalty: int):
    """Start-anchored, end-free banded alignment minimising ``edits * GAP_OPEN_SCALE``
    plus ``open_penalty`` per gap opening and ``sub_penalty`` per substitution."""
    n, m = len(read), len(ref)
    INF = 1 << 60
    scale = GAP_OPEN_SCALE
    ref_sets = [IUPAC_BASES.get(design, frozenset(design)) | {design} for design in ref]
    match_table = [[INF] * (m + 1) for _ in range(n + 1)]
    insert_table = [[INF] * (m + 1) for _ in range(n + 1)]
    delete_table = [[INF] * (m + 1) for _ in range(n + 1)]
    match_table[0][0] = 0
    sub_cost = scale + sub_penalty
    for j in range(1, min(m, band) + 1):
        delete_table[0][j] = j * scale + open_penalty
    for i in range(1, n + 1):
        lo, hi = max(0, i - band), min(m, i + band)
        if lo == 0:
            insert_table[i][0] = i * scale + open_penalty
        base = read[i - 1]
        match_row, insert_row, delete_row = match_table[i], insert_table[i], delete_table[i]
        match_up, insert_up, delete_up = match_table[i - 1], insert_table[i - 1], delete_table[i - 1]
        for j in range(max(1, lo), hi + 1):
            diagonal = min(match_up[j - 1], insert_up[j - 1], delete_up[j - 1])
            match_row[j] = diagonal + (sub_cost if base not in ref_sets[j - 1] else 0)
            insert_row[j] = min(match_up[j] + scale + open_penalty, insert_up[j] + scale, delete_up[j] + scale + open_penalty)
            delete_row[j] = min(match_row[j - 1] + scale + open_penalty, insert_row[j - 1] + scale + open_penalty, delete_row[j - 1] + scale)
    ends = [(min(match_table[n][j], insert_table[n][j], delete_table[n][j]), n, j) for j in range(max(0, m - band), m + 1)]
    ends += [(min(match_table[i][m], insert_table[i][m], delete_table[i][m]), i, m) for i in range(max(0, n - band), n + 1)]
    cost, i, j = min(ends)
    end_read, end_ref = i, j
    state = min((match_table[i][j], 0), (insert_table[i][j], 1), (delete_table[i][j], 2))[1]
    ops = []
    while i > 0 or j > 0:
        if state == 0:
            ops.append(("M", i - 1, j - 1))
            target = match_table[i][j] - (sub_cost if read[i - 1] not in ref_sets[j - 1] else 0)
            i -= 1
            j -= 1
            state = 0 if match_table[i][j] == target else (1 if insert_table[i][j] == target else 2)
        elif state == 1:
            ops.append(("I", i - 1, j))
            value = insert_table[i][j]
            i -= 1
            state = 1 if insert_table[i][j] == value - scale else (0 if match_table[i][j] == value - scale - open_penalty else 2)
        else:
            ops.append(("D", i, j - 1))
            value = delete_table[i][j]
            j -= 1
            state = 2 if delete_table[i][j] == value - scale else (0 if match_table[i][j] == value - scale - open_penalty else 1)
    ops.reverse()
    return tuple(ops), end_read, end_ref, cost // scale


def accepted_gapped_reading(
    sequence: str, template: str, gap_band: int, substitutions: Optional[int] = None
) -> Optional[GappedReading]:
    """The one gapped reading a row may be given, or ``None``.

    Classification and annotation share this: the fewest-openings alignment
    (``gapped_alignment``), taken only when it is strictly cheaper than the
    substitution-only reading (``substitutions``, the caller's bound; by
    default the mismatch count against ``template``), well composed
    (``gapped_explanation_ok``) and verified by at least GAP_VERIFIED_SUFFIX
    bases after its last gap.  A reading that fails any test is not a
    reading of the row at all: the row is then classified at the fixed
    offset and described as substitutions, by both consumers alike.
    """
    if gap_band <= 0 or len(template) != len(sequence):
        return None
    if substitutions is None:
        substitutions = tag_mismatches(sequence, template)
    if substitutions <= 1 or gapped_prescreen(sequence, template, gap_band) >= substitutions:
        return None
    ops, end_read, end_ref, edits = gapped_alignment(sequence, template, gap_band)
    if edits >= substitutions:
        return None
    # Composition is judged on the equal-edit explanation that needs the
    # fewest substitutions (the paralog question), the reading itself is the
    # fewest-openings one (the description and the target mapping).
    guard_ops = fewest_substitution_alignment(sequence, template, gap_band)[0]
    guard_subs = guard_bases = guard_events = 0
    previous = None
    for kind, i, j in guard_ops:
        if kind == "M":
            guard_subs += not _design_matches(sequence[i], template[j])
        else:
            guard_bases += 1
            guard_events += previous != kind
        previous = kind
    if not gapped_explanation_ok(guard_subs, guard_bases, guard_events):
        return None
    reading = _reading_from_ops(sequence, template, ops, end_read, end_ref, edits)
    if reading.last_gap_end > len(sequence) - GAP_VERIFIED_SUFFIX:
        return None
    return reading


def _reading_from_ops(sequence: str, template: str, ops, end_read: int, end_ref: int, edits: int) -> GappedReading:
    """Describe one alignment: substitution and gap counts, the verified-suffix boundary, the target map."""
    substitutions = gap_bases = gap_events = 0
    last_gap_end = -1
    ref_to_read: Dict[int, int] = {}
    previous = None
    for kind, i, j in ops:
        if kind == "M":
            ref_to_read[j] = i
            substitutions += not _design_matches(sequence[i], template[j])
        else:
            gap_bases += 1
            gap_events += previous != kind
            last_gap_end = max(last_gap_end, i + 1 if kind == "I" else i)
        previous = kind
    return GappedReading(edits, substitutions, gap_bases, gap_events, ops, end_read, end_ref, last_gap_end, ref_to_read)


def row_reading(locus: Locus, sequence: str, gap_band: int) -> Optional[GappedReading]:
    """The locus-level decision both consumers share: the accepted gapped reading
    of ``sequence`` against whichever canonical explains it most cheaply, or
    ``None``.  The bound a gapped reading must beat is the cheaper of the two
    substitution-only readings, and the composition guard is applied against
    the canonical the row is read from, so a designed substitution never counts
    against an AltMatch row's indel reading."""
    if gap_band <= 0:
        return None
    substitutions = min(tag_mismatches(sequence, locus.ref_tag), tag_mismatches(sequence, locus.alt_tag))
    best: Optional[GappedReading] = None
    for template in (locus.ref_tag, locus.alt_tag):
        accepted = accepted_gapped_reading(sequence, template, gap_band, substitutions)
        if accepted is not None and (best is None or accepted.edits < best.edits):
            best = accepted
    return best


def _left_align(ops, read: str, ref: str):
    """Shift every gap run left across matching bases (VCF-style left normalisation).

    Equivalent placements of an indel inside a repeat then describe it at
    one position, so the same haplotype yields the same tokens whichever
    equal-cost path the DP happened to take.  A substituted base blocks the
    shift: events never move across each other.
    """
    ops = list(ops)
    k = 0
    while k < len(ops):
        kind = ops[k][0]
        if kind == "M":
            k += 1
            continue
        run_end = k
        while run_end + 1 < len(ops) and ops[run_end + 1][0] == kind:
            run_end += 1
        while k > 0 and ops[k - 1][0] == "M":
            _, pi, pj = ops[k - 1]
            if read[pi] != ref[pj]:
                break
            if kind == "I":
                last_i = ops[run_end][1]
                if read[pi] != read[last_i]:
                    break
                ops[k - 1 : run_end + 1] = [("I", ops[t][1] - 1, pj) for t in range(k, run_end + 1)] + [("M", last_i, pj)]
            else:
                last_j = ops[run_end][2]
                if ref[pj] != ref[last_j]:
                    break
                ops[k - 1 : run_end + 1] = [("D", pi, ops[t][2] - 1) for t in range(k, run_end + 1)] + [("M", pi, last_j)]
            k -= 1
            run_end -= 1
        k = run_end + 1
    return ops


def gapped_explanation_ok(substitutions: int, gap_bases: int, gap_events: int) -> bool:
    """Composition guard for a gapped reading (see the module constants)."""
    return substitutions <= GAP_MAX_SUBSTITUTIONS and gap_events <= GAP_MAX_EVENTS


# ============================================================================
# 3. Panel definition (design files, reference tags, primer-free inference)
# ============================================================================


def _expanded_target(target: str, allele: str, match: re.Match) -> str:
    return (target[: match.start()] + allele + target[match.end() :]).upper()


def _find_primer(genomic: str, primer: str, start: int) -> int:
    """Locate a primer in genomic sequence, tolerating degenerate oligo bases.

    Vendors design primers with IUPAC codes over flanking SNPs; the concrete
    genome then never contains the literal primer string.  Literal find is the
    fast path; the compatibility scan runs only for degenerate primers.
    """
    position = genomic.find(primer, start)
    if position >= 0 or not (set(primer) - set("ACGT")):
        return position
    for position in range(start, len(genomic) - len(primer) + 1):
        if compatible_mismatches(genomic[position : position + len(primer)], primer, 0) == 0:
            return position
    return -1


def _candidate_tags(target: str, alleles: Sequence[str], left_primer: str, tag_length: int) -> List[str]:
    match = re.search(r"\[([^/\]]+)/([^\]]+)\]", target)
    if not match:
        raise MadcError("TargetSequence does not contain a bracketed variant such as [A/G]")
    candidates: List[str] = []
    for allele in alleles:
        genomic = _expanded_target(target, allele, match)
        # A tag is the primer-started window that carries the designed variant,
        # which occupies [variant_start, variant_end) of this expansion.  A
        # primer occurrence elsewhere in the target (a duplicated flank) starts
        # a window that never reaches the variant and is not a candidate.
        variant_start = match.start()
        variant_end = variant_start + max(1, len(allele.replace("-", "")))
        start = 0
        while True:
            position = _find_primer(genomic, left_primer, start)
            if position < 0:
                break
            if position + tag_length <= len(genomic) and position <= variant_start and variant_end <= position + tag_length:
                # Overlay the primer sequence: reads begin with the synthesized
                # oligo, so degenerate primer bases belong in the tag verbatim.
                candidates.append(left_primer + genomic[position + len(left_primer) : position + tag_length])
            start = position + 1

        primer_rc = reverse_complement(left_primer)
        start = 0
        while True:
            position = _find_primer(genomic, primer_rc, start)
            if position < 0:
                break
            lower = position - (tag_length - len(left_primer))
            upper = position + len(left_primer)
            if lower >= 0 and lower <= variant_start and variant_end <= upper:
                window = reverse_complement(genomic[lower:upper])
                candidates.append(left_primer + window[len(left_primer) :])
            start = position + 1
    return sorted(set(candidates))


def _indel_tag_pair(target: str, left_primer: str, tag_length: int, marker: str) -> Tuple[str, str, int]:
    """Reconstruct fixed-length Ref/Alt tags for an indel marker.

    The two allele expansions differ in length; taking tag-length windows at
    the primer position in each yields the vendor's shifted-suffix
    representation.  Requiring both windows to begin with the primer confines
    the indel to the tag body, and the pair must be unique across primer
    occurrences and strands.  Returns (ref_tag, alt_tag, first_difference).
    """
    match = BRACKET_RE.search(target)
    if not match:
        raise MadcError(f"No bracketed variant found in TargetSequence for {marker}")
    first = "" if match.group(1) == "-" else match.group(1)
    second = "" if match.group(2) == "-" else match.group(2)
    pairs = set()
    for strand in ("F", "R"):
        ref_genomic = target[: match.start()] + first + target[match.end():]
        alt_genomic = target[: match.start()] + second + target[match.end():]
        if strand == "R":
            ref_genomic = reverse_complement(ref_genomic)
            alt_genomic = reverse_complement(alt_genomic)
        start = 0
        while True:
            position = _find_primer(ref_genomic, left_primer, start)
            if position < 0:
                break
            if position + tag_length <= len(ref_genomic) and position + tag_length <= len(alt_genomic):
                ref_window = left_primer + ref_genomic[position + len(left_primer) : position + tag_length]
                alt_window = left_primer + alt_genomic[position + len(left_primer) : position + tag_length]
                alt_compatible = compatible_mismatches(
                    alt_genomic[position : position + len(left_primer)], left_primer, 0
                ) == 0
                if alt_compatible and ref_window != alt_window:
                    pairs.add((ref_window, alt_window))
            start = position + 1
    if not pairs:
        raise MadcError(f"Could not reconstruct indel tags for {marker}; check primer orientation")
    if len(pairs) > 1:
        raise MadcError(f"Ambiguous indel-tag reconstruction for {marker}")
    ref_tag, alt_tag = next(iter(pairs))
    first_diff = next(i for i, (a, b) in enumerate(zip(ref_tag, alt_tag)) if a != b)
    return ref_tag, alt_tag, first_diff


def build_panel(primer_csv: Path, probe_design: Path, tag_length: int) -> List[Locus]:
    primer_rows = read_dict_rows(primer_csv, delimiter=",")
    design_rows = read_dict_rows(probe_design)
    require_columns(
        primer_rows,
        [
            "MarkerName",
            "LeftPrimerSequence",
            "RightPrimerSequence",
            "VariantRef",
            "VariantAlt",
            "VariantPositionInTagSequence",
        ],
        primer_csv,
    )
    require_columns(design_rows, ["MarkerName", "TargetSequence"], probe_design)
    designs = {row["MarkerName"]: row for row in design_rows}
    if len(designs) != len(design_rows):
        raise MadcError(f"Duplicate MarkerName values in {probe_design}")

    loci: List[Locus] = []
    seen = set()
    for row in primer_rows:
        marker = row["MarkerName"]
        if marker in seen:
            raise MadcError(f"Duplicate MarkerName in {primer_csv}: {marker}")
        seen.add(marker)
        if marker not in designs:
            raise MadcError(f"Marker {marker} is absent from {probe_design}")
        left = row["LeftPrimerSequence"].upper()
        right = row["RightPrimerSequence"].upper()
        ref_base = row["VariantRef"].upper()
        alt_base = row["VariantAlt"].upper()
        if ref_base == alt_base:
            raise MadcError(f"VariantRef and VariantAlt are identical for {marker}")
        try:
            variant_pos = int(row["VariantPositionInTagSequence"]) - 1
        except ValueError as exc:
            raise MadcError(f"Invalid VariantPositionInTagSequence for {marker}") from exc
        if not 0 <= variant_pos < tag_length:
            raise MadcError(f"Variant position for {marker} is outside a {tag_length}-bp tag")
        if variant_pos < len(left):
            raise MadcError(
                f"Variant position for {marker} lies inside its left primer; Alt reads could never "
                "begin with the designed primer sequence"
            )

        target = designs[marker]["TargetSequence"].upper()
        if len(ref_base) == 1 and len(alt_base) == 1 and "-" not in (ref_base, alt_base):
            bracket = re.search(r"\[([^/\]]+)/([^\]]+)\]", target)
            if not bracket:
                raise MadcError(f"No bracketed variant found in TargetSequence for {marker}")
            candidates = _candidate_tags(target, bracket.groups(), left, tag_length)
            reference_candidates = [tag for tag in candidates if tag[variant_pos] == ref_base]
            if not reference_candidates:
                raise MadcError(
                    f"Could not reconstruct the reference tag for {marker}; check primer orientation and variant position"
                )
            if len(set(reference_candidates)) > 1:
                raise MadcError(f"Ambiguous reference-tag reconstruction for {marker}")
            ref_tag = reference_candidates[0]
            alt_tag = ref_tag[:variant_pos] + alt_base + ref_tag[variant_pos + 1 :]
            if not ref_tag.startswith(left):
                raise MadcError(f"Reconstructed tag for {marker} does not begin with its left primer")
            if alt_tag not in candidates:
                raise MadcError(
                    f"Constructed Alt tag for {marker} does not match any probe-design bracket expansion; "
                    "check VariantRef/VariantAlt and VariantPositionInTagSequence against the design"
                )
        else:
            # Indel marker (dash or multi-base alleles in the primer CSV; the
            # probe-design bracket is anchored).  Tags are same-length windows
            # of the two allele expansions; the pipeline's variant position is
            # the first differing base, which can sit at or beyond the CSV's
            # declared indel start when the flank repeats the indel content.
            ref_tag, alt_tag, first_diff = _indel_tag_pair(target, left, tag_length, marker)
            if first_diff < variant_pos:
                raise MadcError(
                    f"Indel tags for {marker} diverge at position {first_diff + 1}, before the declared "
                    f"VariantPositionInTagSequence {variant_pos + 1}"
                )
            variant_pos = first_diff
            ref_base = ref_tag[first_diff]
            alt_base = alt_tag[first_diff]
        loci.append(
            Locus(
                index=len(loci),
                marker=marker,
                left_primer=left,
                right_primer=right,
                ref_base=ref_base,
                alt_base=alt_base,
                variant_pos=variant_pos,
                ref_tag=ref_tag,
                alt_tag=alt_tag,
            )
        )
    return loci


def load_reference_tags(path: Path, tag_length: Optional[int] = None) -> List[Locus]:
    rows = read_dict_rows(path)
    require_columns(rows, ["MarkerName", "LeftPrimerSequence", "RefSequence", "AltSequence"], path)
    loci: List[Locus] = []
    seen = set()
    for row in rows:
        marker = row["MarkerName"]
        if marker in seen:
            raise MadcError(f"Duplicate MarkerName in {path}: {marker}")
        seen.add(marker)
        ref_tag = row["RefSequence"].upper()
        alt_tag = row["AltSequence"].upper()
        if len(ref_tag) != len(alt_tag):
            raise MadcError(f"RefSequence and AltSequence lengths differ for {marker}")
        if tag_length is not None and len(ref_tag) != tag_length:
            raise MadcError(f"{marker} has {len(ref_tag)} bases; expected {tag_length}")
        differences = [index for index, pair in enumerate(zip(ref_tag, alt_tag)) if pair[0] != pair[1]]
        if not differences:
            raise MadcError(f"RefSequence and AltSequence are identical for {marker}")
        if "VariantPositionInTagSequence" in row and row["VariantPositionInTagSequence"]:
            try:
                variant_pos = int(row["VariantPositionInTagSequence"]) - 1
            except ValueError as exc:
                raise MadcError(f"Invalid VariantPositionInTagSequence for {marker} in {path}") from exc
            if not 0 <= variant_pos < len(ref_tag):
                raise MadcError(
                    f"VariantPositionInTagSequence for {marker} is outside its {len(ref_tag)}-bp tag"
                )
            if variant_pos not in differences:
                raise MadcError(
                    f"VariantPositionInTagSequence for {marker} points at position {variant_pos + 1}, "
                    "where RefSequence and AltSequence do not differ"
                )
        elif len(differences) == 1:
            variant_pos = differences[0]
        else:
            raise MadcError(
                f"{marker} requires VariantPositionInTagSequence because RefSequence and AltSequence differ at "
                f"{len(differences)} positions"
            )
        ref_base = (row.get("VariantRef") or ref_tag[variant_pos]).upper()
        alt_base = (row.get("VariantAlt") or alt_tag[variant_pos]).upper()
        if len(ref_base) != 1 or len(alt_base) != 1:
            raise MadcError(f"Only single-base variant alleles are supported ({marker}: {ref_base}/{alt_base})")
        if ref_base == alt_base:
            raise MadcError(f"VariantRef and VariantAlt are identical for {marker}")
        if not bases_compatible(ref_base, ref_tag[variant_pos]):
            raise MadcError(
                f"VariantRef {ref_base!r} for {marker} is incompatible with RefSequence base "
                f"{ref_tag[variant_pos]!r} at the variant position"
            )
        if not bases_compatible(alt_base, alt_tag[variant_pos]):
            raise MadcError(
                f"VariantAlt {alt_base!r} for {marker} is incompatible with AltSequence base "
                f"{alt_tag[variant_pos]!r} at the variant position"
            )
        left = row["LeftPrimerSequence"].upper()
        if not ref_tag.startswith(left) or not alt_tag.startswith(left):
            raise MadcError(f"Reference sequences for {marker} do not begin with LeftPrimerSequence")
        # A variant inside the left primer is impossible here: the declared
        # position is a Ref/Alt difference, and both tags begin with the same
        # primer, so the difference must lie beyond it.
        loci.append(
            Locus(
                index=len(loci),
                marker=marker,
                left_primer=left,
                right_primer=(row.get("RightPrimerSequence") or "").upper(),
                ref_base=ref_base,
                alt_base=alt_base,
                variant_pos=variant_pos,
                ref_tag=ref_tag,
                alt_tag=alt_tag,
                locus_group=(row.get("LocusGroup") or "").strip(),
            )
        )
    return loci


def write_reference_tags(path: Path, loci: Sequence[Locus]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    # The LocusGroup column appears only when some locus actually uses
    # grouping, so every pre-tiling panel stays byte-identical.
    grouped = any(locus.locus_group for locus in loci)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        header = [
            "MarkerName",
            "LeftPrimerSequence",
            "RightPrimerSequence",
            "RefSequence",
            "AltSequence",
            "VariantRef",
            "VariantAlt",
            "VariantPositionInTagSequence",
        ]
        if grouped:
            header.append("LocusGroup")
        writer.writerow(header)
        for locus in loci:
            row = [
                locus.marker,
                locus.left_primer,
                locus.right_primer,
                locus.ref_tag,
                locus.alt_tag,
                locus.ref_base,
                locus.alt_base,
                locus.variant_pos + 1,
            ]
            if grouped:
                row.append(locus.locus_group)
            writer.writerow(row)


BRACKET_RE = re.compile(r"\[([^/\]]+)/([^\]]+)\]")


def read_probe_targets(path: Path) -> Dict[str, str]:
    """Load MarkerName -> bracket-notation TargetSequence for infer-design."""
    return read_probe_design(path)[0]


def read_probe_design(path: Path) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Load targets plus each marker's locus group (marker itself when absent).

    Markers sharing a LocusGroup are tiles of one physical locus: their
    TargetSequences must be the same region (identical Ref expansions),
    differing only in which variant is bracketed.
    """
    rows = read_dict_rows(path)
    require_columns(rows, ["MarkerName", "TargetSequence"], path)
    targets: Dict[str, str] = {}
    groups: Dict[str, str] = {}
    for row in rows:
        marker = row["MarkerName"]
        if marker in targets:
            raise MadcError(f"Duplicate MarkerName in {path}: {marker}")
        target = row["TargetSequence"].upper()
        for bracket in BRACKET_RE.finditer(target):
            if "/" in bracket.group(2):
                raise MadcError(
                    f"{marker} in {path}: bracket [{bracket.group(1)}/{bracket.group(2)}] names more "
                    "than two alleles; multi-allelic targets are not supported — split the marker "
                    "into biallelic rows or drop it"
                )
        targets[marker] = target
        groups[marker] = (row.get("LocusGroup") or "").strip() or marker
    return targets, groups


def enumerate_variant_windows(
    target: str, tag_length: int
) -> Dict[Tuple[str, int], Tuple[str, str, int]]:
    """Every tag-length window pair of the expanded target that spans the variant.

    The tag must be such a window on one strand — reads are primer-started
    and the assay interrogates the bracketed variant — so this is a closed
    candidate set for tag inference.  Returns ``{(strand, start): (ref_window,
    alt_window, variant_offset)}`` with Ref built from the FIRST bracket
    allele, matching how DArTag designs are generated from the reference
    genome.  Indel brackets are supported: the two allele expansions differ in
    length, and taking fixed windows at the same start yields the vendor's
    shifted-suffix representation.  The divergence point is the shared prefix
    of the full expansions, which also handles anchored-base and homopolymer
    brackets correctly; ``variant_offset`` is the first differing position
    within the window pair.  Empty when the alleles are identical or the
    flanks are too short for any full window.
    """
    match = BRACKET_RE.search(target)
    if not match:
        return {}
    first, second = match.group(1), match.group(2)
    first = "" if first == "-" else first
    second = "" if second == "-" else second
    if first == second:
        return {}
    windows: Dict[Tuple[str, int], Tuple[str, str, int]] = {}
    for strand in ("F", "R"):
        ref_genomic = target[: match.start()] + first + target[match.end():]
        alt_genomic = target[: match.start()] + second + target[match.end():]
        if strand == "R":
            ref_genomic = reverse_complement(ref_genomic)
            alt_genomic = reverse_complement(alt_genomic)
        divergence = 0
        for a, b in zip(ref_genomic, alt_genomic):
            if a != b:
                break
            divergence += 1
        lower = max(0, divergence - tag_length + 1)
        upper = min(divergence, len(ref_genomic) - tag_length, len(alt_genomic) - tag_length)
        for start in range(lower, upper + 1):
            ref_window = ref_genomic[start : start + tag_length]
            alt_window = alt_genomic[start : start + tag_length]
            if ref_window != alt_window:
                windows[(strand, start)] = (ref_window, alt_window, divergence - start)
    return windows


# Worker-side state for the infer-design FASTQ scans.  Each worker rebuilds
# its lookup tables from the probe-design path in the initializer instead of
# receiving them pickled, because the candidate dictionaries are large.
_INFER_STATE: Dict[str, object] = {}


def _infer_exact_init(probe_design: str, tag_length: int, leading_anchor: str = "") -> None:
    targets = read_probe_targets(Path(probe_design))
    lookup: Dict[str, Tuple[Tuple[str, Tuple[str, int]], ...]] = {}
    staging: Dict[str, set] = {}
    for marker, target in targets.items():
        for key, (ref_window, alt_window, _offset) in enumerate_variant_windows(
            target, tag_length
        ).items():
            staging.setdefault(ref_window, set()).add((marker, key))
            staging.setdefault(alt_window, set()).add((marker, key))
    for sequence, hits in staging.items():
        lookup[sequence] = tuple(sorted(hits))
    _INFER_STATE.clear()
    _INFER_STATE.update({"tag_length": tag_length, "exact_lookup": lookup, "leading_anchor": leading_anchor})


def _infer_exact_scan(fastq: str):
    """Count exact read-prefix hits against every candidate window."""
    tag_length = _INFER_STATE["tag_length"]
    lookup = _INFER_STATE["exact_lookup"]
    votes: Counter = Counter()
    for sequence, _quality in fastq_records(Path(fastq), _INFER_STATE.get("leading_anchor", "")):
        if len(sequence) < tag_length:
            continue
        hits = lookup.get(sequence[:tag_length])
        if hits:
            for marker_key in hits:
                votes[marker_key] += 1
    return fastq, votes


def _infer_fuzzy_init(
    probe_design: str,
    tag_length: int,
    max_mismatches: int,
    unresolved: Sequence[str],
    sink_tags: Sequence[str],
    leading_anchor: str = "",
) -> None:
    targets = read_probe_targets(Path(probe_design))
    windows = {
        marker: enumerate_variant_windows(targets[marker], tag_length)
        for marker in unresolved
    }
    _INFER_STATE["leading_anchor"] = leading_anchor
    seed_length = max(12, min(18, tag_length // 3))
    offsets = sorted({*range(0, tag_length - seed_length + 1, seed_length), tag_length - seed_length})
    seeds: Dict[str, set] = {}
    for marker, marker_windows in windows.items():
        for key, (ref_window, alt_window, _offset) in marker_windows.items():
            for window in (ref_window, alt_window):
                for offset in offsets:
                    seeds.setdefault(window[offset : offset + seed_length], set()).add((marker, key))
    # Resolved loci's tags act as competition sinks: a read closer to an
    # already-anchored locus must not vote for an unresolved candidate, or
    # paralogous reads would anchor dead loci to the wrong window.
    for tag in sink_tags:
        for offset in offsets:
            seeds.setdefault(tag[offset : offset + seed_length], set()).add(("", tag))
    _INFER_STATE.clear()
    _INFER_STATE.update(
        {
            "tag_length": tag_length,
            "max_mismatches": max_mismatches,
            "windows": windows,
            "seeds": {seed: tuple(sorted(hits)) for seed, hits in seeds.items()},
            "offsets": offsets,
            "seed_length": seed_length,
            "sinks": tuple(sink_tags),
        }
    )


def _infer_fuzzy_scan(fastq: str):
    """Vote reads onto candidate windows within the identity envelope."""
    tag_length = _INFER_STATE["tag_length"]
    max_mismatches = _INFER_STATE["max_mismatches"]
    windows = _INFER_STATE["windows"]
    seeds = _INFER_STATE["seeds"]
    offsets = _INFER_STATE["offsets"]
    seed_length = _INFER_STATE["seed_length"]
    votes: Counter = Counter()
    for sequence, _quality in fastq_records(Path(fastq), _INFER_STATE.get("leading_anchor", "")):
        if len(sequence) < tag_length:
            continue
        tag = sequence[:tag_length]
        candidates = set()
        for offset in offsets:
            hits = seeds.get(tag[offset : offset + seed_length])
            if hits:
                candidates.update(hits)
        if not candidates:
            continue
        scored = []
        for marker, key in candidates:
            if marker == "":
                distance = compatible_mismatches(tag, key, max_mismatches)
            else:
                ref_window, alt_window, _offset = windows[marker][key]
                distance = min(
                    compatible_mismatches(tag, ref_window, max_mismatches),
                    compatible_mismatches(tag, alt_window, max_mismatches),
                )
            if distance <= max_mismatches:
                scored.append((distance, marker, key))
        if not scored:
            continue
        scored.sort(key=lambda item: item[0])
        best_distance = scored[0][0]
        best = [(marker, key) for distance, marker, key in scored if distance == best_distance]
        if len(best) != 1 or best[0][0] == "":
            continue  # ambiguous between loci, or a resolved-locus sink won
        votes[best[0]] += 1
    return fastq, votes


def choose_anchor_lengths(
    resolved: Mapping[str, Tuple[str, int]], minimum: int, maximum: int
) -> Dict[str, int]:
    """Shortest panel-unique tag prefix per locus, never reaching the variant.

    The anchor substitutes for the unknowable left primer: the trie requires
    it to match exactly, so shorter anchors reject fewer error-carrying reads,
    but anchors must distinguish loci.  Lengths grow one base at a time only
    for colliding groups; loci identical up to their caps keep a shared
    anchor, which the assignment stage already resolves by tag identity.
    """
    lengths = {
        marker: min(max(1, minimum), min(maximum, variant_pos))
        for marker, (_tag, variant_pos) in resolved.items()
    }
    caps = {
        marker: min(maximum, variant_pos) for marker, (_tag, variant_pos) in resolved.items()
    }
    while True:
        groups: Dict[str, List[str]] = defaultdict(list)
        for marker, (tag, _pos) in resolved.items():
            groups[tag[: lengths[marker]]].append(marker)
        grew = False
        for members in groups.values():
            if len(members) == 1:
                continue
            for marker in members:
                if lengths[marker] < caps[marker]:
                    lengths[marker] += 1
                    grew = True
        if not grew:
            return lengths


def panel_fingerprint(loci: Sequence[Locus]) -> str:
    payload = [dataclasses.asdict(locus) for locus in loci]
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


# ============================================================================
# 4. Platform intake (pair merging, window location, tiling and routing, window-scan)
# ============================================================================


# Worker-side state for normalize-reads scans, rebuilt per worker from the
# probe-design path (the k-mer tables are too large to pickle per task).
_NORM_STATE: Dict[str, object] = {}


def _normalize_init(probe_design: str, kmer_size: int) -> None:
    """Build the k-mer location index over both allele expansions per locus.

    Every k-mer maps to (locus group, position), which is what lets an
    unoriented, variable-start read be placed onto its locus.  Tiles of one
    locus (markers sharing a LocusGroup) share one k-mer territory — their
    Ref expansions are required to be identical, so within-group k-mers
    reinforce rather than collide.  A k-mer occurring in two DIFFERENT
    groups, or at two positions within one group's frame (an internal
    repeat), is deleted outright rather than assigned arbitrarily:
    overlapping or paralogous loci must make reads look ambiguous, not
    silently steer them to one of the contenders, and a repeated k-mer
    must not place a read at the wrong copy.
    """
    targets, group_of = read_probe_design(Path(probe_design))
    expansions: Dict[str, Tuple[str, str, int]] = {}
    members_of: Dict[str, List[str]] = defaultdict(list)
    frame_of: Dict[str, str] = {}
    index: Dict[str, Tuple[str, int]] = {}
    collide = set()
    for marker, target in targets.items():
        match = BRACKET_RE.search(target)
        if not match:
            continue
        first, second = match.group(1), match.group(2)
        if len(first) != 1 or len(second) != 1 or first == second:
            continue  # normalize-reads currently supports SNP markers
        ref_g = target[: match.start()] + first + target[match.end():]
        alt_g = target[: match.start()] + second + target[match.end():]
        group = group_of.get(marker, marker)
        frame = frame_of.get(group)
        if frame is None:
            frame_of[group] = ref_g
        elif frame != ref_g:
            raise MadcError(
                f"LocusGroup {group}: member {marker}'s Ref expansion differs from the "
                "group frame; tiles of one locus must share an identical target region "
                "with only the bracketed variant differing"
            )
        expansions[marker] = (ref_g, alt_g, match.start())
        members_of[group].append(marker)
        for g in (ref_g, alt_g):
            for p in range(len(g) - kmer_size + 1):
                kmer = g[p : p + kmer_size]
                prev = index.get(kmer)
                # A k-mer seen in another group, or at a second position
                # within this group's frame, cannot place a read (the Ref
                # and Alt expansions legitimately share k-mers at the same
                # position).
                if prev and (prev[0] != group or prev[1] != p):
                    collide.add(kmer)
                index[kmer] = (group, p)
    for kmer in collide:
        del index[kmer]
    _NORM_STATE.clear()
    _NORM_STATE.update({
        "index": index, "expansions": expansions, "k": kmer_size,
        "members_of": dict(members_of),
        "group_of": {m: group_of.get(m, m) for m in expansions},
    })


def _normalize_locate(sequence: str):
    """Majority-vote marker location: (marker, hits) / (None, None) / ("AMBIG", None).

    A read whose k-mers hit several markers is assigned only when one marker
    holds at least twice the runner-up's hits; error-mutated k-mers
    occasionally collide with another marker, and a single stray hit must not
    discard an otherwise clearly located read.
    """
    index = _NORM_STATE["index"]
    k = _NORM_STATE["k"]
    hits: Dict[str, List[Tuple[int, int]]] = {}
    for rp in range(0, len(sequence) - k + 1, 7):
        got = index.get(sequence[rp : rp + k])
        if got:
            hits.setdefault(got[0], []).append((rp, got[1]))
    if not hits:
        return None, None
    if len(hits) > 1:
        ranked = sorted(hits.items(), key=lambda kv: -len(kv[1]))
        if len(ranked[0][1]) < 2 * len(ranked[1][1]):
            return "AMBIG", None
        return ranked[0][0], ranked[0][1]
    (marker, positions), = hits.items()
    return marker, positions


def _normalize_coverage_scan(fastq: str):
    """Pass 1: per-marker (read start offset, read length) observations."""
    cover: Counter = Counter()
    for sequence, _quality in fastq_records(Path(fastq)):
        for probe in (sequence, reverse_complement(sequence)):
            marker, positions = _normalize_locate(probe)
            if marker not in (None, "AMBIG"):
                rp, tp = positions[0]
                cover[(marker, tp - rp, len(sequence))] += 1
                break
    return fastq, cover


def choose_tag_windows(
    expansions: Mapping[str, Tuple[str, str, int]],
    cover: Mapping[Tuple[str, int, int], int],
    tag_length: int,
    anchor_length: int,
    placement: str,
    str_penalty: float = 0.002,
    groups: Optional[Mapping[str, str]] = None,
) -> Dict[str, Dict[str, object]]:
    """Place each marker's tag window, honoring the anchor-safety invariant.

    The window must contain the variant but keep it out of the first
    ``anchor_length`` bases, so Ref and Alt tags share their assignment
    anchor.  ``auto`` placement maximizes the fraction of observed reads
    that fully cover the window minus a soft tandem-repeat penalty of
    ``str_penalty`` per repeat bp beyond 4 (two-chip honeybee replays showed
    coverage is effectively never exactly tied, so a pure tiebreak cannot
    steer placement away from repeats; the soft penalty trades ~0.02%
    coverage and rescued a locus from 9% to 100% vendor concordance).
    Remaining ties break toward the shortest tandem run, then centering;
    ``center`` needs no read data.  ``str_penalty=0`` restores pure
    coverage ranking.

    When ``groups`` maps markers to locus groups, member windows of one
    group must excise DISTINCT sequences: tiles share a coverage profile,
    so the optimizer can collapse them onto one start, and two markers
    with identical Ref tags make every ref-carrying read an exact tie —
    ties are dropped as ambiguous, silently zeroing the locus's ref depth
    and flipping its calls (measured on GT-seq: 100% -> 3% concordance).
    A colliding member takes its best-ranked start with unclaimed window
    content instead; no distinct start at all is a design error.
    """
    per_marker: Dict[str, List[Tuple[int, int, int]]] = defaultdict(list)
    for (marker, offset, length), count in cover.items():
        per_marker[marker].append((offset, length, count))
    windows: Dict[str, Dict[str, object]] = {}
    claimed: Dict[Tuple[str, str], str] = {}  # (group, window content) -> marker
    for marker, (ref_g, _alt_g, snp) in expansions.items():
        limit = len(ref_g) - tag_length
        lower = max(0, snp - tag_length + 1)
        upper = min(snp - anchor_length, limit)
        if upper < lower:
            continue
        group = groups.get(marker, marker) if groups else marker
        grouped = groups is not None and sum(1 for g in groups.values() if g == group) > 1

        def tandem_run(start: int) -> int:
            return max_tandem_run(ref_g[start : start + tag_length])

        def content(start: int) -> str:
            return ref_g[start : start + tag_length]

        def taken(start: int) -> bool:
            return grouped and (group, content(start)) in claimed

        def claim(start: int) -> None:
            if grouped:
                claimed[(group, content(start))] = marker

        if placement == "center":
            start = min(max(lower, snp - tag_length // 2), upper)
            if taken(start):
                start = next(
                    (c for delta in range(1, upper - lower + 1)
                     for c in (start - delta, start + delta)
                     if lower <= c <= upper and not taken(c)),
                    None,
                )
                if start is None:
                    raise MadcError(
                        f"LocusGroup {group}: no window start gives {marker} a distinct "
                        f"excised sequence from its co-members; tiles this close cannot "
                        f"be told apart at tag length {tag_length}"
                    )
            claim(start)
            windows[marker] = {"start": start, "coverage": None, "max_run": tandem_run(start)}
            continue
        observations = per_marker.get(marker, [])
        total = sum(count for _, _, count in observations) or 1
        scored = []
        for start in range(lower, upper + 1):
            covered = sum(
                count for offset, length, count in observations
                if offset <= start and offset + length >= start + tag_length
            )
            fraction = covered / total
            run = tandem_run(start)
            scored.append(
                (
                    fraction - str_penalty * max(0, run - 4),
                    -run,
                    -abs(start + tag_length // 2 - snp),
                    start,
                    fraction,
                )
            )
        scored.sort(reverse=True)
        pick = next((entry for entry in scored if not taken(entry[3])), None)
        if pick is None:
            raise MadcError(
                f"LocusGroup {group}: no window start gives {marker} a distinct "
                f"excised sequence from its co-members; tiles this close cannot "
                f"be told apart at tag length {tag_length}"
            )
        _, neg_run, _, start, fraction = pick
        claim(start)
        nudged = pick is not scored[0]
        windows[marker] = {"start": start, "coverage": fraction, "max_run": -neg_run}
        if nudged:
            windows[marker]["dedup_nudged_from"] = scored[0][3]
    return windows


def _normalize_extract(arguments):
    """Pass 2: excise each read's tag window and write the normalized FASTQ.

    A read is located (retrying reverse-complemented, with qualities
    reversed to stay aligned), then the window is cut using the k-mer hit
    nearest the window start — distant hits can straddle indels that shift
    the local offset, so the nearest anchor gives the most faithful
    coordinate mapping.  Reads whose span does not fully cover the window
    are counted out rather than padded: a partial window would fabricate
    sequence at the exact bases downstream matching relies on.
    """
    fastq, out_dir, tag_length, window_starts, members_of, routing = arguments
    sample = strip_fastq_suffix(Path(fastq).name)
    stats: Counter = Counter()
    out_path = Path(out_dir) / f"{sample}.fastq"
    with out_path.open("w") as out:
        number = 0
        for sequence, quality in fastq_records(Path(fastq)):
            number += 1
            group, positions = _normalize_locate(sequence)
            if group is None:
                flipped = reverse_complement(sequence)
                group, positions = _normalize_locate(flipped)
                if group not in (None, "AMBIG"):
                    sequence, quality = flipped, quality[::-1]
            if group is None:
                stats["reads_no_panel_kmer"] += 1
                continue
            if group == "AMBIG":
                stats["reads_ambiguous_marker"] += 1
                continue
            candidates = [m for m in members_of.get(group, ()) if m in window_starts]
            if not candidates:
                stats["reads_marker_without_window"] += 1
                continue
            # Routing among a group's member windows.  'nearest' feeds each
            # read to EXACTLY ONE tile — of the member windows the read fully
            # covers, the one whose center is nearest the read's center
            # (deterministic; ties break toward the lower window start) — so
            # a window can never count one fragment twice.  'covered' emits
            # the excision at EVERY member window the read fully covers: a
            # merged full-span fragment genuinely observes each tile's SNP
            # once, and nearest-routing such fragments is winner-take-all
            # (they all share the amplicon center), starving the other tiles.
            # Each window still sees each fragment at most once either way.
            covered = []
            for member in candidates:
                window_start = window_starts[member]
                window_end = window_start + tag_length
                rp, tp = min(positions, key=lambda h: abs(h[1] - window_start))
                offset = tp - rp
                r0, r1 = window_start - offset, window_end - offset
                if r0 < 0 or r1 > len(sequence):
                    continue
                read_center = offset + len(sequence) / 2
                distance = abs((window_start + window_end) / 2 - read_center)
                covered.append(((distance, window_start), r0, r1))
            if not covered:
                stats["reads_window_not_covered"] += 1
                continue
            if routing == "covered":
                emit = sorted(covered)
            else:
                emit = [min(covered)]
            for _, r0, r1 in emit:
                out.write(f"@{sample}_{number}\n{sequence[r0:r1]}\n+\n{quality[r0:r1]}\n")
            stats["reads_extracted"] += 1
            if len(emit) > 1:
                stats["tile_extra_emissions"] += len(emit) - 1
    return sample, stats


@functools.lru_cache(maxsize=64)
def _overlap_shift_order(len1: int, len2: int, min_overlap: int) -> Tuple[int, ...]:
    """Candidate merge shifts ordered longest-overlap-first.

    Overlap length is unimodal in the shift, so first-accept iteration must
    run from the largest overlap down: inside a tandem repeat every in-phase
    shift matches perfectly, and whichever is tried first wins.  Equal
    overlaps prefer the larger shift (the extended geometry).  Read lengths
    are near-constant across a run, so the order is computed once per
    (len1, len2) and cached.
    """
    shifts = range(len1 - min_overlap, -(len2 - min_overlap) - 1, -1)
    return tuple(
        sorted(shifts, key=lambda s: (min(len1, s + len2) - max(0, s), s), reverse=True)
    )


def merge_read_pair(seq1, qual1, seq2, qual2, min_overlap, max_mismatch_rate):
    """Merge one read pair by 3' overlap; None when no acceptable overlap.

    R2 is reverse-complemented so both reads run fragment-forward, then the
    best-scoring acceptable overlap wins (score = overlap - 5*mismatches):
    the true overlap for amplicon/capture pairs is the full shared span at
    near-zero mismatches, and neither a shorter spurious match (which would
    insert phantom sequence) nor a longer one riding the mismatch budget
    (which would delete real sequence) may preempt it.  The merged read is
    always R1's head + the overlap + R2's tail: in the extended geometry
    (insert longer than either read) that is the whole fragment; in the
    engulfed / read-through geometry (insert shorter than a read) R1's
    bases past the fragment end are adapter and fall outside the overlap,
    so they are discarded by construction, while any R2 bases past R1's
    end (reads of unequal length) are fragment and kept.  Disagreeing
    overlap bases take the higher-quality call.
    """
    seq2 = reverse_complement(seq2)
    qual2 = qual2[::-1]
    len1, len2 = len(seq1), len(seq2)
    # shift = position of seq2's start relative to seq1's start.  Candidates
    # are scanned longest-overlap-first and scored as overlap - 5*mismatches
    # (each mismatch forfeits its own position plus a four-base penalty), and
    # the best-scoring acceptable overlap wins.  First-accept in either
    # direction mis-merges tandem repeats: shortest-first takes an in-phase
    # sub-overlap and inserts phantom repeat units, while longest-first lets
    # a longer in-phase overlap ride the mismatch budget past the true
    # zero-mismatch span.  Scanning in descending overlap order lets the
    # scorer stop as soon as no remaining overlap could beat the best score.
    best_score = -1
    best_shift = None
    for shift in _overlap_shift_order(len1, len2, min_overlap):
        lo1, hi1 = max(0, shift), min(len1, shift + len2)
        overlap = hi1 - lo1
        if overlap <= best_score:
            break
        mismatches = 0
        budget = int(overlap * max_mismatch_rate)
        for i in range(lo1, hi1):
            if seq1[i] != seq2[i - shift]:
                mismatches += 1
                if mismatches > budget:
                    break
        else:
            score = overlap - 5 * mismatches
            if score > best_score:
                best_score = score
                best_shift = shift
    if best_shift is None:
        return None
    shift = best_shift
    lo1, hi1 = max(0, shift), min(len1, shift + len2)
    merged = []
    merged_q = []
    for i in range(lo1, hi1):
        j = i - shift
        if seq1[i] == seq2[j] or qual1[i] >= qual2[j]:
            merged.append(seq1[i])
            merged_q.append(max(qual1[i], qual2[j]))
        else:
            merged.append(seq2[j])
            merged_q.append(qual2[j])
    core_seq, core_qual = "".join(merged), "".join(merged_q)
    # R1 head + overlap + R2 tail covers both geometries: rc(R2) ends at
    # the fragment end by construction, so any R2 bases past the overlap
    # are fragment (kept, also when the reads differ in length), while R1
    # bases past the overlap are adapter read-through (dropped by ``hi1``).
    # In the engulfed equal-length case both tails are empty and the
    # overlap IS the fragment.
    return (seq1[:lo1] + core_seq + seq2[hi1 - shift:],
            qual1[:lo1] + core_qual + qual2[hi1 - shift:])


def _merge_pairs_sample(task):
    """Merge one sample's pair of FASTQ files (worker for run_merge_pairs)."""
    r1_path, r2_path, out_path, sample, min_overlap, max_mismatch_rate, unmerged = task
    stats: Counter = Counter()
    with gzip.open(out_path, "wt") as out:
        number = 0
        for mate1, mate2 in itertools.zip_longest(
            fastq_records_named(Path(r1_path)), fastq_records_named(Path(r2_path))
        ):
            number += 1
            if mate1 is None or mate2 is None:
                raise MadcError(f"{sample}: R1 and R2 hold different numbers of reads (record {number})")
            (name1, seq1, qual1), (name2, seq2, qual2) = mate1, mate2
            if mate_id(name1) != mate_id(name2):
                raise MadcError(
                    f"{sample}: mates out of step at record {number}: {name1!r} vs {name2!r}"
                )
            merged = merge_read_pair(seq1, qual1, seq2, qual2, min_overlap, max_mismatch_rate)
            if merged is not None:
                out.write(f"@{sample}_{number}_m\n{merged[0]}\n+\n{merged[1]}\n")
                stats["pairs_merged"] += 1
            elif unmerged == "r1":
                out.write(f"@{sample}_{number}_1\n{seq1}\n+\n{qual1}\n")
                stats["pairs_unmerged_r1_kept"] += 1
            elif unmerged == "pool":
                out.write(f"@{sample}_{number}_1\n{seq1}\n+\n{qual1}\n")
                out.write(f"@{sample}_{number}_2\n{seq2}\n+\n{qual2}\n")
                stats["pairs_unmerged_pooled"] += 1
            else:
                stats["pairs_unmerged_dropped"] += 1
    return sample, stats


def scan_churn(read_accounting: Mapping[str, int]) -> int:
    """Assigned reads that did not survive as untouched evidence: moved by either
    denoising pass, dropped as failed evidence, or dropped as unresolvable N reads
    (those are assigned, then removed by ``resolve_n_tags``)."""
    return (
        read_accounting["reads_moved_by_sample_denoise"]
        + read_accounting["reads_moved_by_global_denoise"]
        + read_accounting["reads_dropped_failed_evidence"]
        + read_accounting.get("reads_n_dropped", 0)
    )


def select_scan_subset(fastqs: Sequence[Path], count: int, seed: Optional[int]) -> List[Path]:
    """The FASTQs a scan runs on: a seeded random draw, or the first ``count`` in sorted order."""
    if seed is None or len(fastqs) <= count:
        return list(fastqs[:count])
    return sorted(random.Random(seed).sample(list(fastqs), count))


# ============================================================================
# 5. Per-sample read assignment (FASTQ, primer index, assignment, cache)
# ============================================================================


class PrimerIndex:
    """Compact trie that reports all best primer-prefix matches."""

    def __init__(self, loci: Sequence[Locus]):
        self.children: List[Dict[str, int]] = [{}]
        self.terminals: List[List[int]] = [[]]
        self.max_length = 0
        self.has_degenerate = False
        for locus in loci:
            node = 0
            self.max_length = max(self.max_length, len(locus.left_primer))
            if any(base not in "ACGT" for base in locus.left_primer):
                self.has_degenerate = True
            for base in locus.left_primer:
                child = self.children[node].get(base)
                if child is None:
                    child = len(self.children)
                    self.children[node][base] = child
                    self.children.append({})
                    self.terminals.append([])
                node = child
            self.terminals[node].append(locus.index)

    def match(self, sequence: str, max_mismatches: int = 0) -> List[int]:
        """Return locus candidates whose full left primer matches the read start.

        Terminals from every matched depth are returned (a primer that is a
        strict prefix of another cannot shadow its locus), and with mismatches
        allowed every primer within the budget is returned regardless of depth
        or mismatch count — an exact short primer must not suppress a longer
        one-mismatch primer.  Matching is IUPAC-aware on the primer side: an
        ambiguity code in the oligo names the bases the pool contains, so any
        member base matches without spending the mismatch budget (a read N
        satisfies nothing and still counts as a mismatch).  The caller
        resolves multi-locus candidates by full-tag identity, which is the
        stronger signal.
        """
        if max_mismatches == 0 and not self.has_degenerate:
            node = 0
            matches: List[int] = []
            for base in sequence[: self.max_length]:
                child = self.children[node].get(base)
                if child is None:
                    break
                node = child
                if self.terminals[node]:
                    matches.extend(self.terminals[node])
            return matches

        states = [(0, 0)]
        candidates: set = set()
        for base in sequence[: self.max_length]:
            next_states: List[Tuple[int, int]] = []
            for node, mismatches in states:
                for expected, child in self.children[node].items():
                    accepted = expected == base or base in IUPAC_BASES.get(expected, frozenset())
                    new_mismatches = mismatches + (0 if accepted else 1)
                    if new_mismatches <= max_mismatches:
                        next_states.append((child, new_mismatches))
                        candidates.update(self.terminals[child])
            if not next_states:
                break
            states = next_states
        return sorted(candidates)


def fastq_records(path: Path, leading_anchor: str = "") -> Iterator[Tuple[str, str]]:
    """``(sequence, quality)`` records.  With ``leading_anchor`` (raw deliveries
    that keep the library's common sequence in front of the amplicon) each
    read starts where the anchor ends, and reads without it are skipped;
    callers that must count those reads cut for themselves (``process_sample``)."""
    search_end = len(leading_anchor) + LEADING_ANCHOR_SEARCH
    with open_text(path, "rt") as handle:
        record_number = 0
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not plus or not quality:
                raise MadcError(f"Truncated FASTQ record {record_number} in {path}")
            if not name.startswith("@") or not plus.startswith("+"):
                raise MadcError(f"Invalid FASTQ record {record_number} in {path}")
            sequence = sequence.strip().upper()
            quality = quality.strip()
            if len(sequence) != len(quality):
                raise MadcError(f"Sequence/quality length mismatch at FASTQ record {record_number} in {path}")
            if leading_anchor:
                cut = sequence.find(leading_anchor, 0, search_end)
                if cut < 0:
                    continue
                cut += len(leading_anchor)
                sequence, quality = sequence[cut:], quality[cut:]
            yield sequence, quality


LEADING_ANCHOR_PREFIX = 15          # bases of a panel tag that identify where the amplicon starts in a read
LEADING_ANCHOR_SAMPLE = 2000        # reads of the first FASTQ inspected by the raw-delivery check
LEADING_ANCHOR_MIN_LENGTH = 8       # a detected anchor must be at least this long
LEADING_ANCHOR_FILES = 8            # FASTQs the check samples, so per-sample barcodes vary and only the shared part survives


def detect_leading_anchor(fastqs: Sequence[Path], tag_prefixes: Set[str], sample_reads: int = LEADING_ANCHOR_SAMPLE) -> Dict[str, object]:
    """Do these FASTQs look like a raw delivery, with a constant sequence in
    front of every amplicon?  A read is *anchored* when a panel tag begins at
    its first base and *shifted* when a panel tag begins further in; the
    bases before a shifted tag are barcode plus common sequence, and the
    anchor is the suffix that at least 90% of those prefixes share.  Reads
    are sampled from several files because an inline barcode is constant
    within a sample and must not be mistaken for the common sequence."""
    anchored = shifted = sampled = 0
    prefixes: List[str] = []
    offsets: Counter = Counter()
    search = LEADING_ANCHOR_PREFIX + LEADING_ANCHOR_SEARCH + 20
    # Spread the sample across the run: neighbouring files often share a
    # plate and part of their barcode, which must not be read as anchor.
    files = list(fastqs)
    if len(files) > LEADING_ANCHOR_FILES:
        step = len(files) / LEADING_ANCHOR_FILES
        files = [files[int(index * step)] for index in range(LEADING_ANCHOR_FILES)]
    per_file = max(1, sample_reads // max(1, len(files)))
    for fastq in files:
        taken = 0
        for sequence, _quality in fastq_records(Path(fastq)):
            sampled += 1
            taken += 1
            if sequence[:LEADING_ANCHOR_PREFIX] in tag_prefixes:
                anchored += 1
            else:
                for offset in range(1, min(len(sequence) - LEADING_ANCHOR_PREFIX, search)):
                    if sequence[offset:offset + LEADING_ANCHOR_PREFIX] in tag_prefixes:
                        shifted += 1
                        prefixes.append(sequence[:offset])
                        offsets[offset] += 1
                        break
            if taken >= per_file:
                break
    anchor = ""
    modal_offset = offsets.most_common(1)[0][0] if offsets else None
    if prefixes:
        # The amplicon starts at the modal offset; a tag found a base or two
        # earlier is a candidate window that happens to begin with the same
        # base as the prefix (chance 1 in 4 per base), so only prefixes of
        # the modal length define the anchor.
        modal_prefixes = [prefix for prefix in prefixes if len(prefix) == modal_offset]
        length = 1
        while modal_prefixes and length <= modal_offset:
            suffix, count = Counter(prefix[-length:] for prefix in modal_prefixes).most_common(1)[0]
            if count < 0.9 * len(modal_prefixes):
                break
            anchor = suffix
            length += 1
    if len(anchor) < LEADING_ANCHOR_MIN_LENGTH:
        anchor = ""
    return {
        "fastqs": [str(fastq) for fastq in files],
        "reads_sampled": sampled,
        "anchored_rate": anchored / max(1, sampled),
        "shifted_rate": shifted / max(1, sampled),
        "offset": modal_offset,
        "anchor": anchor,
    }


def resolve_leading_anchor(args: argparse.Namespace, fastqs: Sequence[Path], tag_prefixes: Set[str], command: str) -> str:
    """Settle ``--leading-anchor`` before any read pass: an explicit sequence
    is used as given; ``auto`` detects it from the first FASTQ; with no
    option the first FASTQ is still checked, and a run that would assign
    almost nothing because the reads carry a constant prefix stops and names
    the anchor to use, rather than producing an empty table."""
    value = getattr(args, "leading_anchor", None)
    if value and value.lower() != "auto":
        return value.upper()
    if not fastqs:
        return ""
    report = detect_leading_anchor([Path(fastq) for fastq in fastqs], tag_prefixes)
    where = Path(fastqs[0]).name + (f" and {len(report['fastqs']) - 1} more" if len(report["fastqs"]) > 1 else "")
    looks_raw = report["shifted_rate"] >= 0.5 and report["anchored_rate"] < 0.05 and report["anchor"]
    if value and value.lower() == "auto":
        if not report["anchor"] or report["shifted_rate"] < 0.5:
            raise MadcError(
                f"{command}: --leading-anchor auto found no constant prefix in {where} "
                f"({report['anchored_rate']:.0%} of {report['reads_sampled']} reads start with a panel tag, "
                f"{report['shifted_rate']:.0%} carry one further in); pass the sequence explicitly or omit the option"
            )
        args.leading_anchor = report["anchor"]
        args.leading_anchor_detected = True
        print(
            f"{command}: leading anchor {report['anchor']} detected in {where} "
            f"(tag start at offset {report['offset']} in {report['shifted_rate']:.0%} of {report['reads_sampled']} reads)",
            file=sys.stderr,
        )
        return report["anchor"]
    if looks_raw:
        raise MadcError(
            f"{command}: {where} looks like a raw delivery: {report['shifted_rate']:.0%} of the first "
            f"{report['reads_sampled']} reads carry a panel tag after a constant prefix ending in {report['anchor']} "
            f"(tag start at offset {report['offset']}), and {report['anchored_rate']:.0%} start with one; rerun with "
            f"--leading-anchor {report['anchor']} (or --leading-anchor auto)"
        )
    return ""


def fastq_records_named(path: Path) -> Iterator[Tuple[str, str, str]]:
    """``(read name, sequence, quality)`` records; the name without the leading @ or comment."""
    with open_text(path, "rt") as handle:
        record_number = 0
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not plus or not quality:
                raise MadcError(f"Truncated FASTQ record {record_number} in {path}")
            if not name.startswith("@") or not plus.startswith("+"):
                raise MadcError(f"Invalid FASTQ record {record_number} in {path}")
            sequence = sequence.strip().upper()
            quality = quality.strip()
            if len(sequence) != len(quality):
                raise MadcError(f"Sequence/quality length mismatch at FASTQ record {record_number} in {path}")
            yield name[1:].split()[0] if name[1:].split() else "", sequence, quality


def mate_id(name: str) -> str:
    """Read name with a trailing /1 or /2 mate suffix removed."""
    return name[:-2] if name.endswith(("/1", "/2")) else name


def _cache_fingerprint(loci: Sequence[Locus], arguments: Mapping[str, object]) -> str:
    # The full loci (tag sequences included) are part of the key: a cache
    # entry is only reusable for the exact panel that produced it, so two
    # window panels over the same FASTQs — e.g. fresh versus frozen windows —
    # can never cross-contaminate even in a shared cache directory.
    payload = {
        "version": CACHE_FORMAT_VERSION,
        "loci": [dataclasses.asdict(locus) for locus in loci],
        "arguments": dict(arguments),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _load_cache(path: Path, sample: Sample, source_stat: os.stat_result) -> Optional[SampleResult]:
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as handle:
            payload = json.load(handle)
        # The fingerprint in the filename covers panel + parameters; size and
        # mtime cover the FASTQ itself, so a re-delivered or edited input
        # invalidates its cache without any content hashing on the hot path.
        if payload["source_size"] != source_stat.st_size or payload["source_mtime_ns"] != source_stat.st_mtime_ns:
            return None
        counts = {(int(marker), sequence): int(count) for marker, sequence, count in payload["counts"]}
        quality_support = {
            (int(marker), sequence): tuple(int(value) for value in values)
            for marker, sequence, *values in payload["quality_support"]
        }
        partial_counts = {
            (int(marker), sequence): int(count)
            for marker, sequence, count in payload["partial_counts"]
        }
        partial_quality_support = {
            (int(marker), sequence): tuple(int(value) for value in values)
            for marker, sequence, *values in payload["partial_quality_support"]
        }
        return SampleResult(
            sample=sample,
            counts=counts,
            quality_support=quality_support,
            partial_counts=partial_counts,
            partial_quality_support=partial_quality_support,
            qc=payload["qc"],
            cache_hit=True,
        )
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
        return None


def _write_cache(path: Path, result: SampleResult, source_stat: os.stat_result) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + f".{os.getpid()}.tmp")
    payload = {
        "source_size": source_stat.st_size,
        "source_mtime_ns": source_stat.st_mtime_ns,
        "qc": result.qc,
        "counts": [[marker, sequence, count] for (marker, sequence), count in result.counts.items()],
        "quality_support": [
            [marker, sequence, *values]
            for (marker, sequence), values in result.quality_support.items()
        ],
        "partial_counts": [
            [marker, sequence, count]
            for (marker, sequence), count in result.partial_counts.items()
        ],
        "partial_quality_support": [
            [marker, sequence, *values]
            for (marker, sequence), values in result.partial_quality_support.items()
        ],
    }
    with gzip.open(temporary, "wt") as handle:
        json.dump(payload, handle, separators=(",", ":"))
    temporary.replace(path)


def select_marker(
    tag: str,
    candidates: Sequence[int],
    loci: Sequence[Locus],
    max_mismatches: int,
    gap_band: int = 0,
) -> Tuple[Optional[int], str]:
    """Resolve primer candidates by identity to the available tag prefix.

    Substitution-only identity decides first and is never overturned.  When
    every candidate rejects the tag and ``gap_band`` is positive, one gapped
    reading is offered under the SAME edit budget (gaps count as edits): a
    start-anchored banded alignment whose end may overhang by up to
    ``gap_band`` bases, so undesigned indel haplotypes — single or multiple
    events — are counted instead of being dropped as off-target.  A gapped
    reading is believed only when it is composed like an indel haplotype
    (few substitutions, few gap events); a co-amplified paralog reaches the
    budget through many substitutions and stays off-target.  A cheap exact
    lower bound (``gapped_prescreen``) skips the alignment for hopeless tags.
    """
    scored = []
    for marker_index in candidates:
        locus = loci[marker_index]
        if tag == locus.ref_tag or tag == locus.alt_tag:
            distance = 0
        else:
            distance = min(
                tag_mismatches(tag, locus.ref_tag[: len(tag)], max_mismatches),
                tag_mismatches(tag, locus.alt_tag[: len(tag)], max_mismatches),
            )
        scored.append((distance, marker_index))
    scored.sort()
    if scored and scored[0][0] <= max_mismatches:
        if len(scored) > 1 and scored[0][0] == scored[1][0]:
            return None, "ambiguous"
        return scored[0][1], "assigned"
    if gap_band <= 0 or not scored:
        return None, "offtarget"
    gapped = []
    for marker_index in candidates:
        locus = loci[marker_index]
        best: Optional[int] = None
        for template in (locus.ref_tag, locus.alt_tag):
            if len(template) != len(tag):
                continue
            if gapped_prescreen(tag, template, gap_band) > max_mismatches:
                continue
            edits, substitutions, gap_bases, gap_events, _, _ = gapped_stats(tag, template, gap_band)
            if edits > max_mismatches or not gapped_explanation_ok(substitutions, gap_bases, gap_events):
                continue
            if best is None or edits < best:
                best = edits
        if best is not None:
            gapped.append((best, marker_index))
    gapped.sort()
    if not gapped:
        return None, "offtarget"
    if len(gapped) > 1 and gapped[0][0] == gapped[1][0]:
        return None, "ambiguous"
    return gapped[0][1], "assigned_gapped"


def resolve_n_tags(
    counts: MutableMapping[Tuple[int, str], int],
    quality_support: MutableMapping[Tuple[int, str], List[int]],
    loci: Sequence[Locus],
    qc: Counter,
) -> None:
    """Reads whose tag carries an N never become rows.

    An N is a missing base call, not a haplotype: a systematically failed
    cycle produces the same N-bearing sequence in every sample and would
    otherwise pass the reproducibility route.  Such reads are credited to
    the one concrete (A/C/G/T-only) sequence at the locus — observed in
    this sample, or a concrete canonical — they are compatible with; when
    none or several qualify they are dropped and counted in
    ``reads_n_dropped``.  An ambiguity-coded canonical is never a
    destination: its row holds exact observations only, and offering it
    would either invent counts for an unobserved template or make a read
    look ambiguous beside the concrete sequence it actually came from.
    """
    n_keys = [key for key in counts if "N" in key[1]]
    if not n_keys:
        return
    concrete_by_marker: Dict[int, List[str]] = defaultdict(list)
    for marker_index, sequence in counts:
        if "N" not in sequence:
            concrete_by_marker[marker_index].append(sequence)
    for marker_index, tag in n_keys:
        locus = loci[marker_index]
        pool = set(concrete_by_marker[marker_index]) | {
            tag_ for tag_ in (locus.ref_tag, locus.alt_tag) if _CONCRETE.issuperset(tag_)
        }
        matches = [s for s in pool if len(s) == len(tag) and compatible_mismatches(tag, s, 0) == 0]
        count = counts.pop((marker_index, tag))
        support = quality_support.pop((marker_index, tag), None)
        if len(matches) == 1:
            destination = (marker_index, matches[0])
            counts[destination] = counts.get(destination, 0) + count
            if support is not None:
                target = quality_support[destination]
                for index, value in enumerate(support):
                    target[index] += value
            qc["reads_n_mapped"] += count
        else:
            qc["reads_n_dropped"] += count


def add_quality_support(
    support: MutableMapping[Tuple[int, str], List[int]],
    key: Tuple[int, str],
    quality: str,
) -> None:
    # min() over the string compares ASCII codes, which is the minimum Phred.
    minimum_quality = ord(min(quality)) - 33
    for quality_index, threshold in enumerate(QUALITY_THRESHOLDS):
        if minimum_quality >= threshold:
            support[key][quality_index] += 1


def process_sample(
    sample: Sample,
    loci: Sequence[Locus],
    primer_index: PrimerIndex,
    settings: AssignmentSettings,
    cache_path: Optional[Path] = None,
) -> SampleResult:
    """Stage 2: stream one sample's FASTQ into exact per-sequence tag counts.

    Each read is routed to locus candidates by its left-primer prefix, then
    accepted for exactly one locus when its tag lies within the identity
    envelope of that locus's Ref/Alt tags, or — when every candidate rejects
    it and ``settings.gap_band`` is positive — under one gapped reading within the
    same edit budget (``select_marker``); ambiguous reads are dropped rather
    than guessed.  Tags carrying an N are resolved at the end of the sample
    (``resolve_n_tags``) and never counted as haplotypes of their own.
    Alongside each count, a per-threshold vector of
    whole-tag-clean read tallies is kept so quality thresholds can change
    later without touching the FASTQs again.  Reads shorter than the tag are
    collected separately for the optional partial-read stage.
    """
    fastq_path = Path(sample.fastq)
    source_stat = fastq_path.stat()
    if cache_path is not None:
        cached = _load_cache(cache_path, sample, source_stat)
        if cached is not None:
            return cached

    qc = Counter()
    exact_counts: Dict[Tuple[int, str], int] = defaultdict(int)
    exact_quality_support: Dict[Tuple[int, str], List[int]] = defaultdict(
        lambda: [0] * len(QUALITY_THRESHOLDS)
    )
    partial_counts: Dict[Tuple[int, str], int] = defaultdict(int)
    partial_quality_support: Dict[Tuple[int, str], List[int]] = defaultdict(
        lambda: [0] * len(QUALITY_THRESHOLDS)
    )
    # Amplicon reads repeat heavily, so the (deterministic) assignment decision
    # is memoized per unique sequence once a primer candidate list was found;
    # primer-less reads fail the trie walk cheaply and are not worth storing.
    full_decisions: Dict[str, Tuple[Optional[int], str]] = {}
    partial_decisions: Dict[str, Tuple[Optional[int], str]] = {}
    for sequence, quality in fastq_records(fastq_path):
        qc["reads_total"] += 1
        if settings.leading_anchor:
            # Raw deliveries that keep the library's common sequence in front
            # of the amplicon: the read starts where that anchor ends.
            cut = sequence.find(settings.leading_anchor, 0, len(settings.leading_anchor) + LEADING_ANCHOR_SEARCH)
            if cut < 0:
                qc["reads_no_leading_anchor"] += 1
                continue
            cut += len(settings.leading_anchor)
            sequence, quality = sequence[cut:], quality[cut:]
        if len(sequence) < settings.tag_length:
            if settings.partial_min_length <= 0 or len(sequence) < settings.partial_min_length:
                qc["reads_too_short"] += 1
                continue
            decision = partial_decisions.get(sequence)
            if decision is None:
                candidates = primer_index.match(sequence, max_mismatches=settings.primer_mismatches)
                if not candidates:
                    qc["reads_partial_no_primer"] += 1
                    continue
                partial_max_mismatches = math.floor(
                    settings.max_offtarget_mismatches * len(sequence) / settings.tag_length
                )
                marker_index, status = select_marker(
                    sequence, candidates, loci, partial_max_mismatches
                )
                decision = (
                    marker_index,
                    "reads_partial" if marker_index is not None else f"reads_partial_{status}",
                )
                partial_decisions[sequence] = decision
            marker_index, qc_key = decision
            qc[qc_key] += 1
            if marker_index is None:
                continue
            key = (marker_index, sequence)
            partial_counts[key] += 1
            add_quality_support(partial_quality_support, key, quality)
            continue
        tag = sequence[:settings.tag_length]
        decision = full_decisions.get(tag)
        if decision is None:
            candidates = primer_index.match(tag, max_mismatches=settings.primer_mismatches)
            if not candidates:
                qc["reads_no_primer"] += 1
                continue
            marker_index, status = select_marker(tag, candidates, loci, settings.max_offtarget_mismatches, settings.gap_band)
            decision = (marker_index, status)
            full_decisions[tag] = decision
        marker_index, status = decision
        if marker_index is None:
            qc[f"reads_{status}"] += 1
            continue
        qc["reads_assigned"] += 1
        if status == "assigned_gapped":
            qc["reads_assigned_gapped"] += 1
        key = (marker_index, tag)
        exact_counts[key] += 1
        add_quality_support(exact_quality_support, key, quality[:settings.tag_length])

    resolve_n_tags(exact_counts, exact_quality_support, loci, qc)
    qc["unique_full_tags_raw"] = len(exact_counts)
    qc["unique_partial_tags_raw"] = len(partial_counts)
    result = SampleResult(
        sample=sample,
        counts=dict(exact_counts),
        quality_support={key: tuple(values) for key, values in exact_quality_support.items()},
        partial_counts=dict(partial_counts),
        partial_quality_support={
            key: tuple(values) for key, values in partial_quality_support.items()
        },
        qc=dict(qc),
        cache_hit=False,
    )
    if cache_path is not None:
        _write_cache(cache_path, result, source_stat)
    return result


def _process_sample_worker(arguments):
    return process_sample(*arguments)


def discover_fastqs(directory: Path, pattern: str) -> List[Path]:
    paths = []
    for path in directory.glob(pattern):
        lower = path.name.lower()
        if path.is_file() and lower.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
            paths.append(path)
    paths.sort(key=lambda path: natural_key(path.name))
    if not paths:
        raise MadcError(f"No FASTQ files matching {pattern!r} found in {directory}")
    return paths


def load_target_metadata(paths: Sequence[Path]) -> Dict[str, Dict[str, str]]:
    metadata: Dict[str, Dict[str, str]] = {}
    for path in paths:
        rows = read_dict_rows(path, delimiter=",")
        if not rows:
            continue
        key = "targetid" if "targetid" in rows[0] else "sample_id" if "sample_id" in rows[0] else None
        if key is None:
            raise MadcError(f"Metadata table {path} needs a targetid or sample_id column")
        for row in rows:
            sample_id = row[key]
            if sample_id in metadata and metadata[sample_id] != row:
                raise MadcError(f"Conflicting metadata for sample {sample_id}")
            metadata[sample_id] = row
    return metadata


def make_samples(
    fastq_paths: Sequence[Path], metadata: Mapping[str, Mapping[str, str]], sample_name_field: str
) -> List[Sample]:
    samples = []
    for fastq in fastq_paths:
        sample_id = strip_fastq_suffix(fastq.name)
        row = dict(metadata.get(sample_id, {}))
        if sample_name_field in ("targetid", "sample_id", "filename"):
            display_name = sample_id
        else:
            display_name = row.get(sample_name_field, "")
            if not display_name:
                raise MadcError(f"No {sample_name_field!r} metadata for FASTQ sample {sample_id}")
        samples.append(Sample(sample_id=sample_id, display_name=display_name, fastq=str(fastq), metadata=row))
    names = [sample.display_name for sample in samples]
    if len(set(names)) != len(names):
        duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
        raise MadcError(
            "Sample column names must be unique; duplicate values from --sample-name-field: "
            + ", ".join(duplicates[:5])
        )
    return samples


def load_sample_lookup(path: Path, sample_ids: Sequence[str]) -> Dict[str, str]:
    """Load a two-column FASTQ-id -> sample-name table for output columns.

    The table is positional (first column: FASTQ-derived sample id, second:
    desired sample name); a first line whose id matches no sample is treated
    as a header and skipped.  Every run sample must be mapped and every
    unknown id is an error — a rename table with holes or typos would
    silently mislabel columns, which is the one mistake a count table can
    never recover from downstream.
    """
    known = set(sample_ids)
    mapping: Dict[str, str] = {}
    unknown: List[str] = []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        sample = handle.read(8192)
        handle.seek(0)
        delimiter = "\t" if "\t" in sample.partition("\n")[0] else ","
        for line_number, row in enumerate(csv.reader(handle, delimiter=delimiter), start=1):
            if not row or not row[0].strip():
                continue
            if len(row) < 2 or not row[1].strip():
                raise MadcError(f"{path} line {line_number}: expected two columns (fastq id, sample name)")
            key, value = row[0].strip(), row[1].strip()
            if key not in known:
                if line_number == 1:
                    continue  # header line
                unknown.append(key)
                continue
            if key in mapping and mapping[key] != value:
                raise MadcError(f"{path} maps sample {key} to two different names")
            mapping[key] = value
    if unknown:
        raise MadcError(f"{path} contains ids matching no FASTQ sample: {', '.join(sorted(unknown)[:5])}")
    missing = sorted(known - set(mapping))
    if missing:
        raise MadcError(
            f"{path} is missing {len(missing)} run sample(s): " + ", ".join(missing[:5])
        )
    duplicates = sorted(name for name, count in Counter(mapping.values()).items() if count > 1)
    if duplicates:
        raise MadcError(f"{path} assigns the same sample name more than once: {', '.join(duplicates[:5])}")
    return mapping


# ============================================================================
# 6. Evidence, denoising and partial-read salvage
# ============================================================================


@dataclasses.dataclass(frozen=True)
class EvidenceRule:
    """Evidence required for a non-canonical sequence to become an output row.

    A sequence must reach ``min_global_count`` project-wide reads and pass at
    least one enabled evidence route:

    - reproducibility route: at least ``reproducible_min_reads`` reads in each
      of at least ``reproducible_min_samples`` distinct samples;
    - depth route: at least ``depth_min_reads`` reads within one sample, so a
      private allele in a single-sample design (for example a diversity panel
      of unreplicated accessions) can still be discovered;
    - quality route: at least ``quality_min_clean_reads`` reads whose every
      tag base is at or above ``discovery_min_quality``.

    A route whose threshold is zero is disabled.  With every route disabled,
    ``min_global_count`` alone gates discovery (exhaustive mode).  The same
    rule decides which sequences denoising must preserve, so discovery,
    protection, and partial-read candidate learning stay consistent.

    ``min_supporting_freq`` optionally requires a sample's reads for the
    sequence to reach that share of the sample's locus total before the
    sample can support the reproducibility or depth route (a
    HaploCharmer-style dosage-aware floor).  Zero disables the floor; note
    that co-amplified paralogs inflate an amplicon pool's denominator, so
    real minor haplotypes can sit at sub-percent frequencies.
    """

    min_global_count: int
    discovery_min_quality: int
    quality_min_clean_reads: int
    reproducible_min_reads: int
    reproducible_min_samples: int
    depth_min_reads: int = 0
    min_supporting_freq: float = 0.0

    @property
    def quality_route_enabled(self) -> bool:
        return self.quality_min_clean_reads > 0

    @property
    def reproducibility_route_enabled(self) -> bool:
        return self.reproducible_min_samples > 0

    @property
    def depth_route_enabled(self) -> bool:
        return self.depth_min_reads > 0

    @property
    def any_route_enabled(self) -> bool:
        return self.quality_route_enabled or self.reproducibility_route_enabled or self.depth_route_enabled

    def evaluate(
        self,
        global_count: int,
        clean_reads: int,
        supporting_samples: int,
        max_sample_count: int,
    ) -> Tuple[bool, str]:
        """Return (passes, route label) for one sequence's aggregated evidence."""
        if global_count < self.min_global_count:
            return False, "below_min_global_count"
        if not self.any_route_enabled:
            return True, "no_route_required"
        passed = []
        if self.reproducibility_route_enabled and supporting_samples >= self.reproducible_min_samples:
            passed.append("reproducibility")
        if self.depth_route_enabled and max_sample_count >= self.depth_min_reads:
            passed.append("depth")
        if self.quality_route_enabled and clean_reads >= self.quality_min_clean_reads:
            passed.append("quality")
        if passed:
            return True, "+".join(passed)
        return False, "failed_evidence"


def collect_evidence(
    per_sample_counts: Sequence[Mapping[Tuple[int, str], int]],
    per_sample_quality: Sequence[Mapping[Tuple[int, str], Sequence[int]]],
    rule: EvidenceRule,
) -> Dict[Tuple[int, str], EvidenceSupport]:
    """Aggregate (global count, clean reads, supporting samples, max sample count) per key."""
    quality_index = QUALITY_THRESHOLDS.index(rule.discovery_min_quality)
    floor = rule.min_supporting_freq
    evidence: Dict[Tuple[int, str], List[int]] = defaultdict(lambda: [0, 0, 0, 0])
    for counts, quality in zip(per_sample_counts, per_sample_quality):
        locus_totals: Dict[int, int] = defaultdict(int)
        if floor > 0.0:
            for (marker_index, _), count in counts.items():
                locus_totals[marker_index] += count
        for key, count in counts.items():
            entry = evidence[key]
            entry[0] += count
            if floor > 0.0 and count < floor * locus_totals[key[0]]:
                # Below the sample's locus-share floor: the global count still
                # accumulates, but this sample cannot support the
                # reproducibility or depth route.
                continue
            if count >= rule.reproducible_min_reads:
                entry[2] += 1
            if count > entry[3]:
                entry[3] = count
        for key, values in quality.items():
            if len(values) > quality_index:
                evidence[key][1] += values[quality_index]
    return {key: EvidenceSupport(entry[0], entry[1], entry[2], entry[3]) for key, entry in evidence.items()}


def denoise_counts(
    counts: Mapping[str, int], max_distance: int, alpha: float, protected: Iterable[str] = ()
) -> Tuple[Dict[str, int], Dict[str, str]]:
    """Abundance-skew denoising inspired by the open UNOISE algorithm.

    A low-abundance sequence is absorbed by a more abundant sequence at Hamming
    distance d when parent/child >= 2 ** (alpha*d + 1).  Among eligible
    parents the closest wins, then the deepest.  Canonical Ref/Alt sequences
    are protected from being absorbed when they are directly observed.
    """
    protected_set = set(protected)
    centroids: Dict[str, int] = {}
    mapping: Dict[str, str] = {}
    # At the default distance of 1 the eligible parents are the sequence's
    # single-substitution neighbours, which are enumerated and looked up
    # rather than found by scanning every centroid: a paralogous locus can
    # carry thousands of unique sequences per sample, and the scan made the
    # pass quadratic in them.  Both paths select the same parent.
    for sequence, abundance in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        if max_distance <= 0 or sequence in protected_set:
            centroids[sequence] = centroids.get(sequence, 0) + abundance
            mapping[sequence] = sequence
            continue
        best_parent: Optional[Tuple[int, int, str]] = None
        minimum_skew = 2 ** (alpha + 1)
        # Scanning a few centroids is cheaper than building 4L neighbours;
        # the lookup pays off once the centroids outnumber the neighbours.
        if max_distance == 1 and len(centroids) > DENOISE_LOOKUP_THRESHOLD * len(sequence):
            for parent in _substitution_neighbours(sequence):
                parent_abundance = centroids.get(parent)
                if parent_abundance is None or parent_abundance < abundance * minimum_skew:
                    continue
                candidate = (1, -parent_abundance, parent)
                if best_parent is None or candidate < best_parent:
                    best_parent = candidate
        else:
            for parent, parent_abundance in centroids.items():
                if parent_abundance < abundance * minimum_skew:
                    continue
                distance = hamming_distance(sequence, parent, stop_after=max_distance)
                if distance == 0 or distance > max_distance:
                    continue
                required_skew = 2 ** (alpha * distance + 1)
                if parent_abundance < abundance * required_skew:
                    continue
                # The CLOSEST parent wins, then the deepest: ranking by depth
                # first could move a distance-1 neighbour of Alt onto a far
                # deeper Ref two edits away — a cross-allele dosage error.
                candidate = (distance, -parent_abundance, parent)
                if best_parent is None or candidate < best_parent:
                    best_parent = candidate
        if best_parent is None:
            centroids[sequence] = abundance
            mapping[sequence] = sequence
        else:
            parent = best_parent[2]
            centroids[parent] += abundance
            mapping[sequence] = parent
    return centroids, mapping


def _substitution_neighbours(sequence: str) -> Iterator[str]:
    """Every sequence at Hamming distance exactly 1 over A, C, G, T and N."""
    for position, base in enumerate(sequence):
        head, tail = sequence[:position], sequence[position + 1:]
        for other in "ACGTN":
            if other != base:
                yield head + other + tail


def supported_sequences(
    sample_results: Sequence[SampleResult],
    loci: Sequence[Locus],
    rule: EvidenceRule,
) -> Dict[int, set[str]]:
    """Return canonical and evidence-passing sequences that denoising must preserve."""
    evidence = collect_evidence(
        [result.counts for result in sample_results],
        [result.quality_support for result in sample_results],
        rule,
    )
    protected: Dict[int, set[str]] = {
        locus.index: {locus.ref_tag, locus.alt_tag} for locus in loci
    }
    for (marker_index, sequence), values in evidence.items():
        if rule.evaluate(*values)[0]:
            protected[marker_index].add(sequence)
    return protected


def _longest_shared_prefix(sequences: Iterable[str]) -> int:
    """Longest prefix shared by any two of ``sequences`` (0 for fewer than two).

    After sorting, the longest prefix shared by any pair equals the longest
    prefix shared by an adjacent pair.

    Follow-up: this compares bases literally, as does the ``startswith`` test
    in ``apply_partial_reads``, so an IUPAC code in a design tag counts as a
    difference and a concrete partial is never a prefix of an ambiguity-coded
    canonical tag (it is unmapped unless a concrete full-length row was
    learned).  That is safe in the never-invent direction but is the one place
    not using the design-match relation (``bases_compatible``) adopted
    elsewhere; both spots should move to it together.
    """
    ordered = sorted(sequences)
    longest = 0
    for left, right in zip(ordered, ordered[1:]):
        shared = 0
        for a, b in zip(left, right):
            if a != b:
                break
            shared += 1
        longest = max(longest, shared)
    return longest


def apply_partial_reads(
    sample_results: Sequence[SampleResult],
    loci: Sequence[Locus],
    mode: str,
    rule: EvidenceRule,
) -> Tuple[List[SampleResult], Dict[str, int]]:
    """Credit trimmed reads to the one full-length allele each is a prefix of.

    Full-length evidence-passing alleles are learned first.  A partial read is
    counted only when it is longer than every prefix shared between two learned
    alleles of its locus and is a prefix of exactly one of them.  The length
    guard is decided per locus, so a locus credits truncated reads for both
    target alleles or for neither and partial assignment can never move an
    allele ratio; because Ref and Alt share their prefix up to the target, a
    partial that does not reach the target base is always ambiguous.  Partial
    reads never establish an allele: a prefix of nothing learned is unmapped.
    (An ``extend`` mode that borrowed the design Ref suffix for such prefixes
    existed through 0.11.0 and was removed after it manufactured a row from a
    primer-dimer-length read.)
    """
    if mode == "off":
        return list(sample_results), {"mode_off": 1}

    candidates = supported_sequences(sample_results, loci, rule)
    shared_prefix = {
        marker_index: _longest_shared_prefix(sequences)
        for marker_index, sequences in candidates.items()
    }
    partial_evidence = collect_evidence(
        [result.partial_counts for result in sample_results],
        [result.partial_quality_support for result in sample_results],
        rule,
    )
    mappings: Dict[Tuple[int, str], Tuple[Optional[str], str]] = {}
    ordered_partials = sorted(
        partial_evidence,
        key=lambda key: (key[0], -len(key[1]), -partial_evidence[key].global_count, key[1]),
    )
    for key in ordered_partials:
        marker_index, partial = key
        if len(partial) <= shared_prefix.get(marker_index, 0):
            mappings[key] = (None, "ambiguous")
            continue
        matches = [sequence for sequence in candidates[marker_index] if sequence.startswith(partial)]
        if len(matches) == 1:
            mappings[key] = (matches[0], "assigned")
            continue
        mappings[key] = (None, "ambiguous" if matches else "unmapped")

    summary = Counter()
    updated = []
    for result in sample_results:
        counts: Dict[Tuple[int, str], int] = defaultdict(int, result.counts)
        quality_support: Dict[Tuple[int, str], List[int]] = defaultdict(
            lambda: [0] * len(QUALITY_THRESHOLDS)
        )
        for key, values in result.quality_support.items():
            quality_support[key] = list(values)
        qc = Counter(result.qc)
        for key, count in result.partial_counts.items():
            destination_sequence, status = mappings[key]
            summary[f"reads_{status}"] += count
            qc[f"partial_reads_{status}"] += count
            if destination_sequence is None:
                continue
            destination = (key[0], destination_sequence)
            counts[destination] += count
            for quality_index_, value in enumerate(result.partial_quality_support.get(key, ())):
                quality_support[destination][quality_index_] += value
        updated.append(
            dataclasses.replace(
                result,
                counts=dict(counts),
                quality_support={key: tuple(values) for key, values in quality_support.items()},
                qc=dict(qc),
            )
        )
    return updated, dict(summary)


def _denoise_sample_result(
    result: SampleResult,
    protected: Mapping[int, set[str]],
    max_distance: int,
    alpha: float,
) -> SampleResult:
    by_marker: Dict[int, Dict[str, int]] = defaultdict(dict)
    for (marker_index, sequence), count in result.counts.items():
        by_marker[marker_index][sequence] = count
    counts: Dict[Tuple[int, str], int] = {}
    quality_support: Dict[Tuple[int, str], List[int]] = defaultdict(
        lambda: [0] * len(QUALITY_THRESHOLDS)
    )
    absorbed_reads = 0
    for marker_index, sequence_counts in by_marker.items():
        centroids, mapping = denoise_counts(
            sequence_counts,
            max_distance=max_distance,
            alpha=alpha,
            protected=protected.get(marker_index, ()),
        )
        for sequence, count in centroids.items():
            counts[(marker_index, sequence)] = count
        for sequence in sequence_counts:
            if mapping[sequence] != sequence:
                absorbed_reads += sequence_counts[sequence]
            destination = (marker_index, mapping[sequence])
            source = result.quality_support.get((marker_index, sequence), ())
            for quality_index, value in enumerate(source):
                quality_support[destination][quality_index] += value
    qc = dict(result.qc)
    qc["unique_tags_before_denoise"] = len(result.counts)
    qc["unique_tags_after_sample_denoise"] = len(counts)
    qc["reads_absorbed_sample_denoise"] = absorbed_reads
    return dataclasses.replace(
        result,
        counts=counts,
        quality_support={key: tuple(values) for key, values in quality_support.items()},
        qc=qc,
    )


_FINALIZE_STATE: Dict[str, object] = {}


def _finalize_init(protected: Mapping[int, set[str]], max_distance: int, alpha: float) -> None:
    """Worker initializer for the per-sample denoise across workers: the
    protected set travels to each worker once, not once per sample."""
    _FINALIZE_STATE.clear()
    _FINALIZE_STATE.update(protected=protected, max_distance=max_distance, alpha=alpha)


def _denoise_sample_worker(result: SampleResult) -> SampleResult:
    return _denoise_sample_result(
        result,
        _FINALIZE_STATE["protected"],  # type: ignore[arg-type]
        max_distance=_FINALIZE_STATE["max_distance"],  # type: ignore[arg-type]
        alpha=_FINALIZE_STATE["alpha"],  # type: ignore[arg-type]
    )


def _pooled_moves(
    locus: Locus,
    sequence_counts: Mapping[str, int],
    protected_sequences: Iterable[str],
    max_distance: int,
    alpha: float,
    gap_band: int,
) -> Dict[str, str]:
    """One locus of the pooled pass: the project-wide skew names a parent for
    every evidence-failing sequence, and the move is kept only when the parent
    is protected or canonical and carries the same target allele."""
    keep = set(protected_sequences) | {locus.ref_tag, locus.alt_tag}
    _, mapping = denoise_counts(sequence_counts, max_distance=max_distance, alpha=alpha, protected=keep)
    groups: Dict[str, str] = {}

    def group(sequence: str) -> str:
        if sequence not in groups:
            groups[sequence] = ALLELE_GROUP.get(classify_sequence(locus, sequence, gap_band), "other")
        return groups[sequence]

    return {
        child: parent
        for child, parent in mapping.items()
        if parent != child and child not in keep and parent in keep and group(child) == group(parent)
    }


_POOLED_STATE: Dict[str, object] = {}


def _pooled_init(
    loci: Sequence[Locus], protected: Mapping[int, set[str]], max_distance: int, alpha: float, gap_band: int
) -> None:
    _POOLED_STATE.clear()
    _POOLED_STATE.update(loci=loci, protected=protected, max_distance=max_distance, alpha=alpha, gap_band=gap_band)


def _pooled_worker(item: Tuple[int, Dict[str, int]]) -> Tuple[int, Dict[str, str]]:
    marker_index, sequence_counts = item
    loci = _POOLED_STATE["loci"]  # type: ignore[assignment]
    protected = _POOLED_STATE["protected"]  # type: ignore[assignment]
    return marker_index, _pooled_moves(
        loci[marker_index], sequence_counts, protected.get(marker_index, ()),  # type: ignore[index,attr-defined]
        _POOLED_STATE["max_distance"], _POOLED_STATE["alpha"], _POOLED_STATE["gap_band"],  # type: ignore[arg-type]
    )


def finalize_counts(
    sample_results: Sequence[SampleResult],
    loci: Sequence[Locus],
    denoise: DenoiseSettings,
    rule: EvidenceRule,
    gap_band: int = GAP_BAND_DEFAULT,
    jobs: int = 1,
) -> List[SampleResult]:
    """Stage 5: denoise counts while preserving every evidence-passing sequence.

    Protection is decided first, from pre-denoise evidence, so a real
    low-abundance allele cannot be absorbed by an abundant neighbor before it
    has had the chance to prove itself.  Abundance-skew denoising then runs
    WITHIN each sample.  Pooled abundance is evidence about the panel, not
    about one sample's genotype: 0.9.x applied the cross-sample mapping to
    every sample regardless of allele, which turned a private 25/25
    heterozygote into Ref=50.  The pooled pass now only folds an
    evidence-failing sequence into the parent the project-wide skew names
    when that parent is protected or canonical AND carries the same target
    allele (Ref/RefMatch, Alt/AltMatch, or Other), so no move can change a
    sample's Ref/Alt dosage while reads that would otherwise be dropped as
    noise keep their depth.  Absorbed reads move into their parent's counts
    (totals are conserved) and the moved amounts are recorded per sample.
    With ``jobs`` > 1 the per-sample pass and the per-locus pooled pass run
    across that many worker processes (0.12.1); protection, the final
    bookkeeping and the table write stay in this process.
    """
    stage_clock = time.perf_counter()

    def _stage(label: str) -> None:
        nonlocal stage_clock
        now = time.perf_counter()
        print(f"finalize: {label} in {now - stage_clock:,.0f} s", file=sys.stderr)
        stage_clock = now

    protected = supported_sequences(sample_results, loci, rule)
    _stage("evidence protection")
    parallel_samples = jobs > 1 and len(sample_results) > 1
    if parallel_samples:
        # Each sample's pass depends only on its own counts and the shared
        # protected set, so the samples are denoised across the workers;
        # map keeps their order and the result equals the sequential pass.
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=jobs, initializer=_finalize_init,
            initargs=(protected, denoise.sample_distance, denoise.alpha),
        ) as executor:
            sample_denoised = list(executor.map(_denoise_sample_worker, sample_results))
    else:
        sample_denoised = [
            _denoise_sample_result(
                result,
                protected,
                max_distance=denoise.sample_distance,
                alpha=denoise.alpha,
            )
            for result in sample_results
        ]
    _stage(f"per-sample denoise ({len(sample_results)} samples, jobs={jobs if parallel_samples else 1})")
    aggregate: Dict[int, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    for result in sample_denoised:
        for (marker_index, sequence), count in result.counts.items():
            aggregate[marker_index][sequence] += count
    moves: Dict[int, Dict[str, str]] = {}
    if jobs > 1 and len(aggregate) > 1:
        # The pooled pass is independent per locus (its aggregate, the
        # protected set and the gapped classifier), so the loci are spread
        # over the workers; the panel and protected set travel once per worker.
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=jobs, initializer=_pooled_init,
            initargs=(loci, protected, denoise.global_distance, denoise.alpha, gap_band),
        ) as executor:
            moves = dict(executor.map(_pooled_worker, aggregate.items(), chunksize=16))
    else:
        for marker_index, sequence_counts in aggregate.items():
            moves[marker_index] = _pooled_moves(
                loci[marker_index], sequence_counts, protected.get(marker_index, ()),
                denoise.global_distance, denoise.alpha, gap_band,
            )
    finalized = []
    for result in sample_denoised:
        counts: Dict[Tuple[int, str], int] = defaultdict(int)
        quality_support: Dict[Tuple[int, str], List[int]] = defaultdict(
            lambda: [0] * len(QUALITY_THRESHOLDS)
        )
        absorbed_reads = 0
        for (marker_index, sequence), count in result.counts.items():
            destination_sequence = moves.get(marker_index, {}).get(sequence, sequence)
            if destination_sequence != sequence:
                absorbed_reads += count
            destination = (marker_index, destination_sequence)
            counts[destination] += count
            source_support = result.quality_support.get((marker_index, sequence), ())
            for quality_index, value in enumerate(source_support):
                quality_support[destination][quality_index] += value
        qc = dict(result.qc)
        qc["unique_tags_after_global_denoise"] = len(counts)
        qc["reads_absorbed_global_denoise"] = absorbed_reads
        finalized.append(
            dataclasses.replace(
                result,
                counts=dict(counts),
                quality_support={key: tuple(values) for key, values in quality_support.items()},
                qc=qc,
            )
        )
    _stage("pooled denoise, moves and per-sample bookkeeping")
    return finalized


# ============================================================================
# 7. Classification and annotation (one reading per row)
# ============================================================================


def _reading(sequence: str, template: str, gap_band: int, guard: bool = True) -> Optional[GappedReading]:
    """The gapped reading of ``sequence`` against ``template``, or ``None`` for a
    substitution-only row (or unequal lengths).

    With ``guard`` (observed rows) this is accepted_gapped_reading, the same
    decision classify_sequence uses, so a row's Variants and its category
    always come from one alignment.  With ``guard=False`` (a design's Alt
    against its Ref, knowledge rather than an observation) it is the
    cheapest explanation whatever its composition.
    """
    if len(template) != len(sequence):
        return None
    substitutions = tag_mismatches(sequence, template)
    if guard:
        return accepted_gapped_reading(sequence, template, gap_band, substitutions)
    if gap_band <= 0 or substitutions <= 1 or gapped_prescreen(sequence, template, gap_band) >= substitutions:
        return None
    ops, end_read, end_ref, edits = gapped_alignment(sequence, template, gap_band)
    if edits >= substitutions:
        return None
    return _reading_from_ops(sequence, template, ops, end_read, end_ref, edits)


def marker_type_label(locus: Locus, sequence: str, gap_band: int = GAP_BAND_DEFAULT, events=None) -> str:
    """What a row carries beyond its canonical: '.', 'snp', 'indel' or 'both'.

    Derived from the row's Variants (variant_events, always against the Ref
    tag): the Alt row is labelled by the design itself ('snp' for a
    designed SNP, 'indel' for a designed indel), every other row by its
    events once the marker's designed events are removed, so an AltMatch of
    a designed indel that adds one SNP is 'snp' while its Variants list the
    inherited indel and the SNP.  'snp' = substitutions only, 'indel' =
    gaps only, 'both' = gaps plus substitutions, '.' = no concrete event
    beyond the canonical (the Ref row; a row that merely realises an
    ambiguity code).  Downstream users can keep or drop undesigned indel
    haplotypes by this label alone.  ``events`` may be passed when already
    computed.
    """
    if sequence == locus.ref_tag:
        return "."
    if events is None:
        events = variant_events(locus, sequence, gap_band)
    if sequence != locus.alt_tag:
        designed = Counter(_designed_events(locus.ref_tag, locus.alt_tag, gap_band))
        residual = []
        for event in events:
            if designed[event]:
                designed[event] -= 1
            else:
                residual.append(event)
        events = residual
    kinds = {event[0] for event in events if event[0] in ("snp", "ins", "del")}
    if not kinds:
        return "."
    if kinds == {"snp"}:
        return "snp"
    return "both" if "snp" in kinds else "indel"


def _describe(ref: str, sequence: str, gap_band: int, gapped: bool):
    """Events of ``sequence`` against ``ref`` (the core of variant_events).

    ``gapped`` is the caller's decision — row_reading for an observed row,
    the unguarded cheapest reading for a design — so the description never
    re-decides what classification decided."""
    if sequence == ref or len(sequence) != len(ref):
        return []
    if gapped:
        ops, end_read, end_ref, _ = gapped_alignment(sequence, ref, gap_band)
    else:
        return [
            ("snp", index + 1, ref[index], sequence[index])
            for index in range(len(ref))
            if not _design_matches(sequence[index], ref[index])
        ]
    ops = _left_align(ops, sequence, ref)
    events = []
    k = 0
    while k < len(ops):
        kind, i, j = ops[k]
        if kind == "M":
            if not _design_matches(sequence[i], ref[j]):
                events.append(("snp", j + 1, ref[j], sequence[i]))
            k += 1
            continue
        run_end = k
        while run_end + 1 < len(ops) and ops[run_end + 1][0] == kind:
            run_end += 1
        if kind == "I":
            inserted = sequence[i : ops[run_end][1] + 1]
            if j > 0:
                events.append(("ins", j, ref[j - 1], ref[j - 1] + inserted))
            else:
                events.append(("ins", 1, ref[0], inserted + ref[0]))
        else:
            deleted = ref[j : ops[run_end][2] + 1]
            if j > 0:
                events.append(("del", j, ref[j - 1] + deleted, ref[j - 1]))
            else:
                following = ops[run_end][2] + 1
                events.append(("del", 1, deleted + ref[following], ref[following]))
        k = run_end + 1
    if end_ref < len(ref):
        events.append(("clip", end_ref + 1, ref[end_ref:], "."))
    if end_read < len(sequence):
        events.append(("overhang", len(ref) + 1, ".", sequence[end_read:]))
    return events


@functools.lru_cache(maxsize=4096)
def _designed_events(ref: str, alt: str, gap_band: int):
    """The marker's own design against Ref: the Alt tag's events (cached per locus)."""
    return tuple(_describe(ref, alt, gap_band, gapped=_reading(alt, ref, gap_band, guard=False) is not None))


def variant_events(locus: Locus, sequence: str, gap_band: int = GAP_BAND_DEFAULT):
    """Every event a row carries against the Ref tag, in 1-based tag coordinates.

    Returns ``(kind, position, ref_allele, alt_allele)`` tuples in tag order:
    ``snp`` (``91, A, T``), ``ins`` and ``del`` anchored VCF-style on the
    preceding reference base (``55, C, CAAAA``; ``70, TAG, T``) and left
    normalised, then the window geometry a fixed-length row cannot show:
    ``clip`` (the reference bases an insertion pushed past the window, with
    their first position) and ``overhang`` (the downstream bases a deletion
    pulled in, positioned one past the tag).  The frame is always Ref, so an
    AltMatch row lists the designed substitution, and every row of a marker
    whose own Alt differs from Ref by a gap is read with gaps (the inherited
    designed indel appears in its list rather than as a run of shifted
    substitutions).  Positions are on the tag as printed (design strand); a
    consumer that knows the strand flips them, as it does for substitutions
    already.
    """
    designed = _designed_events(locus.ref_tag, locus.alt_tag, gap_band)
    force_gapped = any(kind in ("ins", "del") for kind, *_ in designed)
    gapped = force_gapped or row_reading(locus, sequence, gap_band) is not None
    return _describe(locus.ref_tag, sequence, gap_band, gapped)


def format_variants(events) -> str:
    """The Variants cell: ``;``-joined tokens, ``pos:REF>ALT`` per event plus ``clip=k`` / ``overhang=k``; '.' when empty."""
    tokens = []
    for kind, position, ref_allele, alt_allele in events:
        if kind == "clip":
            tokens.append(f"clip={len(ref_allele)}")
        elif kind == "overhang":
            tokens.append(f"overhang={len(alt_allele)}")
        else:
            tokens.append(f"{position}:{ref_allele}>{alt_allele}")
    return ";".join(tokens) or "."


def classify_sequence(locus: Locus, sequence: str, gap_band: int = GAP_BAND_DEFAULT) -> str:
    """Classify an observed tag against the canonical Ref/Alt definitions.

    The Ref/Alt categories require exact sequence equality; per-row read
    counts are always exact-sequence counts, so an ambiguity-coded canonical
    tag never accumulates concrete observations.  The target base is read at
    the position the CHEAPER explanation of the tag puts it: the fixed
    offset when a substitution-only reading is as good as any gapped one, or
    the aligned position of the one accepted gapped reading
    (accepted_gapped_reading: fewest gap openings, strictly cheaper, well
    composed, verified by GAP_VERIFIED_SUFFIX bases after its last gap —
    clustered errors at the read end must not masquerade as a shift).  That
    reading is the same one variant_events describes, so a row's category
    and its Variants never disagree about where the target sits.  A target inside a deletion, or
    pushed past the read, is missing data and classifies as Other.  The
    variant-base comparison is IUPAC-aware on the design side only: an
    ambiguous observed base such as N is missing data, not evidence for
    either allele, and classifies as Other.
    """
    if sequence == locus.ref_tag:
        return "Ref"
    if sequence == locus.alt_tag:
        return "Alt"
    index = locus.variant_pos
    if gap_band > 0:
        # The cheapest substitution-only reading over both canonicals is the
        # bound a gapped reading must beat; the reading itself is the one
        # accepted_gapped_reading gives, shared with the row's Variants.
        best = row_reading(locus, sequence, gap_band)
        if best is not None:
            mapped = best.ref_to_read.get(locus.variant_pos)
            if mapped is None or mapped >= len(sequence):
                return "Other"
            index = mapped
    base = sequence[index]
    if len(IUPAC_BASES.get(base, frozenset(base))) != 1:
        return "Other"
    if bases_compatible(base, locus.ref_base):
        return "RefMatch"
    if bases_compatible(base, locus.alt_base):
        return "AltMatch"
    return "Other"


# ============================================================================
# 8. MADC input/output and QC
# ============================================================================


def build_dart_preamble(
    samples: Sequence[Sample], project: str, provenance: str = ""
) -> List[List[str]]:
    """Vendor-style MADC preamble: seven starred rows above the header.

    Replicates the raw-MADC layout so preamble-expecting consumers (HapApp,
    HapApp_utils) accept generated files: per-sample cells carry the
    service/project label, plate barcode, plate number, plate row, plate
    column, plate name, and finally the targetid row, which sits directly
    above the header row's sample names.  The targetid row carries the
    FASTQ-derived sample ids, preserving the fastq-id-above-sample-name
    adjacency of the vendor format.  Every cell we have no value for is
    '*', exactly as in vendor files — including the plate cells we cannot
    know, which are unknowable rather than blank.

    The single exception is the file's very first cell, which carries the
    whole ``provenance`` string (pipeline, version, parameters) so the file
    records how it was made even when the QC sidecar gets separated from
    it.  The string is prefixed with ``* `` so consumers that skip preamble
    lines by their leading star (HapApp_utils parsers do) still skip it.
    """
    width = len(MADC_METADATA_COLUMNS)
    stars = ["*"] * width
    star_cells = ["*"] * len(samples)
    return [
        [f"* {provenance}" if provenance else "*"] + ["*"] * (width - 1) + [project] * len(samples),
        stars + star_cells,
        stars + star_cells,
        stars + star_cells,
        stars + star_cells,
        stars + star_cells,
        stars + [sample.sample_id for sample in samples],
    ]


def is_madc_header(row: Sequence[str]) -> bool:
    """The header row names AlleleID: first in vendor and MADCATS files, after
    a leading Code_version column in HapApp's fixed-allele-ID output."""
    return any(cell.strip() == "AlleleID" for cell in row)


def find_madc_header(rows: Sequence[Sequence[str]]) -> int:
    for index, row in enumerate(rows):
        if is_madc_header(row):
            return index
    raise MadcError("Could not find an AlleleID header row in MADC file")


def madc_identity_columns(header: Sequence[str], path: Path) -> Tuple[int, int, int]:
    """Indices of AlleleID, CloneID and AlleleSequence in a MADC header row."""
    try:
        return header.index("AlleleID"), header.index("CloneID"), header.index("AlleleSequence")
    except ValueError as exc:
        raise MadcError(f"MADC is missing a required identifier column: {path}") from exc


def madc_sample_start(header: Sequence[str]) -> int:
    """First per-sample column: after the statistics block, or after the sequence columns of a short layout."""
    if "RatioAvgCountRefAvgCountSnp" in header:
        return header.index("RatioAvgCountRefAvgCountSnp") + 1
    if "ClusterConsensusSequence" in header:
        return header.index("ClusterConsensusSequence") + 1
    return header.index("AlleleSequence") + 1


def madc_read_counts(cells: Sequence[str], path: Path) -> List[int]:
    """A MADC row's per-sample counts.  Blank cells are 0; HapApp writes counts of merged rows as floats (12.0)."""
    try:
        return list(map(int, cells))
    except ValueError:
        counts = []
        for value in cells:
            try:
                counts.append(int(float(value or 0)))
            except ValueError as exc:
                raise MadcError(f"Non-numeric read count {value!r} in {path}") from exc
        return counts


class CountOverflowCheck:
    """Signed 8-bit overflow in a MADC's counts.

    Some vendor MADCs store counts in a signed byte: counts of 128-255 come
    out as count - 256 (negative) and counts of 256 or more wrap to small
    positives with no trace.  A negative read count is impossible, so any is
    reported, and a maximum of exactly 127 beside negatives points to overflow.
    Adding 256 to the negatives would recover only the 128-255 band, so the
    counts are reported, never corrected.
    """

    def __init__(self) -> None:
        self.negative_cells = 0
        self.markers: Set[str] = set()
        self.minimum: Optional[int] = None
        self.maximum: Optional[int] = None

    def add(self, clone_id: str, counts: Sequence[int]) -> None:
        if not counts:
            return
        low, high = min(counts), max(counts)
        self.minimum = low if self.minimum is None else min(self.minimum, low)
        self.maximum = high if self.maximum is None else max(self.maximum, high)
        if low < 0:
            self.negative_cells += sum(value < 0 for value in counts)
            self.markers.add(clone_id)

    def summary(self) -> Dict[str, object]:
        return {
            "negative_cells": self.negative_cells,
            "markers_with_negative_cells": len(self.markers),
            "minimum": self.minimum,
            "maximum": self.maximum,
            "signed_8bit_signature": bool(self.negative_cells) and self.maximum == 127,
        }

    def warning(self, path: Path) -> Optional[str]:
        if not self.negative_cells:
            return None
        message = (
            f"warning: {path} has {self.negative_cells:,} negative read counts across {len(self.markers):,} "
            f"markers (minimum {self.minimum})"
        )
        if self.maximum == 127:
            message += (
                ", and its largest count is exactly 127: the counts look like signed 8-bit overflow, where "
                "128-255 are stored as count - 256 and counts of 256 or more wrap silently"
            )
        return message + ". Counts are used as stored; regenerate them from the FASTQs or ask the vendor."


def read_madc(path: Path):
    with open_text(path, "rt") as handle:
        rows = list(csv.reader(handle))
    header_index = find_madc_header(rows)
    header = rows[header_index]
    allele_index, clone_index, sequence_index = madc_identity_columns(header, path)
    sample_start = madc_sample_start(header)
    sample_names = list(header[sample_start:])
    for preamble in reversed(rows[:header_index]):
        values = [value.strip() for value in preamble[sample_start:]]
        if len(values) != len(sample_names) or not values:
            continue
        # A numeric row of unique values is a targetid row; plate, replicate,
        # and service-number rows are numeric too but repeat values.
        if len(set(values)) != len(values):
            continue
        numeric = sum(bool(re.fullmatch(r"\d+", value)) for value in values)
        if numeric >= max(1, math.ceil(0.8 * len(values))):
            sample_names = values
            break
    records = []
    overflow = CountOverflowCheck()
    for row in rows[header_index + 1 :]:
        if not row or not row[allele_index].strip():
            continue
        padded = list(row) + [""] * max(0, len(header) - len(row))
        counts = madc_read_counts(padded[sample_start : sample_start + len(sample_names)], path)
        overflow.add(padded[clone_index], counts)
        records.append(
            {
                "allele_id": padded[allele_index],
                "clone_id": padded[clone_index],
                "sequence": padded[sequence_index].upper(),
                "counts": counts,
            }
        )
    warning = overflow.warning(path)
    if warning:
        print(warning, file=sys.stderr)
    return {
        "samples": sample_names, "records": records, "header_index": header_index, "count_overflow": overflow.summary()
    }


def allele_type(allele_id: str) -> str:
    value = allele_id.split("|", 1)[1] if "|" in allele_id else allele_id
    value = re.sub(r"_?\d+$", "", value)
    # Match categories must precede Ref/Alt because "RefMatch" starts with
    # "Ref" (and likewise for AltMatch).
    for category in ("RefMatch", "AltMatch", "Other", "Ref", "Alt"):
        if value.startswith(category):
            return category
    return value


def write_madc(
    output: Path,
    loci: Sequence[Locus],
    results: Sequence[SampleResult],
    rule: EvidenceRule,
    options: OutputOptions = OutputOptions(),
    gap_band: int = GAP_BAND_DEFAULT,
) -> Dict[str, object]:
    """Stage 6: emit the MADC table and account for every counted read.

    Both canonical rows are forced for every locus, at zero counts if
    necessary, so the file always presents the designed biallelic axis; a
    zero-count Ref beside a read-bearing RefMatch is the expected signature
    of germplasm that carries the target base on a diverged haplotype.
    Novel sequences must pass the evidence rule, with the deciding route
    recorded per row class.  Vendor conventions are followed: canonical rows
    use their own sequence as ClusterConsensusSequence, Match/Other rows use
    the Ref tag, CallRate is canonical sample coverage, and AvgCount fields
    are nonzero means.  Statistics columns the pipeline does not compute are
    '.' rather than blank.  Reads in sequences that pass no route are
    tallied in the exclusion buckets so the table's coverage of assigned
    reads is always visible.
    """
    output.parent.mkdir(parents=True, exist_ok=True)
    observed: Dict[int, Counter] = defaultdict(Counter)
    for result in results:
        for (marker_index, sequence), count in result.counts.items():
            observed[marker_index][sequence] += count
    evidence = collect_evidence(
        [result.counts for result in results],
        [result.quality_support for result in results],
        rule,
    )
    category_counts = Counter()
    marker_type_counts: Counter = Counter()
    variant_counts: Counter = Counter()
    rows_with_variants = 0
    sidecar_rows: List[List[object]] = []
    route_counts = Counter()
    rows_written = 0
    excluded = {
        "below_min_global_count": {"rows": 0, "reads": 0},
        "failed_evidence": {"rows": 0, "reads": 0},
    }
    reads_in_emitted_rows = 0
    per_sample_emitted = [0] * len(results)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        for row in options.preamble or ():
            writer.writerow(row)
        header_columns = [
            ANNOTATION_COLUMNS.get(column, column) if options.annotations else column
            for column in MADC_METADATA_COLUMNS
        ]
        writer.writerow(header_columns + [result.sample.display_name for result in results])
        for locus in loci:
            forced = {locus.ref_tag: "Ref", locus.alt_tag: "Alt"}
            sequences = set(forced)
            for sequence, count in observed.get(locus.index, {}).items():
                if sequence in forced:
                    continue
                passed, route = rule.evaluate(*evidence[(locus.index, sequence)])
                if passed:
                    sequences.add(sequence)
                    route_counts[route] += 1
                else:
                    excluded[route]["rows"] += 1
                    excluded[route]["reads"] += count
            rows = []
            for sequence in sequences:
                category = forced.get(sequence) or classify_sequence(locus, sequence, gap_band)
                rows.append((ALLELE_TYPE_ORDER.get(category, 99), category, sequence))
            rows.sort(key=lambda row: (row[0], row[2]))
            type_serials = Counter()
            ref_counts = [result.counts.get((locus.index, locus.ref_tag), 0) for result in results]
            alt_counts = [result.counts.get((locus.index, locus.alt_tag), 0) for result in results]
            canonical_call_rate = sum((ref + alt) > 0 for ref, alt in zip(ref_counts, alt_counts)) / len(results)
            avg_ref = mean_nonzero(ref_counts)
            avg_alt = mean_nonzero(alt_counts)
            for _, category, sequence in rows:
                counts = [result.counts.get((locus.index, sequence), 0) for result in results]
                type_serials[category] += 1
                suffix = category
                if options.number_match_ids and category in ("RefMatch", "AltMatch", "Other"):
                    suffix += f"_{type_serials[category]:04d}"
                allele_id = f"{locus.marker}|{suffix}"
                if category in ("Ref", "Alt"):
                    call_rate = format_number(canonical_call_rate)
                    avg_count_ref = format_number(avg_ref)
                    avg_count_snp = format_number(avg_alt)
                    consensus = sequence
                else:
                    call_rate = "0"
                    average = mean_nonzero(counts)
                    avg_count_ref = format_number(average)
                    avg_count_snp = format_number(average)
                    consensus = locus.ref_tag
                # Statistics columns the pipeline does not compute emit '.'
                # (missing value), never a blank cell.  With ``annotations``
                # two never-computed columns carry the row's receipt instead:
                # PICSnp its Variants (events against the Ref tag) and AvgPIC
                # its MarkerType label; the sidecar gets one line per event.
                label = variants = "."
                if options.annotations:
                    events = variant_events(locus, sequence, gap_band)
                    label = marker_type_label(locus, sequence, gap_band, events)
                    marker_type_counts[label] += 1
                    variants = format_variants(events)
                    if events:
                        rows_with_variants += 1
                    for kind, position, ref_allele, alt_allele in events:
                        variant_counts[kind] += 1
                        sidecar_rows.append([locus.marker, allele_id, sequence, position, ref_allele, alt_allele, kind])
                metadata = [
                    allele_id,
                    locus.marker,
                    sequence,
                    consensus,
                    call_rate,
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    ".",
                    variants,
                    label,
                    avg_count_ref,
                    avg_count_snp,
                    ".",
                ]
                writer.writerow(metadata + counts)
                rows_written += 1
                category_counts[category] += 1
                reads_in_emitted_rows += sum(counts)
                for sample_index, value in enumerate(counts):
                    per_sample_emitted[sample_index] += value
    if options.annotations and options.variants_path is not None:
        options.variants_path.parent.mkdir(parents=True, exist_ok=True)
        with options.variants_path.open("w", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(VARIANTS_SIDECAR_COLUMNS)
            writer.writerows(sidecar_rows)
    elif options.variants_path is not None and options.variants_path.exists():
        options.variants_path.unlink()  # a sidecar from an earlier run must not outlive the table it described
    return {
        "rows": rows_written,
        "categories": dict(category_counts),
        "marker_types": dict(marker_type_counts),
        "variants": {
            "rows_with_events": rows_with_variants,
            "events": dict(variant_counts),
            "sidecar": str(options.variants_path) if (options.annotations and options.variants_path is not None) else None,
        },
        "novel_rows_by_evidence_route": dict(route_counts),
        "excluded_novel": excluded,
        "reads_in_emitted_rows": reads_in_emitted_rows,
        "per_sample_reads_in_emitted_rows": per_sample_emitted,
    }


def validate_madc(generated_path: Path, expected_path: Path) -> Dict[str, object]:
    generated = read_madc(generated_path)
    expected = read_madc(expected_path)
    generated_samples = {name: index for index, name in enumerate(generated["samples"])}
    expected_samples = {name: index for index, name in enumerate(expected["samples"])}
    common_samples = sorted(set(generated_samples) & set(expected_samples), key=natural_key)
    generated_rows = {(row["clone_id"], row["sequence"]): row for row in generated["records"]}
    expected_rows = {(row["clone_id"], row["sequence"]): row for row in expected["records"]}
    common_rows = set(generated_rows) & set(expected_rows)
    generated_only = set(generated_rows) - set(expected_rows)
    expected_only = set(expected_rows) - set(generated_rows)
    classification_matches = sum(
        allele_type(generated_rows[key]["allele_id"]) == allele_type(expected_rows[key]["allele_id"])
        for key in common_rows
    )

    n = exact_cells = 0
    absolute_error = sx = sy = sxx = syy = sxy = 0.0
    union_n = union_exact = 0
    union_absolute_error = union_sx = union_sy = union_sxx = union_syy = union_sxy = 0.0
    generated_total = expected_total = 0
    for key in common_rows:
        generated_row = generated_rows[key]
        expected_row = expected_rows[key]
        for sample in common_samples:
            x = generated_row["counts"][generated_samples[sample]]
            y = expected_row["counts"][expected_samples[sample]]
            n += 1
            exact_cells += x == y
            absolute_error += abs(x - y)
            sx += x
            sy += y
            sxx += x * x
            syy += y * y
            sxy += x * y
            generated_total += x
            expected_total += y
            if x > 0 or y > 0:
                union_n += 1
                union_exact += x == y
                union_absolute_error += abs(x - y)
                union_sx += x
                union_sy += y
                union_sxx += x * x
                union_syy += y * y
                union_sxy += x * y

    expected_category_counts = Counter(allele_type(row["allele_id"]) for row in expected["records"])
    generated_category_counts = Counter(allele_type(row["allele_id"]) for row in generated["records"])
    category_overlap = {}
    for category in ALLELE_TYPE_ORDER:
        generated_keys = {
            key for key, row in generated_rows.items() if allele_type(row["allele_id"]) == category
        }
        expected_keys = {key for key, row in expected_rows.items() if allele_type(row["allele_id"]) == category}
        overlap = len(generated_keys & expected_keys)
        category_overlap[category] = {
            "common": overlap,
            "generated": len(generated_keys),
            "expected": len(expected_keys),
            "precision": overlap / len(generated_keys) if generated_keys else None,
            "recall": overlap / len(expected_keys) if expected_keys else None,
        }
    report = {
        "generated": str(generated_path.resolve()),
        "expected": str(expected_path.resolve()),
        "sample_overlap": {
            "common": len(common_samples),
            "generated": len(generated_samples),
            "expected": len(expected_samples),
        },
        "count_overflow": {"generated": generated["count_overflow"], "expected": expected["count_overflow"]},
        "allele_rows": {
            "common": len(common_rows),
            "generated": len(generated_rows),
            "expected": len(expected_rows),
            "precision": len(common_rows) / len(generated_rows) if generated_rows else 0.0,
            "recall": len(common_rows) / len(expected_rows) if expected_rows else 0.0,
            "generated_only": len(generated_only),
            "expected_only": len(expected_only),
            "classification_agreement": classification_matches / len(common_rows) if common_rows else None,
            "generated_categories": dict(generated_category_counts),
            "expected_categories": dict(expected_category_counts),
            "by_category": category_overlap,
        },
        "common_row_count_cells": {
            "cells": n,
            "exact_fraction": exact_cells / n if n else None,
            "mean_absolute_error": absolute_error / n if n else None,
            "pearson_r": pearson_from_sums(n, sx, sy, sxx, syy, sxy),
            "generated_total": generated_total,
            "expected_total": expected_total,
            "total_ratio": generated_total / expected_total if expected_total else None,
        },
        "nonzero_union_count_cells": {
            "cells": union_n,
            "exact_fraction": union_exact / union_n if union_n else None,
            "mean_absolute_error": union_absolute_error / union_n if union_n else None,
            "pearson_r": pearson_from_sums(
                union_n, union_sx, union_sy, union_sxx, union_syy, union_sxy
            ),
        },
    }
    return report


# ============================================================================
# 9. Commands and command line
# ============================================================================


def run_merge_pairs(args: argparse.Namespace) -> int:
    """Stitch overlapping paired-end reads into single fragment reads.

    The paired-end analog of normalize-reads' role for AgriSeq: platforms
    whose insert is shorter than the summed read lengths (GT-seq PE100 on
    ~130 bp amplicons, read-through capture) carry each fragment twice with
    an overlap; merging reconstructs the fragment once, with read-through
    adapter tails removed by the engulfed-overlap geometry.  Pairs with no
    acceptable overlap (true long inserts) fall back per --unmerged: keep
    the R1 (default), keep both mates (pool — for dual-anchored designs
    such as Flex-Seq whose mates cover opposite locus halves), or drop.
    Samples are independent, so --jobs merges them in parallel processes.
    """
    fastq_dir = Path(args.fastq_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    r1_files = sorted(fastq_dir.glob(f"*{args.r1_suffix}"))
    if not r1_files:
        raise MadcError(f"No files matching *{args.r1_suffix} in {fastq_dir}")
    tasks = []
    for r1_path in r1_files:
        sample = r1_path.name[: -len(args.r1_suffix)]
        r2_path = fastq_dir / f"{sample}{args.r2_suffix}"
        if not r2_path.exists():
            raise MadcError(f"Missing mate file for {sample}: {r2_path.name}")
        tasks.append(
            (str(r1_path), str(r2_path), str(out_dir / f"{sample}.fastq.gz"),
             sample, args.min_overlap, args.max_mismatch_rate, args.unmerged)
        )
    totals: Counter = Counter()
    if args.jobs == 1:
        results = [_merge_pairs_sample(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.jobs) as executor:
            results = list(executor.map(_merge_pairs_sample, tasks))
    for _sample, stats in results:
        totals.update(stats)
    grand = sum(totals.values())
    print(
        f"merge-pairs: {totals['pairs_merged']:,}/{grand:,} pairs merged "
        f"({totals['pairs_merged']/max(1, grand):.1%}) across {len(r1_files)} samples -> {out_dir}",
        file=sys.stderr,
    )
    return 0


def run_window_scan(args: argparse.Namespace) -> int:
    """Sweep candidate tag lengths and report the three selection curves.

    Window length trades three measurable quantities, and until now the
    trade was settled per platform against external truth (vendor calls or
    a reference VCF).  This diagnostic computes the same decision from the
    project's own reads:

    - retention: fraction of located reads fully covering the best window
      at each length.  Predicted for a dense length grid from ONE coverage
      pass (read span classes make it arithmetic), and measured exactly for
      each candidate length by a real normalize pass.
    - cleanliness: of assigned reads, the fraction NOT absorbed by
      denoising and NOT dropped as failed evidence — per-base error and
      indel exposure accumulate with length, so this falls (and cliffs on
      homopolymer chemistries).
    - information: evidence-passing novel (Match) rows per 100 markers —
      microhaplotype diversity captured; rises with length until the local
      linked variation is exhausted, then saturates.

    Information leads the eye but cannot lead the decision: systematic
    errors (Ion homopolymer indel families recur across samples) pass the
    reproducibility route and inflate the information curve with false
    Match rows, and an information-weighted score follows them long —
    measured on Ion oyster data, where such a score recommended 81 bp
    against the externally validated 41-bp knee.  The recommendation is
    therefore SURVIVAL-FIRST: survival = retention x assignment x
    cleanliness; candidates within ``--info-band`` of the best survival
    are eligible, and the eligible candidate with the most information
    wins (ties break short).  On the four benchmarked platforms this rule
    reproduces every externally validated choice.  The report also names
    the information-saturation length (smallest candidate reaching >= 95%
    of maximum information).  Each candidate runs the real
    normalize-reads -> generate machinery on a sample subset
    (``--scan-samples``: the first files in sorted order, or a random draw
    with ``--scan-seed``; ``--generate-args`` forwards evidence flags), so
    every number in the table is the pipeline's
    own accounting, not a model.  Dense-grid prediction skips the
    tile-dedup guard (prediction only); candidate runs enforce it, and a
    candidate the guard refuses (tiles inseparable at that length) is
    reported as failed rather than aborting the sweep.

    ``--scan-draws N`` repeats the whole sweep on N independent seeded sample
    draws (seeds ``--scan-seed``, +1, +2, ...) and reports each draw's
    recommendation and their agreement: a recommendation that flips between
    draws is a property of the sample subset, not of the panel, and a user
    without external truth has no other way to see that.
    """

    def call(argv: List[str]) -> int:
        try:
            return main(argv)
        except SystemExit as exc:  # parser.exit from a failed sub-run
            return int(exc.code or 0)

    probe_design = Path(args.probe_design)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    lengths = sorted({int(x) for x in args.lengths.split(",") if x.strip()})
    if not lengths:
        raise MadcError("--lengths must name at least one candidate tag length")
    fastqs = discover_fastqs(Path(args.fastq_dir), args.fastq_glob)
    if args.scan_draws == 1:
        subset = select_scan_subset(fastqs, args.scan_samples, args.scan_seed)
        _scan_once(args, call, probe_design, lengths, subset, out_dir)
        return 0
    base_seed = 0 if args.scan_seed is None else args.scan_seed
    draws = []
    for draw in range(args.scan_draws):
        seed = base_seed + draw
        subset = select_scan_subset(fastqs, args.scan_samples, seed)
        print(f"\nwindow-scan: draw {draw + 1}/{args.scan_draws} (seed {seed})", file=sys.stderr)
        report = _scan_once(args, call, probe_design, lengths, subset, out_dir / f"draw_{draw + 1}")
        draws.append({"draw": draw + 1, "seed": seed, "fastqs_scanned": report["inputs"]["fastqs_scanned"],
                      "recommended_tag_length": report["recommended_tag_length"],
                      "information_saturation_length": report["information_saturation_length"],
                      "eligible_lengths": report["eligible_lengths"], "candidates": report["candidates"]})
    votes = Counter(d["recommended_tag_length"] for d in draws)
    top = max(votes.values())
    consensus = min(length for length, n in votes.items() if n == top)
    stability = votes[consensus] / len(draws)
    by_length = {}
    for L in lengths:
        rows = [c for d in draws for c in d["candidates"] if c["tag_length"] == L]
        if rows:
            by_length[str(L)] = {
                "draws_completed": len(rows),
                "draws_eligible": sum(1 for r in rows if r["eligible"]),
                "draws_recommended": votes.get(L, 0),
                "survival_mean": sum(r["survival"] for r in rows) / len(rows),
                "survival_min": min(r["survival"] for r in rows),
                "survival_max": max(r["survival"] for r in rows),
                "novel_rows_per_100_markers_mean": sum(r["novel_rows_per_100_markers"] for r in rows) / len(rows),
            }
    summary = {
        "pipeline": "MADCATS window-scan",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parameters": {"lengths": lengths, "scan_samples": args.scan_samples, "scan_draws": args.scan_draws,
                       "base_seed": base_seed, "info_band": args.info_band},
        "inputs": {"probe_design": str(probe_design.resolve()), "fastqs_available": len(fastqs)},
        "draws": draws,
        "by_length": by_length,
        "recommended_tag_length": consensus,
        "recommendation_votes": {str(L): n for L, n in sorted(votes.items())},
        "recommendation_stability": stability,
    }
    _write_json(out_dir / "window_scan.qc.json", summary)
    with (out_dir / "window_scan.tsv").open("w") as out:
        out.write("tag_length\tdraws_eligible\tdraws_recommended\tsurvival_mean\tsurvival_min\tsurvival_max\tnovel_rows_per_100_markers_mean\n")
        for L, row in by_length.items():
            out.write(f"{L}\t{row['draws_eligible']}\t{row['draws_recommended']}\t{row['survival_mean']:.4f}\t"
                      f"{row['survival_min']:.4f}\t{row['survival_max']:.4f}\t{row['novel_rows_per_100_markers_mean']:.2f}\n")
    print(
        f"\nwindow-scan: {args.scan_draws} draws of {args.scan_samples} samples; recommendations "
        + ", ".join(f"{L} x{n}" for L, n in sorted(votes.items()))
        + f"; consensus {consensus} ({stability:.0%} of draws)"
        + ("" if stability == 1.0 else " — the recommendation depends on the sample subset: scan more samples or arbitrate externally")
        + f"; QC: {out_dir / 'window_scan.qc.json'}",
        file=sys.stderr,
    )
    return 0


def _scan_once(
    args: argparse.Namespace, call, probe_design: Path, lengths: Sequence[int], subset: Sequence[Path], out_dir: Path
) -> Dict[str, object]:
    """One sweep over ``lengths`` on one sample subset; writes the TSV and QC into
    ``out_dir`` and returns the QC report."""
    out_dir.mkdir(parents=True, exist_ok=True)
    scan_in = out_dir / "scan_input"
    scan_in.mkdir(exist_ok=True)
    for fastq in subset:
        link = scan_in / fastq.name
        if not link.exists():
            os.symlink(fastq.resolve(), link)
    # One coverage pass serves every length on the dense prediction grid.
    _normalize_init(str(probe_design), args.kmer_size)
    expansions = _NORM_STATE["expansions"]
    members_of = _NORM_STATE["members_of"]
    group_of = _NORM_STATE["group_of"]
    cover: Counter = Counter()
    print(f"window-scan: coverage pass over {len(subset)} FASTQ(s)", file=sys.stderr)
    for fastq in subset:
        _, c = _normalize_coverage_scan(str(fastq))
        cover.update(c)
    member_cover: Counter = Counter()
    group_reads: Counter = Counter()
    for (group, offset, length), count in cover.items():
        group_reads[group] += count
        for member in members_of.get(group, ()):
            member_cover[(member, offset, length)] += count

    def predicted_retention(tag_length: int) -> Optional[float]:
        windows = choose_tag_windows(
            expansions, member_cover, tag_length, args.anchor_length, "auto",
            str_penalty=args.str_penalty,
        )
        num = den = 0.0
        for marker, info in windows.items():
            group = group_of.get(marker, marker)
            weight = group_reads.get(group, 0) / max(1, len(members_of.get(group, (marker,))))
            if info["coverage"] is None or not weight:
                continue
            num += info["coverage"] * weight
            den += weight
        return num / den if den else None

    grid_lo = min(21, lengths[0])
    grid = sorted({*range(grid_lo, lengths[-1] + 1, args.grid_step), *lengths})
    dense = {L: predicted_retention(L) for L in grid}

    rows = []
    for L in lengths:
        ldir = out_dir / f"len_{L}"
        norm_dir = ldir / "normalized"
        code = call([
            "normalize-reads", "--probe-design", str(probe_design),
            "--fastq-dir", str(scan_in), "--tag-length", str(L),
            "--anchor-length", str(args.anchor_length), "--kmer-size", str(args.kmer_size),
            "--str-penalty", str(args.str_penalty), "--group-routing", args.group_routing,
            "--jobs", str(args.jobs), "--out-dir", str(norm_dir),
        ])
        if code != 0:
            print(f"window-scan: length {L}: normalize-reads failed (exit {code}); skipped", file=sys.stderr)
            continue
        madc_path = ldir / "scan_MADC.csv"
        code = call([
            "generate", "--fastq-dir", str(norm_dir), "--fastq-glob", "*.fastq",
            "--reference-tags", str(norm_dir / "window_reference_tags.csv"),
            "--tag-length", str(L), "--jobs", str(args.jobs),
            "--cache-dir", str(ldir / "cache"), "--output", str(madc_path),
            *shlex.split(args.generate_args),
        ])
        if code != 0:
            print(f"window-scan: length {L}: generate failed (exit {code}); skipped", file=sys.stderr)
            continue
        with (norm_dir / "normalize_reads.qc.json").open() as handle:
            nqc = json.load(handle)
        with open(str(madc_path) + ".qc.json") as handle:
            gqc = json.load(handle)
        totals = nqc["totals"]
        grand = sum(v for k, v in totals.items() if k != "tile_extra_emissions")
        ra = gqc["read_accounting"]
        assigned = ra["reads_assigned_full_length"]
        churn = scan_churn(ra)
        novel = sum(gqc["output_summary"]["novel_rows_by_evidence_route"].values())
        rows.append({
            "tag_length": L,
            "retention_predicted": dense.get(L),
            "retention_measured": totals.get("reads_extracted", 0) / max(1, grand),
            "assignment_rate": assigned / max(1, ra["reads_total"]),
            "clean_rate": 1.0 - churn / max(1, assigned),
            "novel_rows_per_100_markers": 100.0 * novel / max(1, gqc["loci"]),
            "loci_with_windows": gqc["loci"],
        })
    if not rows:
        raise MadcError("window-scan: no candidate length completed")

    max_info = max(r["novel_rows_per_100_markers"] for r in rows)
    for r in rows:
        r["survival"] = r["retention_measured"] * r["assignment_rate"] * r["clean_rate"]
    best_survival = max(r["survival"] for r in rows)
    for r in rows:
        r["eligible"] = r["survival"] >= (1.0 - args.info_band) * best_survival
    eligible = [r for r in rows if r["eligible"]]
    recommended = min(
        (r for r in eligible
         if r["novel_rows_per_100_markers"] == max(e["novel_rows_per_100_markers"] for e in eligible)),
        key=lambda r: r["tag_length"],
    )["tag_length"]
    saturation = next(
        (r["tag_length"] for r in rows
         if max_info and r["novel_rows_per_100_markers"] >= 0.95 * max_info),
        rows[0]["tag_length"],
    )

    tsv_path = out_dir / "window_scan.tsv"
    with tsv_path.open("w") as out:
        out.write("tag_length\tretention_predicted\tretention_measured\tassignment_rate\t"
                  "clean_rate\tnovel_rows_per_100_markers\tsurvival\teligible\n")
        for r in rows:
            predicted = "" if r["retention_predicted"] is None else f"{r['retention_predicted']:.4f}"
            out.write(f"{r['tag_length']}\t{predicted}\t"
                      f"{r['retention_measured']:.4f}\t{r['assignment_rate']:.4f}\t"
                      f"{r['clean_rate']:.4f}\t{r['novel_rows_per_100_markers']:.2f}\t"
                      f"{r['survival']:.4f}\t{int(r['eligible'])}\n")
        out.write("\n# dense predicted-retention grid (coverage pass only)\n")
        out.write("grid_length\tretention_predicted\n")
        for L in grid:
            if dense[L] is not None:
                out.write(f"{L}\t{dense[L]:.4f}\n")
    report = {
        "pipeline": "MADCATS window-scan",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parameters": {
            "lengths": lengths, "scan_samples": len(subset), "grid_step": args.grid_step,
            "anchor_length": args.anchor_length, "kmer_size": args.kmer_size,
            "str_penalty": args.str_penalty, "group_routing": args.group_routing,
        },
        "inputs": {"probe_design": str(probe_design.resolve()), "fastqs_scanned": [f.name for f in subset]},
        "candidates": rows,
        "dense_retention_grid": {str(L): dense[L] for L in grid if dense[L] is not None},
        "info_band": args.info_band,
        "eligible_lengths": sorted(r["tag_length"] for r in eligible),
        "recommended_tag_length": recommended,
        "information_saturation_length": saturation,
    }
    _write_json(out_dir / "window_scan.qc.json", report)

    def bar(value: float) -> str:
        return "#" * int(round(value * 30))

    print(f"\nwindow-scan over {len(subset)} sample(s), {rows[0]['loci_with_windows']}+ markers:", file=sys.stderr)
    print(f"{'len':>5s} {'retain':>7s} {'assign':>7s} {'clean':>7s} {'surviv':>7s} {'info/100':>9s}", file=sys.stderr)
    for r in rows:
        mark = " <- recommended" if r["tag_length"] == recommended else ("" if r["eligible"] else "  (below survival band)")
        print(f"{r['tag_length']:>5d} {r['retention_measured']:>7.1%} {r['assignment_rate']:>7.1%} "
              f"{r['clean_rate']:>7.1%} {r['survival']:>7.3f} {r['novel_rows_per_100_markers']:>9.1f}{mark}", file=sys.stderr)
        print(f"      surviv |{bar(r['survival'])}", file=sys.stderr)
        if max_info:
            print(f"      info   |{bar(r['novel_rows_per_100_markers'] / max_info)}", file=sys.stderr)
    print(f"recommended tag length: {recommended} (survival-first, info tiebreak within "
          f"{args.info_band:.0%} band; information saturates at {saturation}); "
          f"table: {tsv_path}; QC: {out_dir / 'window_scan.qc.json'}", file=sys.stderr)
    return report


def _normalize_reused_windows(
    args: argparse.Namespace, expansions: Mapping[str, Tuple[str, str, int]]
) -> Tuple[Dict[str, Dict[str, object]], List[Locus], List[Dict[str, object]]]:
    """Frozen-window mode: replay a previous run's window panel verbatim
    instead of learning placement from this run's coverage.  Coverage-
    learned windows drift between runs (read geometry differs run to
    run), and a moved window makes every microhaplotype at that marker
    a new sequence — indistinguishable from novel biology downstream.
    Reuse is therefore validated hard: each panel row must sit at
    exactly one position in the design expansions (both tags matching
    pins the window to the variant, since Ref and Alt differ only
    there), and any mismatch aborts the run rather than silently
    excising a partial or shifted panel.  Returns the windows, the loci,
    and the design markers absent from the reused panel (skipped)."""
    reused = load_reference_tags(Path(args.reuse_tags), tag_length=args.tag_length)
    problems: List[str] = []
    windows: Dict[str, Dict[str, object]] = {}
    loci: List[Locus] = []
    for locus in sorted(reused, key=lambda l: l.marker):
        expansion = expansions.get(locus.marker)
        if expansion is None:
            problems.append(f"{locus.marker}: not a SNP target in --probe-design")
            continue
        ref_g, alt_g, snp = expansion
        # The variant sits at target position `snp` and at tag position
        # `variant_pos`, so the window start is their difference — no
        # substring search, hence no risk of a repeat-shifted false match.
        start = snp - locus.variant_pos
        end = start + len(locus.ref_tag)
        if start < 0 or ref_g[start:end] != locus.ref_tag or alt_g[start:end] != locus.alt_tag:
            problems.append(f"{locus.marker}: reused tags do not match the design target around the variant")
            continue
        windows[locus.marker] = {
            "start": start, "coverage": None, "max_run": max_tandem_run(locus.ref_tag),
        }
        loci.append(
            Locus(
                index=len(loci), marker=locus.marker,
                left_primer=locus.left_primer, right_primer=locus.right_primer,
                ref_base=locus.ref_base, alt_base=locus.alt_base,
                variant_pos=locus.variant_pos, ref_tag=locus.ref_tag, alt_tag=locus.alt_tag,
                locus_group=locus.locus_group,
            )
        )
    # Same distinct-window invariant fresh placement enforces: two group
    # members with identical Ref tags tie every ref-carrying read into
    # an ambiguity drop.  A frozen panel is replayed verbatim, so a twin
    # here is a panel defect to report, not to silently nudge.
    seen_content: Dict[Tuple[str, str], str] = {}
    for locus in loci:
        key = (locus.group, locus.ref_tag)
        other = seen_content.get(key)
        if other is not None and locus.locus_group:
            problems.append(
                f"LocusGroup {locus.group}: members {other} and {locus.marker} reuse an "
                f"identical window — ref-carrying reads would all tie and be dropped"
            )
        seen_content[key] = locus.marker
    if problems:
        raise MadcError(
            "--reuse-tags panel is inconsistent with --probe-design; refusing a partial reuse "
            "(subset the panel CSV deliberately if that is intended): " + "; ".join(problems)
        )
    skipped = [
        {"marker": marker, "reason": "not_in_reused_panel"} for marker in expansions if marker not in windows
    ]
    return windows, loci, skipped


def _normalize_coverage_pass(args: argparse.Namespace, fastqs: Sequence[Path], probe_design: Path) -> Counter:
    """An optional coverage pass records, per marker, how observed reads land
    on the target ((start offset, length) classes), so window choice can
    maximize the fraction of reads that will fully cover the excised window
    instead of guessing.  Coverage is observed per locus GROUP; every tile
    of a group sees the same fragment-position profile."""
    cover: Counter = Counter()
    if args.window_placement == "auto":
        print(f"normalize-reads: coverage pass over {len(fastqs)} FASTQ(s)", file=sys.stderr)
        if args.jobs == 1:
            for fastq in fastqs:
                _, c = _normalize_coverage_scan(str(fastq))
                cover.update(c)
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=args.jobs, initializer=_normalize_init,
                initargs=(str(probe_design), args.kmer_size),
            ) as executor:
                for _, c in executor.map(_normalize_coverage_scan, map(str, fastqs)):
                    cover.update(c)
    members_lookup = _NORM_STATE["members_of"]
    member_cover: Counter = Counter()
    for (group, offset, length), count in cover.items():
        for member in members_lookup.get(group, ()):
            member_cover[(member, offset, length)] += count
    return member_cover


def _normalize_learned_windows(
    args: argparse.Namespace, fastqs: Sequence[Path], probe_design: Path, expansions: Mapping[str, Tuple[str, str, int]]
) -> Tuple[Dict[str, Dict[str, object]], List[Locus]]:
    """Fresh placement: choose each marker's window from observed coverage,
    then materialize the chosen windows as an ordinary reference-tag panel:
    the excised window plays the role a primer-defined tag plays for
    DArTag, with a panel-unique prefix as the assignment "primer"
    downstream generation matches reads against.  Anchor lengths are
    per locus: every marker keeps the --anchor-length minimum, and only
    colliding prefix groups grow (choose_anchor_lengths, the same
    allocator infer-design uses) — at panel scale prefix collisions are
    a certainty (a 17.5K genome-wide panel measured 111 colliding
    12-mers), so a single global length would either fail or tax every
    marker's window placement for a few markers' collisions."""
    member_cover = _normalize_coverage_pass(args, fastqs, probe_design)
    windows = choose_tag_windows(
        expansions, member_cover, args.tag_length, args.anchor_length, args.window_placement,
        str_penalty=args.str_penalty, groups=_NORM_STATE["group_of"],
    )
    window_tags = {}
    for marker in sorted(windows):
        ref_g, alt_g, snp = expansions[marker]
        start = windows[marker]["start"]
        window_tags[marker] = (
            ref_g[start : start + args.tag_length],
            alt_g[start : start + args.tag_length],
            snp - start,
        )
    anchor_lengths = choose_anchor_lengths(
        {m: (t[0], t[2]) for m, t in window_tags.items()},
        args.anchor_length,
        args.tag_length,
    )
    loci: List[Locus] = []
    group_names = _NORM_STATE["group_of"]
    for marker in sorted(window_tags):
        ref_tag, alt_tag, variant_offset = window_tags[marker]
        group = group_names.get(marker, marker)
        loci.append(
            Locus(
                index=len(loci), marker=marker,
                left_primer=ref_tag[: anchor_lengths[marker]], right_primer="",
                ref_base=ref_tag[variant_offset], alt_base=alt_tag[variant_offset],
                variant_pos=variant_offset, ref_tag=ref_tag, alt_tag=alt_tag,
                locus_group=group if group != marker else "",
            )
        )
    return windows, loci


def _predicted_full_coverage(member_cover: Mapping[Tuple[str, int, int], int], marker: str, start: int, tag_length: int) -> Optional[float]:
    """Fraction of a marker's located reads that fully cover the window at ``start``
    (the quantity auto placement maximises), or ``None`` without observations."""
    observations = [(offset, length, count) for (m, offset, length), count in member_cover.items() if m == marker]
    total = sum(count for _, _, count in observations)
    if not total:
        return None
    covered = sum(count for offset, length, count in observations if offset <= start and offset + length >= start + tag_length)
    return covered / total


def _normalize_freeze_health(
    args: argparse.Namespace,
    fastqs: Sequence[Path],
    probe_design: Path,
    expansions: Mapping[str, Tuple[str, str, int]],
    windows: Dict[str, Dict[str, object]],
    loci: Sequence[Locus],
) -> Optional[Dict[str, object]]:
    """How well a frozen panel fits this run's reads.  Freezing is what keeps
    haplotype identity comparable across runs, but read geometry drifts run
    to run, so the coverage pass that fresh placement would use is run anyway
    and every frozen window is scored against it: its predicted full coverage,
    and the coverage the window fresh placement would choose today.  The
    user sees what freezing costs this run and when a panel refresh is due.
    ``--window-placement center`` skips the pass (no read data needed).  The
    fresh comparison keeps the frozen panel's own anchor geometry (its
    shortest anchor is the --anchor-length it was built with; only colliding
    groups grow), so the two placements answer the same question."""
    if args.window_placement != "auto" or not loci:
        return None
    member_cover = _normalize_coverage_pass(args, fastqs, probe_design)
    frozen_expansions = {marker: expansions[marker] for marker in windows}
    anchor_length = min(len(locus.left_primer) for locus in loci)
    fresh = choose_tag_windows(
        frozen_expansions, member_cover, args.tag_length, anchor_length, "auto",
        str_penalty=args.str_penalty, groups=_NORM_STATE["group_of"],
    )
    weighted_frozen = weighted_fresh = weight_total = 0.0
    losses = []
    for marker, info in windows.items():
        coverage = _predicted_full_coverage(member_cover, marker, info["start"], args.tag_length)
        info["coverage"] = coverage
        choice = fresh.get(marker)
        info["fresh_start"] = choice["start"] if choice else None
        info["fresh_coverage"] = choice["coverage"] if choice else None
        if coverage is None or choice is None or choice["coverage"] is None:
            continue
        reads = sum(count for (m, _, _), count in member_cover.items() if m == marker)
        weighted_frozen += coverage * reads
        weighted_fresh += choice["coverage"] * reads
        weight_total += reads
        losses.append((marker, choice["coverage"] - coverage, info["start"], choice["start"]))
    worst = sorted(losses, key=lambda item: -item[1])[:10]
    health = {
        "markers_scored": len(losses),
        "read_weighted_full_coverage_frozen": weighted_frozen / weight_total if weight_total else None,
        "read_weighted_full_coverage_fresh": weighted_fresh / weight_total if weight_total else None,
        "markers_losing_over_0.10": sum(1 for _, loss, _, _ in losses if loss > 0.10),
        "markers_moved_by_fresh_placement": sum(1 for _, _, frozen_start, fresh_start in losses if frozen_start != fresh_start),
        "worst_markers": [
            {"marker": marker, "coverage_loss": loss, "frozen_start": frozen_start, "fresh_start": fresh_start}
            for marker, loss, frozen_start, fresh_start in worst if loss > 0
        ],
    }
    if weight_total:
        print(
            f"normalize-reads: frozen windows fully cover {health['read_weighted_full_coverage_frozen']:.1%} of located "
            f"reads (fresh placement would cover {health['read_weighted_full_coverage_fresh']:.1%}); "
            f"{health['markers_losing_over_0.10']} marker(s) lose more than 10 points — see freeze_health in the QC",
            file=sys.stderr,
        )
    return health


def _normalize_shared_anchor_groups(loci: Sequence[Locus]) -> List[List[str]]:
    """Markers whose windows are identical all the way to their growth cap
    (the variant) keep a shared anchor: no prefix length can separate
    them, and generation already resolves shared-anchor reads by full-tag
    identity, dropping true ties as ambiguous.  Surface them rather than
    refuse — they are the same paralog signal the k-mer index reports."""
    anchor_groups: Dict[str, List[str]] = defaultdict(list)
    for locus in loci:
        anchor_groups[locus.left_primer].append(locus.marker)
    shared_anchor_groups = sorted(ms for ms in anchor_groups.values() if len(ms) > 1)
    if shared_anchor_groups:
        preview = "; ".join(", ".join(group) for group in shared_anchor_groups[:3])
        print(
            f"normalize-reads: {len(shared_anchor_groups)} marker group(s) share an anchor "
            f"(windows identical to the variant; resolved by tag identity downstream): {preview}",
            file=sys.stderr,
        )
    return shared_anchor_groups


def _normalize_extract_all(
    args: argparse.Namespace,
    fastqs: Sequence[Path],
    out_dir: Path,
    probe_design: Path,
    windows: Mapping[str, Mapping[str, object]],
) -> Tuple[Counter, Dict[str, Dict[str, int]]]:
    """Excise every read's window into the per-sample FASTQs; totals and per-sample stats."""
    print(f"normalize-reads: extracting {len(fastqs)} FASTQ(s)", file=sys.stderr)
    window_starts = {marker: windows[marker]["start"] for marker in windows}
    totals: Counter = Counter()
    per_sample = {}
    members_of = _NORM_STATE["members_of"]
    tasks = [
        (str(fastq), str(out_dir), args.tag_length, window_starts, members_of,
         args.group_routing)
        for fastq in fastqs
    ]
    if args.jobs == 1:
        results = [_normalize_extract(task) for task in tasks]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.jobs, initializer=_normalize_init,
            initargs=(str(probe_design), args.kmer_size),
        ) as executor:
            results = list(executor.map(_normalize_extract, tasks))
    for sample, stats in results:
        per_sample[sample] = dict(stats)
        totals.update(stats)
    return totals, per_sample


def _normalize_report(
    args: argparse.Namespace,
    probe_design: Path,
    fastqs: Sequence[Path],
    tags_path: Path,
    expansions: Mapping[str, Tuple[str, str, int]],
    windows: Mapping[str, Mapping[str, object]],
    loci: Sequence[Locus],
    skipped: Sequence[Dict[str, object]],
    shared_anchor_groups: Sequence[Sequence[str]],
    totals: Counter,
    per_sample: Mapping[str, Dict[str, int]],
) -> Dict[str, object]:
    return {
        "pipeline": "MADCATS normalize-reads",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parameters": {
            "tag_length": args.tag_length, "anchor_length": args.anchor_length,
            "kmer_size": args.kmer_size,
            "window_placement": "reuse" if args.reuse_tags else args.window_placement,
            "reuse_tags": str(Path(args.reuse_tags).resolve()) if args.reuse_tags else None,
            "str_penalty": args.str_penalty,
            "group_routing": args.group_routing,
        },
        "inputs": {"probe_design": str(probe_design.resolve()), "fastqs": len(fastqs)},
        "reference_tags": str(tags_path.resolve()),
        "totals": dict(totals),
        "markers": {
            marker: {
                "window_start": windows[marker]["start"],
                "variant_position": expansions[marker][2] - windows[marker]["start"] + 1,
                "predicted_full_coverage": windows[marker]["coverage"],
                "max_tandem_run_bp": windows[marker]["max_run"],
                **({"fresh_window_start": windows[marker]["fresh_start"], "fresh_predicted_full_coverage": windows[marker]["fresh_coverage"]}
                   if "fresh_start" in windows[marker] else {}),
            }
            for marker in sorted(windows)
        },
        "skipped_markers": list(skipped),
        "locus_groups": {
            "groups": len(group_sizes := Counter(l.group for l in loci)),
            "tiled_groups": sum(1 for size in group_sizes.values() if size > 1),
            "window_dedup_nudged": sorted(
                marker for marker, info in windows.items() if "dedup_nudged_from" in info
            ),
        },
        "anchor_summary": {
            "length_histogram": dict(sorted(Counter(len(l.left_primer) for l in loci).items())),
            "shared_anchor_groups": [list(g) for g in shared_anchor_groups],
        },
        "per_sample": dict(per_sample),
    }


def run_normalize_reads(args: argparse.Namespace) -> int:
    """Normalize variable-start, bidirectional amplicon reads into pseudo-tags.

    Platforms such as AgriSeq sequence amplicons from both ends with primer
    positions outside the submitted target context, so reads are neither
    oriented nor primer-started.  Each read is located by unique k-mer
    against the probe-design targets, oriented forward, and a fixed tag-
    length window containing the variant is excised with its qualities.  The
    command emits pipeline-ready per-sample FASTQs, the matching window
    reference-tags CSV, and a QC ledger; generation then runs unchanged
    (``--primer-mismatches 1``, the default since 0.9.0, matters most on
    indel-prone chemistries).
    """
    probe_design = Path(args.probe_design)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fastqs = discover_fastqs(Path(args.fastq_dir), args.fastq_glob)
    _normalize_init(str(probe_design), args.kmer_size)
    expansions = dict(_NORM_STATE["expansions"])
    targets = read_probe_targets(probe_design)
    skipped = [
        {"marker": marker, "reason": "non_snp_or_invalid_bracket"}
        for marker in targets if marker not in expansions
    ]
    # Windows come from a frozen panel (replayed verbatim) or are learned
    # from this run's coverage; either way they become a reference-tag panel.
    freeze_health = None
    if args.reuse_tags:
        windows, loci, not_reused = _normalize_reused_windows(args, expansions)
        skipped.extend(not_reused)
        freeze_health = _normalize_freeze_health(args, fastqs, probe_design, expansions, windows, loci)
    else:
        windows, loci = _normalize_learned_windows(args, fastqs, probe_design, expansions)
    shared_anchor_groups = _normalize_shared_anchor_groups(loci)
    tags_path = Path(args.output_tags) if args.output_tags else out_dir / "window_reference_tags.csv"
    write_reference_tags(tags_path, loci)
    if loci:
        load_reference_tags(tags_path, tag_length=args.tag_length)  # round-trip self-check

    totals, per_sample = _normalize_extract_all(args, fastqs, out_dir, probe_design, windows)
    report = _normalize_report(
        args, probe_design, fastqs, tags_path, expansions, windows, loci, skipped,
        shared_anchor_groups, totals, per_sample,
    )
    if freeze_health is not None:
        report["freeze_health"] = freeze_health
    report_path = Path(args.report) if args.report else out_dir / "normalize_reads.qc.json"
    _write_json(report_path, report)
    # tile_extra_emissions counts additional per-window copies of already-
    # extracted fragments (covered routing), not input reads — keep it out
    # of the read-fate denominator.
    grand = sum(count for key, count in totals.items() if key != "tile_extra_emissions")
    extra = totals.get("tile_extra_emissions", 0)
    print(
        f"Normalized {totals['reads_extracted']:,}/{grand:,} reads "
        f"({totals['reads_extracted']/grand:.1%}) across {len(loci)} markers to {out_dir}"
        + (f" (+{extra:,} extra tile emissions)" if extra else "")
        + f"; panel: {tags_path}; QC: {report_path}",
        file=sys.stderr,
    )
    return 0


def _resolve_generate_defaults(args: argparse.Namespace) -> None:
    """Defaults that depend on other options, resolved once the parser is done."""
    if args.min_partial_length is None:
        args.min_partial_length = max(1, args.tag_length - PARTIAL_TAIL_DEFAULT)
    if args.partial_read_mode != "off" and args.min_partial_length >= args.tag_length:
        raise MadcError("--min-partial-length must be shorter than --tag-length when partial reads are enabled")


def _define_panel(args: argparse.Namespace) -> List[Locus]:
    """Stage 1: panel definition — the canonical tags every read is measured
    against, from either input route."""
    if args.reference_tags:
        return load_reference_tags(Path(args.reference_tags), tag_length=args.tag_length)
    if not args.primers or not args.probe_design:
        raise MadcError("Provide --reference-tags, or provide both --primers and --probe-design")
    return build_panel(Path(args.primers), Path(args.probe_design), args.tag_length)


def _evidence_rule(args: argparse.Namespace) -> EvidenceRule:
    """One evidence rule is shared by protection, discovery, and emission
    (stages 4-6), so no stage can disagree about which sequences are real."""
    return EvidenceRule(
        min_global_count=args.min_global_count,
        discovery_min_quality=args.discovery_min_quality,
        quality_min_clean_reads=args.quality_min_clean_reads,
        reproducible_min_reads=args.reproducible_min_reads,
        reproducible_min_samples=args.reproducible_min_samples,
        depth_min_reads=args.depth_min_reads,
        min_supporting_freq=args.min_supporting_freq,
    )


def _assignment_settings(args: argparse.Namespace) -> AssignmentSettings:
    """The read-pass settings, which are also the cache key."""
    return AssignmentSettings(
        tag_length=args.tag_length,
        max_offtarget_mismatches=math.floor((1.0 - args.min_identity) * args.tag_length + 1e-9),
        primer_mismatches=args.primer_mismatches,
        partial_min_length=args.min_partial_length,
        gap_band=args.max_indel_length if args.gapped_rescue == "on" else 0,
        leading_anchor=(getattr(args, "leading_anchor", None) or "").upper(),
    )


def _discover_samples(args: argparse.Namespace) -> Tuple[List[Sample], Path, List[Path]]:
    """The run's samples: one per FASTQ, named from the target manifest or the lookup."""
    fastq_directory = Path(args.fastq_dir)
    fastqs = discover_fastqs(fastq_directory, args.fastq_glob)
    if args.target_metadata:
        metadata_paths = [Path(path) for path in args.target_metadata]
    else:
        metadata_paths = sorted(fastq_directory.glob("targets_*.csv"))
    metadata = load_target_metadata(metadata_paths)
    samples = make_samples(fastqs, metadata, args.sample_name_field)
    if args.sample_lookup:
        # The lookup replaces output column names wholesale (validated as
        # complete and collision-free), while sample_id stays FASTQ-derived —
        # it keys the cache and fills the preamble's targetid row.
        lookup = load_sample_lookup(Path(args.sample_lookup), [s.sample_id for s in samples])
        samples = [dataclasses.replace(s, display_name=lookup[s.sample_id]) for s in samples]
    return samples, fastq_directory, metadata_paths


def _warn_about_evidence_routes(rule: EvidenceRule, sample_count: int) -> None:
    if rule.reproducibility_route_enabled and sample_count < rule.reproducible_min_samples:
        remaining = " and no other evidence route is enabled" if not (
            rule.depth_route_enabled or rule.quality_route_enabled
        ) else ""
        print(
            f"warning: --reproducible-min-samples {rule.reproducible_min_samples} exceeds the "
            f"{sample_count} available sample(s); the reproducibility route can never pass{remaining}",
            file=sys.stderr,
        )
    if not rule.any_route_enabled:
        print(
            "note: every evidence route is disabled; all sequences at or above "
            f"--min-global-count {rule.min_global_count} will be emitted (exhaustive mode)",
            file=sys.stderr,
        )


def _report_sample(index: int, total: int, result: SampleResult) -> None:
    print(
        f"[{index}/{total}] {result.sample.sample_id}: "
        f"{result.qc.get('reads_assigned', 0):,} assigned"
        + (" (cache)" if result.cache_hit else ""),
        file=sys.stderr,
    )


def _assign_samples(
    samples: Sequence[Sample],
    loci: Sequence[Locus],
    settings: AssignmentSettings,
    cache_dir: Optional[str],
    jobs: int,
) -> List[SampleResult]:
    """Stage 2: per-sample read processing — the only FASTQ-touching, expensive
    stage — parallel across samples and skipped entirely on cache hits."""
    primer_index = PrimerIndex(loci)
    cache_fingerprint = _cache_fingerprint(loci, settings.cache_arguments())[:16]
    worker_arguments = []
    for sample in samples:
        cache_path = Path(cache_dir) / f"{sample.sample_id}.{cache_fingerprint}.json.gz" if cache_dir else None
        worker_arguments.append((sample, loci, primer_index, settings, cache_path))
    print(f"Processing {len(samples)} samples across {len(loci)} loci with {jobs} worker(s)...", file=sys.stderr)
    if jobs == 1:
        results = []
        for index, arguments in enumerate(worker_arguments, start=1):
            result = _process_sample_worker(arguments)
            results.append(result)
            _report_sample(index, len(samples), result)
        return results
    by_id = {}
    with concurrent.futures.ProcessPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(_process_sample_worker, arguments): arguments[0].sample_id for arguments in worker_arguments}
        for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            by_id[result.sample.sample_id] = result
            _report_sample(index, len(samples), result)
    return [by_id[sample.sample_id] for sample in samples]


def _generate_provenance(args: argparse.Namespace) -> str:
    """The whole provenance string lives in the file's first cell: enough to
    rebuild the run (pipeline, version, every generation parameter) without
    the QC sidecar.  No timestamp — the table must stay byte-reproducible for
    regression gating."""
    return " ".join(
        [
            "MADCATS",
            f"v{__version__}",
            f"tag_length={args.tag_length}",
            f"min_identity={args.min_identity}",
            f"primer_mismatches={args.primer_mismatches}",
            f"gapped_rescue={args.gapped_rescue}",
            f"max_indel_length={args.max_indel_length}",
            f"annotations={args.annotations}",
            f"partial_read_mode={args.partial_read_mode}",
            f"min_partial_length={args.min_partial_length}",
            *([f"leading_anchor={args.leading_anchor.upper()}"] if getattr(args, "leading_anchor", None) else []),
            f"denoise_distance={args.denoise_distance}",
            f"global_denoise_distance={args.global_denoise_distance}",
            f"denoise_alpha={args.denoise_alpha}",
            f"min_global_count={args.min_global_count}",
            f"discovery_min_quality={args.discovery_min_quality}",
            f"quality_min_clean_reads={args.quality_min_clean_reads}",
            f"reproducible_min_reads={args.reproducible_min_reads}",
            f"reproducible_min_samples={args.reproducible_min_samples}",
            f"depth_min_reads={args.depth_min_reads}",
            f"min_supporting_freq={args.min_supporting_freq}",
        ]
    )


def _output_options(args: argparse.Namespace, samples: Sequence[Sample], output: Path) -> OutputOptions:
    preamble = None
    if args.dart_preamble:
        preamble = build_dart_preamble(
            samples, args.preamble_project or output.name.split(".")[0], _generate_provenance(args)
        )
    return OutputOptions(
        number_match_ids=args.number_match_ids,
        preamble=preamble,
        annotations=args.annotations == "on",
        variants_path=output.with_suffix(output.suffix + ".variants.csv"),
    )


def _read_accounting(
    results: Sequence[SampleResult], partial_summary: Mapping[str, int], output_summary: Mapping[str, object]
) -> Dict[str, int]:
    """The ledger from raw records through emitted rows."""
    total = lambda key: sum(result.qc.get(key, 0) for result in results)
    return {
        "reads_total": total("reads_total"),
        "reads_no_leading_anchor": total("reads_no_leading_anchor"),
        "reads_assigned_full_length": total("reads_assigned"),
        "reads_assigned_gapped": total("reads_assigned_gapped"),
        "reads_n_mapped": total("reads_n_mapped"),
        "reads_n_dropped": total("reads_n_dropped"),
        "partial_reads_added_to_counts": partial_summary.get("reads_assigned", 0),
        "reads_in_final_count_table": sum(sum(result.counts.values()) for result in results),
        "reads_moved_by_sample_denoise": total("reads_absorbed_sample_denoise"),
        "reads_moved_by_global_denoise": total("reads_absorbed_global_denoise"),
        "reads_in_emitted_rows": output_summary["reads_in_emitted_rows"],
        "reads_dropped_below_min_global_count": output_summary["excluded_novel"]["below_min_global_count"]["reads"],
        "reads_dropped_failed_evidence": output_summary["excluded_novel"]["failed_evidence"]["reads"],
    }


def _generate_parameters(args: argparse.Namespace, settings: AssignmentSettings, output: Path) -> Dict[str, object]:
    return {
        "tag_length": args.tag_length,
        "minimum_identity": args.min_identity,
        "maximum_offtarget_mismatches": settings.max_offtarget_mismatches,
        "primer_mismatches": args.primer_mismatches,
        "gapped_rescue": args.gapped_rescue,
        "max_indel_length": args.max_indel_length,
        "annotations": args.annotations,
        "partial_read_mode": args.partial_read_mode,
        "minimum_partial_length": args.min_partial_length,
        "leading_anchor": (getattr(args, "leading_anchor", None) or "").upper() or None,
        "leading_anchor_detected": bool(getattr(args, "leading_anchor_detected", False)),
        "sample_denoise_distance": args.denoise_distance,
        "global_denoise_distance": args.global_denoise_distance,
        "denoise_alpha": args.denoise_alpha,
        "minimum_global_count": args.min_global_count,
        "discovery_min_quality": args.discovery_min_quality,
        "quality_min_clean_reads": args.quality_min_clean_reads,
        "reproducible_min_reads": args.reproducible_min_reads,
        "reproducible_min_samples": args.reproducible_min_samples,
        "depth_min_reads": args.depth_min_reads,
        "min_supporting_freq": args.min_supporting_freq,
        "dart_preamble": bool(args.dart_preamble),
        "preamble_project": (args.preamble_project or output.name.split(".")[0]) if args.dart_preamble else None,
        "sample_lookup": str(Path(args.sample_lookup).resolve()) if args.sample_lookup else None,
    }


def _write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def run_generate(args: argparse.Namespace) -> int:
    """The pipeline, stage by stage (see the module docstring)."""
    _resolve_generate_defaults(args)
    loci = _define_panel(args)
    rule = _evidence_rule(args)
    samples, fastq_directory, metadata_paths = _discover_samples(args)
    resolve_leading_anchor(
        args, [Path(sample.fastq) for sample in samples],
        {tag[:LEADING_ANCHOR_PREFIX] for locus in loci for tag in (locus.ref_tag, locus.alt_tag)}, "generate",
    )
    settings = _assignment_settings(args)
    _warn_about_evidence_routes(rule, len(samples))
    results = _assign_samples(samples, loci, settings, args.cache_dir, args.jobs)
    # Stage 3: optional partial-read salvage; stages 4-5: evidence and denoising
    # with evidence-passing sequences protected from absorption.
    results, partial_summary = apply_partial_reads(results, loci, mode=args.partial_read_mode, rule=rule)
    denoise = DenoiseSettings(args.denoise_distance, args.global_denoise_distance, args.denoise_alpha)
    results = finalize_counts(results, loci, denoise, rule, gap_band=settings.gap_band, jobs=args.jobs)
    # Stage 6: evidence-gated emission with forced canonical rows.
    output = Path(args.output)
    output_summary = write_madc(output, loci, results, rule, _output_options(args, samples, output), settings.gap_band)
    per_sample_emitted = output_summary.pop("per_sample_reads_in_emitted_rows")
    for result, emitted in zip(results, per_sample_emitted):
        result.qc["reads_in_emitted_rows"] = emitted
    # Stage 7: QC/provenance sidecar — inputs, parameters, panel fingerprints,
    # per-sample QC, and the read-accounting ledger.
    qc_path = Path(args.qc_output) if args.qc_output else output.with_suffix(output.suffix + ".qc.json")
    qc = {
        "pipeline": "MADCATS",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "output": str(output.resolve()),
        "inputs": {
            "fastq_directory": str(fastq_directory.resolve()),
            "reference_tags": str(Path(args.reference_tags).resolve()) if args.reference_tags else None,
            "primers": str(Path(args.primers).resolve()) if args.primers else None,
            "probe_design": str(Path(args.probe_design).resolve()) if args.probe_design else None,
            "target_metadata": [str(path.resolve()) for path in metadata_paths],
        },
        "panel_fingerprint_sha256": panel_fingerprint(loci),
        "parameters": _generate_parameters(args, settings, output),
        "loci": len(loci),
        "samples": len(samples),
        "read_accounting": _read_accounting(results, partial_summary, output_summary),
        "output_summary": output_summary,
        "partial_read_summary": partial_summary,
        "sample_qc": [
            {"sample_id": result.sample.sample_id, "fastq": result.sample.fastq, **result.qc}
            for result in results
        ],
    }
    _write_json(qc_path, qc)
    print(f"Wrote {output_summary['rows']:,} allele rows to {output}", file=sys.stderr)
    print(f"Wrote QC/provenance to {qc_path}", file=sys.stderr)
    return 0


def run_validate(args: argparse.Namespace) -> int:
    report = validate_madc(Path(args.generated), Path(args.expected))
    if args.report:
        report_path = Path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with report_path.open("w") as handle:
            json.dump(report, handle, indent=2, sort_keys=True)
            handle.write("\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0


def run_design(args: argparse.Namespace) -> int:
    """Reconstruct the panel from vendor design files and save it as a CSV.

    A thin wrapper around ``build_panel``: the value of the subcommand is
    that the reconstructed tags become a reviewable, platform-neutral
    artifact (the same CSV format every other route produces) instead of an
    in-memory step of generation.
    """
    loci = build_panel(Path(args.primers), Path(args.probe_design), args.tag_length)
    output = Path(args.output)
    write_reference_tags(output, loci)
    print(f"Wrote {len(loci):,} reference-tag definitions to {output}", file=sys.stderr)
    return 0


def _infer_candidate_windows(
    targets: Mapping[str, str], tag_length: int
) -> Tuple[Dict[str, Dict[Tuple[str, int], Tuple[str, str, int]]], List[Dict[str, object]], set]:
    """Stage 1: enumerate each marker's closed candidate set — every
    tag-length window of the expanded target, on both strands, that spans
    the variant.  Inference never invents sequence; it only chooses among
    these design-derived candidates, which is what keeps a wrong answer
    detectable (it must still be a real window of the submitted target).

    Returns the candidates per marker, the markers excluded with a reason,
    and the markers flagged as paralogs: sequence-identical candidates shared
    between markers are the paralog signature — reads there vote for several
    loci at once, so both sides carry a review flag whatever the outcome."""
    windows_by_marker: Dict[str, Dict[Tuple[str, int], Tuple[str, str, int]]] = {}
    excluded: List[Dict[str, object]] = []
    for marker, target in targets.items():
        match = BRACKET_RE.search(target)
        if not match:
            excluded.append({"marker": marker, "reason": "no_bracketed_variant"})
            continue
        first = "" if match.group(1) == "-" else match.group(1)
        second = "" if match.group(2) == "-" else match.group(2)
        if first == second:
            excluded.append({"marker": marker, "reason": "identical_alleles"})
            continue
        windows = enumerate_variant_windows(target, tag_length)
        if not windows:
            excluded.append({"marker": marker, "reason": "flanks_too_short_for_tag_length"})
            continue
        windows_by_marker[marker] = windows
    owners: Dict[str, set] = defaultdict(set)
    for marker, windows in windows_by_marker.items():
        for ref_window, alt_window, _offset in windows.values():
            owners[ref_window].add(marker)
            owners[alt_window].add(marker)
    paralog_flagged = set()
    for markers in owners.values():
        if len(markers) > 1:
            paralog_flagged.update(markers)
    return windows_by_marker, excluded, paralog_flagged


def _run_fastq_scans(fastqs: Sequence[Path], jobs: int, scan, initializer, initargs) -> Tuple[Counter, Dict[str, Counter]]:
    """Run one vote-collecting scan over every FASTQ, serially or in a pool whose
    workers rebuild the lookup tables from the initializer arguments."""
    aggregated: Counter = Counter()
    per_fastq: Dict[str, Counter] = {}
    if jobs == 1:
        initializer(*initargs)
        for index, fastq in enumerate(fastqs, start=1):
            name, votes = scan(str(fastq))
            per_fastq[name] = votes
            aggregated.update(votes)
            print(f"[{index}/{len(fastqs)}] scanned {Path(name).name}", file=sys.stderr)
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=jobs, initializer=initializer, initargs=initargs
        ) as executor:
            futures = [executor.submit(scan, str(fastq)) for fastq in fastqs]
            for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                name, votes = future.result()
                per_fastq[name] = votes
                aggregated.update(votes)
                print(f"[{index}/{len(fastqs)}] scanned {Path(name).name}", file=sys.stderr)
    return aggregated, per_fastq


def _infer_resolve_exact(
    windows_by_marker: Mapping[str, Mapping[Tuple[str, int], Tuple[str, str, int]]],
    exact_votes: Mapping[Tuple[str, Tuple[str, int]], int],
    dominance_ratio: float,
    paralog_flagged: set,
    excluded: List[Dict[str, object]],
) -> Tuple[Dict[str, Dict[str, object]], List[str]]:
    """Stage 2: decide each marker from exact read-prefix votes.  Exact
    evidence is authoritative — a read literally beginning with a candidate
    window is the assay demonstrating its register — but the winner must
    dominate the runner-up, because primer-started reads all agree on one
    window and a close race is the signature of paralogs or fragmentation.
    Returns the resolved loci and the markers with no exact support."""
    resolved: Dict[str, Dict[str, object]] = {}
    unresolved: List[str] = []
    for marker, windows in windows_by_marker.items():
        support = {key: exact_votes.get((marker, key), 0) for key in windows}
        ranked = sorted(support.items(), key=lambda item: (-item[1], item[0]))
        best_key, best = ranked[0]
        runner = ranked[1][1] if len(ranked) > 1 else 0
        if best == 0:
            unresolved.append(marker)
            continue
        if runner > 0 and best < dominance_ratio * runner:
            excluded.append(
                {"marker": marker, "reason": "ambiguous_windows", "best": best, "runner_up": runner}
            )
            continue
        # Fragmented reads give most runner-up windows a few stray votes, so
        # only genuinely close races earn a review flag.
        flags = ["close_runner_up"] if runner > 0 and best < 2 * dominance_ratio * runner else []
        if marker in paralog_flagged:
            flags.append("paralog_shared_candidates")
        resolved[marker] = {
            "route": "exact", "window": best_key, "support": best,
            "runner_up": runner, "supporting_samples": None, "flags": flags,
        }
    return resolved, unresolved


def _infer_resolve_fuzzy(
    args: argparse.Namespace,
    probe_design: Path,
    fastqs: Sequence[Path],
    max_mismatches: int,
    windows_by_marker: Mapping[str, Mapping[Tuple[str, int], Tuple[str, str, int]]],
    resolved: Dict[str, Dict[str, object]],
    unresolved: Sequence[str],
    paralog_flagged: set,
    excluded: List[Dict[str, object]],
) -> None:
    """Stage 3: a fuzzy tier for loci with zero exact support — germplasm can
    be diverged enough that no read matches any candidate verbatim.  Reads
    vote within the identity envelope instead, with already-resolved tags
    included as competition sinks, and every fuzzy anchor carries a review
    flag: this tier can be wrong in ways the exact tier cannot.  Adds to
    ``resolved`` and ``excluded`` in place."""
    if unresolved and args.fuzzy_min_reads > 0:
        if len(fastqs) < args.fuzzy_min_samples:
            print(
                f"warning: --fuzzy-min-samples {args.fuzzy_min_samples} exceeds the "
                f"{len(fastqs)} available FASTQ(s); no fuzzy anchor can qualify",
                file=sys.stderr,
            )
        sink_tags = []
        for marker, record in resolved.items():
            ref_window, alt_window, _offset = windows_by_marker[marker][record["window"]]
            sink_tags.extend((ref_window, alt_window))
        print(
            f"infer-design: fuzzy pass for {len(unresolved)} unresolved marker(s) "
            f"within {max_mismatches} mismatches, {len(sink_tags):,} resolved sinks",
            file=sys.stderr,
        )
        fuzzy_votes, fuzzy_per_fastq = _run_fastq_scans(
            fastqs,
            args.jobs,
            _infer_fuzzy_scan,
            _infer_fuzzy_init,
            (str(probe_design), args.tag_length, max_mismatches, tuple(sorted(unresolved)), tuple(sink_tags), (getattr(args, "leading_anchor", None) or "").upper()),
        )
        samples_with_vote: Counter = Counter()
        for votes in fuzzy_per_fastq.values():
            for key in votes:
                samples_with_vote[key] += 1
        for marker in unresolved:
            support = {key: fuzzy_votes.get((marker, key), 0) for key in windows_by_marker[marker]}
            ranked = sorted(support.items(), key=lambda item: (-item[1], item[0]))
            best_key, best = ranked[0]
            runner = ranked[1][1] if len(ranked) > 1 else 0
            samples = samples_with_vote.get((marker, best_key), 0)
            if best == 0:
                excluded.append({"marker": marker, "reason": "no_reads_within_identity"})
                continue
            if best < args.fuzzy_min_reads or samples < args.fuzzy_min_samples:
                excluded.append(
                    {"marker": marker, "reason": "fuzzy_support_below_threshold",
                     "reads": best, "samples": samples}
                )
                continue
            if runner > 0 and best < args.dominance_ratio * runner:
                excluded.append(
                    {"marker": marker, "reason": "ambiguous_windows", "best": best, "runner_up": runner}
                )
                continue
            # The fuzzy route is the risky tier — wrong anchors tend to win by
            # large margins, so the route itself is the review flag.
            flags = ["fuzzy_anchor"]
            if runner > 0 and best < 2 * args.dominance_ratio * runner:
                flags.append("close_runner_up")
            if marker in paralog_flagged:
                flags.append("paralog_shared_candidates")
            resolved[marker] = {
                "route": "fuzzy", "window": best_key, "support": best,
                "runner_up": runner, "supporting_samples": samples, "flags": flags,
            }
    elif unresolved:
        for marker in unresolved:
            excluded.append({"marker": marker, "reason": "no_exact_support_fuzzy_disabled"})


def _infer_anchor_and_build(
    args: argparse.Namespace,
    windows_by_marker: Mapping[str, Mapping[Tuple[str, int], Tuple[str, str, int]]],
    resolved: Dict[str, Dict[str, object]],
    excluded: List[Dict[str, object]],
) -> Tuple[List[Locus], Dict[str, int]]:
    """Stage 4: assign per-locus anchors and build the panel.  The anchor is
    the unknowable left primer's substitute in downstream generation, and
    it must not contain the variant (Ref and Alt reads have to share it),
    so a window whose variant sits at position zero can never be anchored
    and is excluded rather than guessed (removed from ``resolved``)."""
    anchor_input: Dict[str, Tuple[str, int]] = {}
    for marker, record in sorted(resolved.items()):
        ref_window, _alt_window, offset = windows_by_marker[marker][record["window"]]
        if offset == 0:
            excluded.append({"marker": marker, "reason": "variant_at_tag_start"})
            del resolved[marker]
            continue
        anchor_input[marker] = (ref_window, offset)
    anchor_lengths = choose_anchor_lengths(
        anchor_input, args.anchor_min_length, args.anchor_max_length
    )
    loci: List[Locus] = []
    for marker in sorted(resolved):
        record = resolved[marker]
        ref_window, alt_window, offset = windows_by_marker[marker][record["window"]]
        loci.append(
            Locus(
                index=len(loci),
                marker=marker,
                left_primer=ref_window[: anchor_lengths[marker]],
                right_primer="",
                ref_base=ref_window[offset],
                alt_base=alt_window[offset],
                variant_pos=offset,
                ref_tag=ref_window,
                alt_tag=alt_window,
            )
        )
    return loci, anchor_lengths


def _infer_report(
    args: argparse.Namespace,
    probe_design: Path,
    fastq_directory: Path,
    fastqs: Sequence[Path],
    max_mismatches: int,
    targets: Mapping[str, str],
    windows_by_marker: Mapping[str, Mapping[Tuple[str, int], Tuple[str, str, int]]],
    resolved: Mapping[str, Dict[str, object]],
    anchor_lengths: Mapping[str, int],
    loci: Sequence[Locus],
    excluded: Sequence[Dict[str, object]],
) -> Dict[str, object]:
    return {
        "pipeline": "MADCATS infer-design",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "probe_design": str(probe_design.resolve()),
            "fastq_directory": str(fastq_directory.resolve()),
            "fastqs": len(fastqs),
        },
        "parameters": {
            "tag_length": args.tag_length,
            "minimum_identity": args.min_identity,
            "maximum_mismatches": max_mismatches,
            "dominance_ratio": args.dominance_ratio,
            "fuzzy_min_reads": args.fuzzy_min_reads,
            "fuzzy_min_samples": args.fuzzy_min_samples,
            "anchor_min_length": args.anchor_min_length,
            "anchor_max_length": args.anchor_max_length,
        },
        "summary": {
            "markers_in_design": len(targets),
            "resolved": len(loci),
            "resolved_exact": sum(1 for r in resolved.values() if r["route"] == "exact"),
            "resolved_fuzzy": sum(1 for r in resolved.values() if r["route"] == "fuzzy"),
            "flagged_for_review": sum(1 for r in resolved.values() if r["flags"]),
            "excluded": len(excluded),
            "excluded_by_reason": dict(Counter(entry["reason"] for entry in excluded)),
        },
        "per_locus": [
            {
                "marker": marker,
                "route": record["route"],
                "strand": record["window"][0],
                "window_start": record["window"][1],
                "variant_position": windows_by_marker[marker][record["window"]][2] + 1,
                "support": record["support"],
                "runner_up_support": record["runner_up"],
                "supporting_samples": record["supporting_samples"],
                "anchor_length": anchor_lengths[marker],
                "flags": record["flags"],
            }
            for marker, record in sorted(resolved.items())
        ],
        "excluded_loci": sorted(excluded, key=lambda entry: str(entry["marker"])),
    }


def run_infer_design(args: argparse.Namespace) -> int:
    """Infer a reference-tag panel from the probe design plus FASTQ evidence.

    For when the primer file is unavailable.  The tag must be a tag-length
    window of the expanded target containing the SNP, so inference selects
    among a closed candidate set.  Exact read-prefix votes decide first and
    are authoritative; an identity-envelope fuzzy pass runs only for loci
    with zero exact support and can never overrule exact evidence.  The
    output is an ordinary reference-tags CSV plus a QC report recording each
    locus's route, support, margin, and flags for review.
    """
    probe_design = Path(args.probe_design)
    targets = read_probe_targets(probe_design)
    fastq_directory = Path(args.fastq_dir)
    fastqs = discover_fastqs(fastq_directory, args.fastq_glob)
    if args.anchor_max_length < args.anchor_min_length:
        raise MadcError("--anchor-max-length must be at least --anchor-min-length")
    max_mismatches = math.floor((1.0 - args.min_identity) * args.tag_length + 1e-9)

    windows_by_marker, excluded, paralog_flagged = _infer_candidate_windows(targets, args.tag_length)
    resolve_leading_anchor(
        args, fastqs,
        {window[:LEADING_ANCHOR_PREFIX] for windows in windows_by_marker.values() for ref_window, alt_window, _offset in windows.values() for window in (ref_window, alt_window)},
        "infer-design",
    )
    print(
        f"infer-design: exact pass over {len(fastqs)} FASTQ(s), "
        f"{sum(len(w) for w in windows_by_marker.values()):,} candidate windows "
        f"across {len(windows_by_marker):,} markers",
        file=sys.stderr,
    )
    exact_votes, _ = _run_fastq_scans(
        fastqs, args.jobs, _infer_exact_scan, _infer_exact_init, (str(probe_design), args.tag_length, (getattr(args, "leading_anchor", None) or "").upper())
    )
    resolved, unresolved = _infer_resolve_exact(
        windows_by_marker, exact_votes, args.dominance_ratio, paralog_flagged, excluded
    )
    _infer_resolve_fuzzy(
        args, probe_design, fastqs, max_mismatches, windows_by_marker, resolved, unresolved, paralog_flagged, excluded
    )
    loci, anchor_lengths = _infer_anchor_and_build(args, windows_by_marker, resolved, excluded)
    output = Path(args.output)
    write_reference_tags(output, loci)
    if loci:
        load_reference_tags(output, tag_length=args.tag_length)  # round-trip self-check

    report = _infer_report(
        args, probe_design, fastq_directory, fastqs, max_mismatches, targets,
        windows_by_marker, resolved, anchor_lengths, loci, excluded,
    )
    report_path = Path(args.report) if args.report else output.with_suffix(output.suffix + ".qc.json")
    _write_json(report_path, report)
    summary = report["summary"]
    print(
        f"Inferred {summary['resolved']:,} loci "
        f"({summary['resolved_exact']:,} exact, {summary['resolved_fuzzy']:,} fuzzy, "
        f"{summary['flagged_for_review']:,} flagged) to {output}; "
        f"{summary['excluded']:,} excluded — details in {report_path}",
        file=sys.stderr,
    )
    return 0


def positive_int(value: str) -> int:
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def nonnegative_int(value: str) -> int:
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be at least 0")
    return parsed


def probability(value: str) -> float:
    parsed = float(value)
    if not 0 < parsed <= 1:
        raise argparse.ArgumentTypeError("must be greater than 0 and at most 1")
    return parsed


def proper_fraction(value: str) -> float:
    parsed = float(value)
    if not 0 <= parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 0 and below 1")
    return parsed


def partial_read_mode(value: str) -> str:
    """``--partial-read-mode``: on/off, still accepting the 0.12 spelling ``assign``."""
    if value == "assign":
        print("warning: --partial-read-mode assign is now spelled 'on'", file=sys.stderr)
        return "on"
    return value


def ratio_at_least_one(value: str) -> float:
    parsed = float(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be at least 1")
    return parsed


def _default_jobs() -> int:
    return max(1, min(4, os.cpu_count() or 1))


def _add_jobs_argument(parser: argparse.ArgumentParser, what: str = "samples") -> None:
    parser.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(), help=f"{what} processed in parallel (default: up to 4)"
    )


def _add_fastq_arguments(parser: argparse.ArgumentParser, directory_help: str = "folder containing one FASTQ(.gz) per sample") -> None:
    parser.add_argument("--fastq-dir", required=True, help=directory_help)
    parser.add_argument("--fastq-glob", default="*", help="FASTQ filename glob within --fastq-dir (default: *)")


def _add_design_command(subparsers) -> None:
    design = subparsers.add_parser(
        "design",
        help="reconstruct each marker's Ref/Alt tags from DArT primer and probe-design files",
        description=(
            "Reconstruct each marker's canonical Ref and Alt tags from a DArT primer-sequence CSV and "
            "probe-design file, and write them as a reference-tag CSV for generate. Uses no sequencing data."
        ),
    )
    design.add_argument("--primers", required=True, help="DArT primer-sequence CSV")
    design.add_argument(
        "--probe-design", required=True, help="probe-design TSV or CSV with one bracketed variant per marker, e.g. [C/T]"
    )
    design.add_argument(
        "--tag-length",
        type=positive_int,
        default=TAG_LENGTH_DEFAULT,
        help=f"tag length to build, e.g. 54, 81 or 109 (default: {TAG_LENGTH_DEFAULT})",
    )
    design.add_argument("--output", required=True, help="reference-tag CSV to write")
    design.set_defaults(func=run_design)


def _add_infer_design_command(subparsers) -> None:
    infer = subparsers.add_parser(
        "infer-design",
        help="infer the reference-tag panel from the probe design and the reads when the primer file is unavailable",
        description=(
            "Infer each marker's Ref and Alt tags from the probe design and the reads, for panels whose primer "
            "file is not available. Reads vote for each locus's tag window: exact prefix matches decide first, "
            "and a fuzzy pass runs only for loci with no exact votes. Loci that cannot be resolved are left out "
            "and listed, with the reason, in the QC report."
        ),
    )
    infer.add_argument("--probe-design", required=True, help="probe-design TSV or CSV with bracketed variants")
    _add_fastq_arguments(infer)
    infer.add_argument(
        "--leading-anchor",
        default=None,
        help="constant sequence preceding the amplicon in raw deliveries, or 'auto' to detect it (see generate --leading-anchor); reads without it do not vote",
    )
    infer.add_argument(
        "--tag-length",
        type=positive_int,
        default=TAG_LENGTH_DEFAULT,
        help="reported tag length of the panel, e.g. 54, 81, or 109 (default: 81)",
    )
    infer.add_argument(
        "--min-identity",
        type=probability,
        default=0.80,
        help="minimum identity for fuzzy-pass votes, as in generate (default: 0.80)",
    )
    infer.add_argument(
        "--dominance-ratio",
        type=ratio_at_least_one,
        default=10.0,
        help="a winning window must have at least this many times the runner-up's votes (default: 10)",
    )
    infer.add_argument(
        "--fuzzy-min-reads",
        type=nonnegative_int,
        default=30,
        help="fuzzy-pass votes needed to resolve a locus with no exact votes; 0 disables the fuzzy pass (default: 30)",
    )
    infer.add_argument(
        "--fuzzy-min-samples",
        type=nonnegative_int,
        default=2,
        help="distinct samples that must contribute fuzzy votes (default: 2)",
    )
    infer.add_argument(
        "--anchor-min-length",
        type=positive_int,
        default=15,
        help="shortest pseudo-primer anchor to try (default: 15)",
    )
    infer.add_argument(
        "--anchor-max-length",
        type=positive_int,
        default=30,
        help="longest pseudo-primer anchor; an anchor always ends before the variant (default: 30)",
    )
    _add_jobs_argument(infer)
    infer.add_argument("--output", required=True, help="reference-tag CSV to write")
    infer.add_argument("--report", help="inference QC JSON (default: OUTPUT.qc.json)")
    infer.set_defaults(func=run_infer_design)


def _add_merge_pairs_command(subparsers) -> None:
    merge = subparsers.add_parser(
        "merge-pairs",
        help="merge overlapping paired-end reads into single fragment reads",
        description=(
            "Merge each overlapping R1/R2 pair into one fragment read, before normalize-reads or generate. R2 is "
            "reverse-complemented and the best-scoring overlap is used (overlap length with a penalty for "
            "each mismatch). Adapter read-through falls outside the "
            "overlap and never reaches the merged read, so no adapter trimming is needed. Use this only when "
            "mates overlap."
        ),
    )
    merge.add_argument("--fastq-dir", required=True, help="folder holding paired FASTQ(.gz) files")
    merge.add_argument("--r1-suffix", default="_R1.fastq.gz", help="filename suffix of the R1 mate (default: _R1.fastq.gz)")
    merge.add_argument("--r2-suffix", default="_R2.fastq.gz", help="filename suffix of the R2 mate (default: _R2.fastq.gz)")
    merge.add_argument(
        "--min-overlap", type=positive_int, default=15,
        help="shortest overlap, in bases, that lets a pair merge (default: 15)",
    )
    merge.add_argument(
        "--max-mismatch-rate", type=float, default=0.1,
        help="mismatch fraction tolerated inside the overlap (default: 0.1)",
    )
    merge.add_argument(
        "--unmerged", choices=("r1", "pool", "drop"), default="r1",
        help=(
            "pairs that do not merge: keep the R1 (default), keep both mates (pool; for dual-anchored "
            "designs whose mates cover opposite halves of the locus), or drop the pair"
        ),
    )
    _add_jobs_argument(merge)
    merge.add_argument("--out-dir", required=True, help="folder for merged per-sample FASTQs")
    merge.set_defaults(func=run_merge_pairs)


def _add_window_scan_command(subparsers) -> None:
    scan = subparsers.add_parser(
        "window-scan",
        help="try several tag lengths on a subset of samples and recommend one",
        description=(
            "Run normalize-reads and generate at each candidate tag length on a subset of samples, measure "
            "read retention, cleanliness and haplotype information, and recommend a length. Use it before "
            "the first normalize-reads run of a new panel."
        ),
    )
    scan.add_argument("--probe-design", required=True, help="probe-design TSV or CSV with bracketed variants")
    _add_fastq_arguments(scan)
    scan.add_argument(
        "--lengths", required=True,
        help="comma-separated candidate tag lengths to run in full, e.g. 41,61,81,110",
    )
    scan.add_argument(
        "--scan-seed", type=int, default=None,
        help="draw the --scan-samples subset at random with this seed instead of taking the first files in sorted order",
    )
    scan.add_argument(
        "--generate-args", default="",
        help="extra generate flags for every candidate run, quoted as one string (e.g. \"--min-supporting-freq 0.1\")",
    )
    scan.add_argument(
        "--scan-samples", type=positive_int, default=6,
        help="number of samples each candidate runs on: the first in sorted order, or a --scan-seed draw (default: 6)",
    )
    scan.add_argument(
        "--scan-draws", type=positive_int, default=1,
        help=(
            "independent sample draws to sweep (seeds --scan-seed, +1, +2, ...; default: 1); with more than one, "
            "each draw's recommendation and their agreement are reported, so you can see whether the "
            "recommendation depends on which samples were scanned"
        ),
    )
    scan.add_argument(
        "--grid-step", type=positive_int, default=10,
        help="step for the dense predicted-retention grid from the coverage pass (default: 10)",
    )
    scan.add_argument(
        "--info-band", type=probability, default=0.10,
        help=(
            "survival tolerance for the information tiebreak: candidates within this fraction of "
            "the best survival (retention x assignment x cleanliness) are eligible, and the most "
            "informative eligible length is recommended (default: 0.10)"
        ),
    )
    scan.add_argument(
        "--anchor-length", type=positive_int, default=12,
        help="as in normalize-reads; used for the coverage pass and every candidate run (default: 12)",
    )
    scan.add_argument(
        "--kmer-size", type=positive_int, default=21,
        help="as in normalize-reads; used for the coverage pass and every candidate run (default: 21)",
    )
    scan.add_argument(
        "--str-penalty", type=float, default=0.002,
        help="as in normalize-reads; used for the coverage pass and every candidate run (default: 0.002)",
    )
    scan.add_argument(
        "--group-routing", choices=("nearest", "covered"), default="nearest",
        help="passed to the candidate normalize-reads runs; use covered for merged, tiled input (default: nearest)",
    )
    _add_jobs_argument(scan)
    scan.add_argument("--out-dir", required=True, help="folder for per-length runs and the scan report")
    scan.set_defaults(func=run_window_scan)


def _add_normalize_reads_command(subparsers) -> None:
    normalize = subparsers.add_parser(
        "normalize-reads",
        help="orient and window reads that do not start at a primer (e.g. AgriSeq), ready for generate",
        description=(
            "For platforms such as AgriSeq (Ion Torrent) whose reads do not start at a fixed primer. Each read "
            "is located against the probe-design targets, oriented forward, and cut to a fixed-length window "
            "containing the variant. Writes one FASTQ per sample plus a window reference-tags CSV, ready for "
            "generate. Pass the first run's CSV to later runs with --reuse-tags."
        ),
    )
    normalize.add_argument("--probe-design", required=True, help="probe-design TSV or CSV with bracketed variants")
    _add_fastq_arguments(normalize, directory_help="folder containing one raw FASTQ(.gz) per sample")
    normalize.add_argument(
        "--tag-length",
        type=positive_int,
        default=41,
        help="window length to cut; shorter windows lose fewer reads to indels on indel-prone chemistries (default: 41)",
    )
    normalize.add_argument(
        "--window-placement",
        choices=("auto", "center"),
        default="auto",
        help=(
            "auto learns each marker's window from observed read coverage (extra pass); center needs no read data "
            "(default: auto).  With --reuse-tags the windows are fixed, and auto instead runs the pass to score "
            "the frozen windows against this run's reads (freeze_health in the QC); center skips it"
        ),
    )
    normalize.add_argument(
        "--anchor-length",
        type=positive_int,
        default=12,
        help=(
            "minimum assignment-anchor prefix length; anchors grow per locus only where window "
            "prefixes collide, and the variant is always kept out of the anchor (default: 12)"
        ),
    )
    normalize.add_argument(
        "--reuse-tags",
        help=(
            "window reference-tags CSV from a previous normalize-reads run to reuse verbatim, freezing "
            "window placement across runs of the same panel so microhaplotype sequences stay comparable "
            "(--anchor-length is ignored; --window-placement auto scores the frozen windows, center skips that)"
        ),
    )
    normalize.add_argument(
        "--group-routing",
        choices=("nearest", "covered"),
        default="nearest",
        help=(
            "how a read covering several windows of one LocusGroup is routed: 'nearest' (default) sends it "
            "to the one tile whose window center is nearest the read's center; 'covered' cuts a window at "
            "every tile the read fully covers. Use covered with merged full-span fragments: they all share "
            "the amplicon's center, so with nearest one tile would get almost every read. In both modes "
            "each window counts each fragment at most once"
        ),
    )
    normalize.add_argument(
        "--str-penalty",
        type=float,
        default=0.002,
        help=(
            "auto placement: coverage penalty per tandem-repeat bp (unit sizes 1-3) beyond 4 in a candidate "
            "window; steers windows away from repeats where reads carrying one allele are lost more often "
            "(0 disables; default: 0.002)"
        ),
    )
    normalize.add_argument(
        "--kmer-size",
        type=positive_int,
        default=21,
        help="k-mer length for locating reads against the design targets (default: 21)",
    )
    _add_jobs_argument(normalize)
    normalize.add_argument("--out-dir", required=True, help="folder for the normalized per-sample FASTQs")
    normalize.add_argument("--output-tags", help="window reference-tags CSV (default: OUT_DIR/window_reference_tags.csv)")
    normalize.add_argument("--report", help="QC/provenance JSON (default: OUT_DIR/normalize_reads.qc.json)")
    normalize.set_defaults(func=run_normalize_reads)


def _add_generate_command(subparsers) -> None:
    """``generate``'s options, grouped by the pipeline stage each one feeds."""
    generate = subparsers.add_parser(
        "generate",
        help="count alleles in every sample and write the MADC table",
        description=(
            "Assign reads to markers, count exact tag sequences per sample, report new microhaplotypes that "
            "pass the evidence rules, remove sequencing errors by denoising, and write a raw MADC table plus "
            "a QC/provenance JSON. Run one generate over the whole cohort: the reproducibility route counts "
            "supporting samples within the run. The raw MADC goes through HapApp before BIGr and other "
            "downstream tools."
        ),
    )
    panel = generate.add_argument_group("panel and samples (stage 1)")
    _add_fastq_arguments(panel)
    panel.add_argument("--primers", help="DArT primer-sequence CSV (use with --probe-design)")
    panel.add_argument("--probe-design", help="DArT probe-design TSV or CSV (use with --primers)")
    panel.add_argument(
        "--reference-tags",
        help=(
            "reference-tag CSV from design, infer-design or normalize-reads, or curated by hand; "
            "use instead of --primers/--probe-design"
        ),
    )
    panel.add_argument(
        "--target-metadata",
        nargs="*",
        help="target-manifest CSV(s); defaults to targets_*.csv in the FASTQ folder",
    )
    panel.add_argument(
        "--sample-name-field",
        default="targetid",
        help="unique metadata field used for sample columns (default: targetid/FASTQ stem)",
    )
    panel.add_argument(
        "--sample-lookup",
        help=(
            "two-column CSV/TSV mapping FASTQ-derived sample ids to sample names; replaces the output "
            "column names (every run sample must be mapped) while the preamble targetid row keeps the "
            "FASTQ ids, matching the raw-MADC layout"
        ),
    )
    panel.add_argument(
        "--tag-length",
        type=positive_int,
        default=TAG_LENGTH_DEFAULT,
        help=f"tag length; must match the panel, e.g. 41 after normalize-reads (default: {TAG_LENGTH_DEFAULT})",
    )

    assignment = generate.add_argument_group("read assignment (stage 2)")
    assignment.add_argument(
        "--min-identity",
        type=probability,
        default=0.80,
        help=(
            "minimum identity to the marker's Ref or Alt tag, counting substitutions only; a gapped "
            "alignment uses the same edit limit (default: 0.80)"
        ),
    )
    assignment.add_argument(
        "--leading-anchor",
        default=None,
        help=(
            "constant sequence that precedes the amplicon in every read of a raw delivery (DArT raw "
            "FASTQs carry an inline barcode and the common sequence GATGGATCCAGTGC); each read is cut "
            f"to start right after its first occurrence within the first {LEADING_ANCHOR_SEARCH} bases "
            "beyond the anchor's length, and reads without it are counted as reads_no_leading_anchor; "
            "'auto' detects the anchor from the first FASTQ (default: none; preprocessed deliveries are "
            "primer-started, and a run whose reads carry an undeclared constant prefix stops and names it)"
        ),
    )
    assignment.add_argument(
        "--primer-mismatches",
        type=int,
        choices=(0, 1),
        default=1,
        help=(
            "mismatches allowed when matching a read's start to a primer/anchor (default: 1). With 0, "
            "reads from Match haplotypes whose distinguishing SNP falls inside the anchor are lost"
        ),
    )
    assignment.add_argument(
        "--gapped-rescue",
        choices=("on", "off"),
        default="on",
        help=(
            "give reads that fail the substitution-only identity check for every marker one gapped "
            "alignment under the same edit limit, so undesigned indels are counted instead of dropped "
            "(default: on)"
        ),
    )
    assignment.add_argument(
        "--max-indel-length",
        type=positive_int,
        default=GAP_BAND_DEFAULT,
        help=f"widest gap (alignment band) a gapped alignment may use (default: {GAP_BAND_DEFAULT})",
    )
    assignment.add_argument(
        "--partial-read-mode",
        type=partial_read_mode,
        choices=("on", "off"),
        default="on",
        help=(
            "handling of reads shorter than --tag-length: on adds such a read to the one full-length "
            "allele whose start it matches exactly, only when that allele already passes the evidence rules "
            "and the read is longer than every prefix shared between that locus's alleles, so partial reads "
            "cannot create an allele or change an allele ratio (default: on); off drops them"
        ),
    )
    assignment.add_argument(
        "--min-partial-length",
        type=positive_int,
        default=None,
        help=(
            "minimum length for a partial read to be counted with --partial-read-mode on "
            "(default: --tag-length minus 9)"
        ),
    )
    assignment.add_argument(
        "--cache-dir",
        help=(
            "folder for a reusable per-sample cache, so re-runs skip the read pass. The cache is rebuilt "
            "when the panel or --tag-length, --min-identity, --leading-anchor, --primer-mismatches, "
            "--gapped-rescue, --max-indel-length or --min-partial-length changes"
        ),
    )
    _add_jobs_argument(assignment)

    evidence = generate.add_argument_group("evidence (stage 4: which new sequences are reported)")
    evidence.add_argument(
        "--min-global-count",
        type=positive_int,
        default=1,
        help="minimum project-wide reads for a newly discovered allele (default: 1)",
    )
    evidence.add_argument(
        "--reproducible-min-reads",
        type=positive_int,
        default=10,
        help="reproducibility route: reads per sample counted toward --reproducible-min-samples (default: 10)",
    )
    evidence.add_argument(
        "--reproducible-min-samples",
        type=nonnegative_int,
        default=2,
        help=(
            "reproducibility route: a new allele needs at least --reproducible-min-reads reads in at "
            "least this many samples; 0 disables the route (default: 2)"
        ),
    )
    evidence.add_argument(
        "--depth-min-reads",
        type=nonnegative_int,
        default=30,
        help=(
            "depth route: reads within one sample that are enough on their own to report a new allele, so "
            "alleles found in only one sample can still be reported when samples are not replicated; "
            "0 disables the route (default: 30)"
        ),
    )
    evidence.add_argument(
        "--discovery-min-quality",
        type=int,
        choices=QUALITY_THRESHOLDS,
        default=20,
        help="Phred score every tag base must reach for a read to count as high-quality on the quality route (default: 20)",
    )
    evidence.add_argument(
        "--quality-min-clean-reads",
        type=nonnegative_int,
        default=0,
        help=(
            "quality route: project-wide high-quality reads that are enough on their own to report a new "
            "allele; 0 disables the route (default: 0). Off by default because one systematically "
            "low-quality base position would block real alleles; reproducibility is the main route"
        ),
    )
    evidence.add_argument(
        "--min-supporting-freq",
        type=proper_fraction,
        default=0.0,
        help=(
            "minimum share of a sample's reads at the locus before that sample can count toward the "
            "reproducibility or depth route; 0 disables the floor (default: 0). Use 0.10 at very high depth "
            "(above ~500x per sample per marker, e.g. AgriSeq), or about 0.01 for unreplicated designs that "
            "rely on the depth route. Otherwise leave it off: co-amplified paralogs lower real alleles' "
            "shares, so a floor can drop them"
        ),
    )

    denoising = generate.add_argument_group("denoising (stage 5)")
    denoising.add_argument(
        "--denoise-distance",
        type=int,
        choices=(0, 1, 2),
        default=1,
        help=(
            "per-sample denoising: the most differing bases allowed between a rare sequence and the more "
            "abundant parent it is merged into; 0 disables (default: 1)"
        ),
    )
    denoising.add_argument(
        "--global-denoise-distance",
        type=int,
        choices=(0, 1, 2),
        default=1,
        help=(
            "pooled same-allele merge: a sequence that failed the evidence rules and differs at no more "
            "than this many bases from a protected parent carrying the same target allele is merged into "
            "it in every sample; 0 disables (default: 1)"
        ),
    )
    denoising.add_argument(
        "--denoise-alpha",
        type=float,
        default=2.0,
        help="UNOISE-style abundance-skew exponent: a parent must be at least 2^(alpha*d + 1) times more abundant (default: 2.0)",
    )

    output = generate.add_argument_group("output (stages 6-7)")
    output.add_argument("--output", required=True, help="MADC CSV to write")
    output.add_argument("--qc-output", help="QC/provenance JSON (default: OUTPUT.qc.json)")
    output.add_argument(
        "--annotations",
        choices=("on", "off"),
        default="on",
        help=(
            "annotate every row in two otherwise unused statistics columns and a sidecar file: PICSnp becomes "
            "Variants (the row's events against the Ref tag in 1-based tag coordinates, e.g. "
            "55:C>CAAAA;91:A>T;clip=4), AvgPIC becomes MarkerType ('.', 'snp', 'indel' or 'both' relative "
            "to the nearest canonical), and OUTPUT.variants.csv lists one event per line; 'off' keeps the "
            "exact vendor header and writes no sidecar (default: on)"
        ),
    )
    output.add_argument(
        "--number-match-ids",
        action="store_true",
        help="append a stable serial to each RefMatch/AltMatch/Other ID (e.g. |RefMatch_0001) so every AlleleID is unique",
    )
    output.add_argument(
        "--dart-preamble",
        action="store_true",
        help=(
            "prepend the vendor-style seven-row starred preamble (service label, plate rows, and a "
            "targetid row with the FASTQ-derived sample ids) that HapApp expects"
        ),
    )
    output.add_argument(
        "--preamble-project",
        help="project/service label for the first preamble row (default: the output filename before its first dot)",
    )
    generate.set_defaults(func=run_generate)


def _add_validate_command(subparsers) -> None:
    validate = subparsers.add_parser(
        "validate",
        help="compare a generated MADC with a reference MADC",
        description=(
            "Compare a generated MADC with a reference MADC, usually the vendor's delivery for the same "
            "FASTQs: row precision and recall by allele category, classification agreement on shared rows, "
            "and count agreement across shared cells. The report is printed as JSON. The reference MADC is "
            "used only for this comparison."
        ),
    )
    validate.add_argument("--generated", required=True, help="MADC table written by generate")
    validate.add_argument(
        "--expected",
        required=True,
        help="reference MADC to compare against, with or without the vendor's starred preamble rows",
    )
    validate.add_argument("--report", help="also write the JSON report to this file")
    validate.set_defaults(func=run_validate)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="madcats",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "MADCATS (Marker-Anchored Discovery and Counting of Alleles from Targeted\n"
            "Sequencing) turns targeted-sequencing FASTQs into raw MADC tables: the\n"
            "microhaplotypes at each marker and their read counts in every sample. It\n"
            "counts the alleles each panel was designed for and discovers new ones."
        ),
        epilog=(
            "typical workflow:\n"
            "  1. build the panel           design, or infer-design without a primer file\n"
            "  2. prepare reads if needed   merge-pairs, window-scan, normalize-reads\n"
            "  3. count                     generate\n"
            "  4. clean (optional)          filter-paralogs, with a reference genome\n"
            "  5. check (optional)          validate\n"
            "  6. genotype (optional)       call, from the MADC to VCFs\n"
            "\n"
            "Run 'madcats <command> --help' for a command's options.\n"
            "Documentation: https://madcats.readthedocs.io"
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_design_command(subparsers)
    _add_infer_design_command(subparsers)
    _add_merge_pairs_command(subparsers)
    _add_window_scan_command(subparsers)
    _add_normalize_reads_command(subparsers)
    _add_generate_command(subparsers)
    _add_validate_command(subparsers)
    import madcats_paralogs  # here, not at the top: it imports this module

    madcats_paralogs.add_filter_paralogs_command(subparsers)
    import madcats_call

    madcats_call.add_call_command(subparsers)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (MadcError, OSError, csv.Error) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    # Run as a script, this file is __main__; madcats_paralogs's "import madcats"
    # must find it rather than load a second copy with its own MadcError.
    sys.modules.setdefault("madcats", sys.modules[__name__])
    raise SystemExit(main())
