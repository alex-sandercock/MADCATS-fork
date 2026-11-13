"""``madcats filter-paralogs``: clean a finished MADC with a reference genome.

In panels where paralogous copies co-amplify with a target, a paralog read
that differs from the Ref tag only at non-target bases becomes a RefMatch or
AltMatch row of that marker.  This module places each marker on the genome
by its Ref tag and, for CHROM_POS CloneIDs, by that position (see
``--clone-positions``); aligns every haplotype; and drops or flags the rows
whose best alignment lies at another locus (paralogs) or ties between loci
(ambiguous).

It is kept apart from ``madcats.py`` because it is the only command with a
dependency outside the standard library: the ``mappy`` aligner (minimap2),
imported on first use (``pip install 'madcats[genome]'``), so ``madcats``
runs without it.  The MADC reading helpers it shares with ``validate``
(header detection, sample columns, allele categories, the count-overflow
check) stay in ``madcats``.

What is minimap2's and what is this module's.  Alignment is minimap2's,
through its Python binding mappy (MIT licence; https://github.com/lh3/minimap2):

    Li H (2018) Minimap2: pairwise alignment for nucleotide sequences.
        Bioinformatics 34(18):3094-3100.  doi:10.1093/bioinformatics/bty191
    Li H (2021) New strategies to improve minimap2 alignment accuracy.
        Bioinformatics 37(23):4572-4574.  doi:10.1093/bioinformatics/btab705

The index (k=21, w=11) and the alignment scores are those of minimap2's
short-read preset (``-x sr``).  Three of the preset's mapping options are
changed for the genome-wide search (``GENOME_WIDE_MAP_OPTIONS``).  The layout
of an index file read by ``_index_part_contigs`` is that written by
minimap2's ``mm_idx_dump`` (index.c).  Everything else is MADCATS's own, with
no published counterpart: the score recomputed from an alignment
(``hit_score``), the rule for when two alignments are one locus, the
realignment of every haplotype at its own marker's site with short seeds,
marker placement with the CloneID position, and the labels.  The label rule
began as "a read whose best alignment score exceeds its second best", as
bwa mem reports them, in a prototype that is not part of MADCATS; here the
comparison is between the marker's own site and anywhere else.  The
measurements behind the choices are in docs/DEVELOPMENT.md ("Genome paralog
filter").

Reading map: scoring and alignment (``hit_score``, ``GenomeAligner``,
``WindowIndex``, ``map_sequences``); marker placement and labels
(``place_marker``, ``score_at_sites``, ``genome_label``); the two passes over
the table (``_read_haplotypes``, ``_write_filtered_table``) and the function
that joins them (``filter_paralogs``); the command line
(``add_filter_paralogs_command``, ``run_filter_paralogs``).
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import datetime as dt
import re
import struct
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Mapping, NamedTuple, Optional, Sequence, Set, Tuple

from madcats import (
    CountOverflowCheck,
    MadcError,
    __version__,
    _default_jobs,
    _write_json,
    allele_type,
    is_madc_header,
    madc_identity_columns,
    madc_read_counts,
    madc_sample_start,
    open_text,
    positive_int,
)

# In flag mode the label goes into the column of this vendor statistic, which is renamed.
REPLACED_STATISTICS_COLUMN = "PICRef"
GENOME_LABEL_COLUMN = "GenomeLabel"

# The alignment scores of minimap2's short-read preset (-x sr): match, mismatch, gap open, gap
# extend, long-gap open, long-gap extend; and minimap2's general score for an ambiguous base.
# They are passed to the aligner, although they are its preset's values, so that hit_score and
# the aligner cannot come to differ.
SHORT_READ_SCORING = (2, 8, 12, 2, 24, 1, 1)
# The index of that preset: minimizers of 21 bases in windows of 11.
SHORT_READ_INDEX_K_W = (21, 11)
# Changes to the preset's mapping options for the genome-wide search (the index is unchanged):
# chains of a single seed (the preset needs 2), a lower least chain score (25) and up to 50
# secondary hits (20), so 51 hits in all.  An 81-base allele with several SNPs then still finds
# its copies.  Chosen on the lettuce MADCATS MADC, where the rows kept or dropped agreed with a
# bwa mem classification for 98.2% of rows, against 97.3% with the preset's values
# (docs/DEVELOPMENT.md).  The realignment at a marker's site keeps the preset's values.
GENOME_WIDE_MAP_OPTIONS = {"min_cnt": 1, "min_chain_score": 20, "best_n": 50}
# Distinct genome sites remembered per haplotype, best first.  Ref tags keep all of theirs.
MAX_SITES_PER_HAPLOTYPE = 8
# Sequences handed to a worker thread at a time.
MAPPING_BATCH = 2000

# A CloneID of the form CHROM_POS; the position may be zero-padded.
CLONE_POSITION_RE = re.compile(r"(.+)_0*(\d+)")
# A CloneID position up to this many bases outside a site still counts as at that site.
CLONE_POSITION_SLACK = 10
# More than this share of CloneID positions away from their Ref tags' sites brings a warning
# that the positions may refer to another assembly.
CLONE_DISAGREEMENT_WARNING_SHARE = 0.1

# The realignment of a haplotype at its own marker's site: the site with this many genome bases
# either side, indexed with minimizers of 11 bases in windows of 3, so that an allele too
# diverged to seed in the genome-wide search (21 bases) still aligns where it belongs.  The
# rest of the genome is searched with the longer seeds only, so a copy elsewhere that is as
# diverged as that may go unfound.
SITE_FLANK = 20
SITE_INDEX_K, SITE_INDEX_W = 11, 3

# CIGAR operations as mappy numbers them.
CIGAR_INSERTION, CIGAR_DELETION = 1, 2
# A minimap2 index file: the bytes each part starts with, and the flag of a part stored
# without its sequences (MM_I_NO_SEQ).
INDEX_PART_MAGIC = b"MMI\x02"
INDEX_WITHOUT_SEQUENCES = 2

# The labels sidecar's Action column holds what was done with the row: kept, flagged,
# protected or dropped.
LABELS_SIDECAR_COLUMNS = [
    "CloneID", "AlleleID", "AlleleSequence", "Label", "Action",
    "TargetScore", "TargetFrom", "OffTargetScore", "OffTargetSites", "OffTargetHit",
]
MARKERS_SIDECAR_COLUMNS = ["CloneID", "Placement", "RefTagSites", "CloneIDCheck", "Sites"]


class GenomeHit(NamedTuple):
    """One alignment to the genome: a 0-based, half-open interval, its strand (1 or -1) and its score.

    Also used for a genome window that stands in for a site, with strand 0 and score 0.
    """

    contig: str
    start: int
    end: int
    strand: int
    score: int


class GenomePlacement(NamedTuple):
    """A sequence's distinct genome sites, best score first (hits overlapping a better one are merged into it)."""

    sites: Tuple[GenomeHit, ...]

    @property
    def equal_best(self) -> Tuple[GenomeHit, ...]:
        """The sites that share the best score."""
        return tuple(site for site in self.sites if site.score == self.sites[0].score)


class MarkerPlacement(NamedTuple):
    """Where a marker's own tag lies in the genome.

    ``placed_by`` (the markers sidecar's Placement column): ``ref_tag`` (the
    Ref tag has one best site), ``ref_tag+clone_id`` (several equal-best
    sites; the CloneID position picks one), ``ref_tag_multicopy`` (several,
    none picked: every site is the marker's), ``clone_id_site`` (no
    equal-best site is at the CloneID position; the Ref tag's alignment
    within a window around it), ``clone_id`` (the Ref tag did not align; the
    window around the CloneID position stands in) or ``unplaced``.
    ``ref_tag_site_count`` is the number of equal-best genome sites of the Ref
    tag, whichever of them became ``sites``.  ``clone_id_check`` compares the
    CloneID position with the Ref-tag sites: ``agree``, ``disagree``,
    ``no_position`` (the CloneID is not CHROM_POS), ``contig_absent``,
    ``fallback``, or ``ignored`` (``--clone-positions ignore``).
    """

    sites: Tuple[GenomeHit, ...]
    placed_by: str
    ref_tag_site_count: int
    clone_id_check: str


class GenomeLabel(NamedTuple):
    """A haplotype's label and the two scores that decided it.

    The fields after the label are the labels sidecar's columns of the same
    names.  ``target_score`` is the best alignment at the marker's own site
    and ``target_from`` says which search gave it: ``genome`` (the
    genome-wide search), ``site`` (the realignment at the site, which scored
    higher or was the only one to align there) or ``.`` (none).
    ``off_target_score`` is the best alignment anywhere else,
    ``off_target_sites`` the number of distinct sites elsewhere with that
    score, among the at most ``MAX_SITES_PER_HAPLOTYPE`` sites kept, and
    ``off_target_hit`` the first of them.
    """

    label: str
    target_score: Optional[int]
    target_from: str
    off_target_score: Optional[int]
    off_target_sites: int
    off_target_hit: Optional[GenomeHit]


def hit_score(cigar: Sequence[Sequence[int]], edit_distance: int, matches: int, block_length: int) -> int:
    """The score of an alignment under ``SHORT_READ_SCORING``, from what mappy reports of it.

    mappy returns an alignment's CIGAR, edit distance (NM), number of matching
    bases and block length, but not minimap2's own alignment score (AS), so
    the score is worked out here.  As minimap2 defines them (``mm_update_extra``
    in align.c), the block length counts aligned columns and gap bases but
    not columns with an ambiguous base, and NM counts all three kinds of
    difference; hence mismatches = block - matches - gap bases, and ambiguous
    columns = NM - (block - matches).  A gap costs the cheaper of minimap2's
    two affine pieces (Li 2018).  Both differences are held at zero or above,
    since a gap base that is itself ambiguous is not in the block length.
    """
    (match_score, mismatch_penalty, gap_open_penalty, gap_extend_penalty, long_gap_open_penalty,
     long_gap_extend_penalty, ambiguous_penalty) = SHORT_READ_SCORING
    gap_bases = gap_cost = 0
    for length, operation in cigar:
        if operation in (CIGAR_INSERTION, CIGAR_DELETION):
            gap_bases += length
            gap_cost += min(
                gap_open_penalty + gap_extend_penalty * length, long_gap_open_penalty + long_gap_extend_penalty * length
            )
    differences = block_length - matches
    mismatches = max(0, differences - gap_bases)
    ambiguous_bases = max(0, edit_distance - differences)
    return match_score * matches - mismatch_penalty * mismatches - ambiguous_penalty * ambiguous_bases - gap_cost


def _same_locus(hit: GenomeHit, site: GenomeHit) -> bool:
    """Whether two alignments are at one locus: their intervals share at least half of the shorter one.

    Adjacent tandem copies, which share a few bases, stay distinct.
    """
    if hit.contig != site.contig:
        return False
    overlap = min(hit.end, site.end) - max(hit.start, site.start)
    return overlap > 0 and 2 * overlap >= min(hit.end - hit.start, site.end - site.start)


def distinct_sites(
    hits: Iterable[GenomeHit], max_sites: Optional[int] = MAX_SITES_PER_HAPLOTYPE
) -> GenomePlacement:
    """A sequence's distinct genome sites, best first, at most ``max_sites`` of them (None: all).

    A hit at the same locus as a better-scoring one (``_same_locus``) is
    merged into it.
    """
    sites: List[GenomeHit] = []
    for hit in sorted(hits, key=lambda hit: (-hit.score, hit.contig, hit.start)):
        if any(_same_locus(hit, site) for site in sites):
            continue
        sites.append(hit)
        if len(sites) == max_sites:
            break
    return GenomePlacement(tuple(sites))


def _index_part_contigs(path: Path) -> Optional[List[List[str]]]:
    """The contig names in each part of a minimap2 index file, or None when the file is not one.

    minimap2 writes an index in one part for every ``-I`` bases of genome,
    and mappy loads only the first part.  The file is read as minimap2's
    ``mm_idx_dump`` (index.c) writes it.  Each part is: 4 bytes ``MMI\\2``;
    five 32-bit words (window, k-mer length, bucket bits, number of
    sequences, flags); for every sequence its name, preceded by the name's
    length in one byte, and its length in 4 bytes; for every bucket two
    counts, of 8-byte positions and of 16-byte hash entries; and, unless the
    flags say the sequences are left out, the sequences at 4 bits a base in
    32-bit words.
    """
    file_size = path.stat().st_size
    with path.open("rb", buffering=0) as handle:

        def read_exactly(count: int) -> bytes:
            data = handle.read(count)
            if len(data) != count:
                raise MadcError(f"{path} ends inside a minimap2 index part; rebuild the index")
            return data

        def read_word() -> int:
            return struct.unpack("<I", read_exactly(4))[0]

        parts: List[List[str]] = []
        while handle.tell() < file_size:
            if read_exactly(4) != INDEX_PART_MAGIC:
                if not parts:
                    return None
                raise MadcError(f"{path} has unreadable data after its first index part; rebuild the index")
            _, _, bucket_bits, sequence_count, flags = struct.unpack("<5I", read_exactly(20))
            contig_names, total_bases = [], 0
            for _ in range(sequence_count):
                contig_names.append(read_exactly(read_exactly(1)[0]).decode())
                total_bases += read_word()
            for _ in range(1 << bucket_bits):
                handle.seek(8 * read_word(), 1)
                handle.seek(16 * read_word(), 1)
            if not flags & INDEX_WITHOUT_SEQUENCES:
                handle.seek(4 * ((total_bases + 7) // 8), 1)
            parts.append(contig_names)
        if handle.tell() != file_size:
            raise MadcError(f"{path} ends inside a minimap2 index part; rebuild the index")
        return parts


def _genome_hit(hit, contig: Optional[str] = None, offset: int = 0) -> GenomeHit:
    """A mappy alignment as a GenomeHit.  ``contig`` and ``offset`` put a hit within a window on the genome."""
    return GenomeHit(
        hit.ctg if contig is None else contig, offset + hit.r_st, offset + hit.r_en, hit.strand,
        hit_score(hit.cigar, hit.NM, hit.mlen, hit.blen),
    )


def _load_mappy():
    try:
        import mappy
    except ImportError as exc:
        raise MadcError(
            "filter-paralogs needs the optional minimap2 aligner: pip install 'madcats[genome]' (or pip install mappy)"
        ) from exc
    return mappy


class GenomeAligner:
    """The genome-wide search: minimap2's short-read preset through mappy, with
    ``SHORT_READ_SCORING`` and the mapping options of ``GENOME_WIDE_MAP_OPTIONS``."""

    def __init__(self, aligner, mappy_module) -> None:
        self._aligner = aligner
        self._mappy = mappy_module
        self.contigs = frozenset(aligner.seq_names)
        self.description = {
            "mappy": getattr(mappy_module, "__version__", None),
            "preset": "sr",
            "k": aligner.k,
            "w": aligner.w,
            "scoring": list(SHORT_READ_SCORING),
            "map_options": dict(GENOME_WIDE_MAP_OPTIONS),
            "site_realignment": {"k": SITE_INDEX_K, "w": SITE_INDEX_W, "pad": SITE_FLANK},
        }

    @classmethod
    def open(cls, genome: Optional[Path], index: Optional[Path], jobs: int) -> "GenomeAligner":
        """Load ``index`` when it exists; otherwise index ``genome`` (saving to ``index`` when given)."""
        mappy = _load_mappy()
        if index is not None and index.exists():
            source, index_to_write = index, None
        elif genome is not None:
            source, index_to_write = genome, index
        else:
            raise MadcError(f"Genome index not found: {index} (pass --genome to build it)")
        if not source.exists():
            raise MadcError(f"Genome file not found: {source}")
        # The file may be an index whichever option named it; an index is loaded, never rebuilt.
        parts = _index_part_contigs(source)
        if parts is not None:
            index_to_write = None
        if index_to_write is not None:
            index_to_write.parent.mkdir(parents=True, exist_ok=True)
        if parts is not None and len(parts) > 1:
            raise MadcError(
                f"{source} is a minimap2 index in {len(parts)} parts ({sum(map(len, parts)):,} contigs), and mappy "
                f"searches only the first ({len(parts[0]):,} contigs), so the rest of the genome would be missed. "
                "Rebuild it as one part: give --genome and a new --genome-index path, or run minimap2 -x sr "
                "with -I above the genome size"
            )
        verb = "Loading index" if parts is not None else "Indexing genome"
        print(f"{verb} {source}" + (f" (saving {index_to_write})" if index_to_write else ""), file=sys.stderr)
        aligner = mappy.Aligner(
            str(source), preset="sr", n_threads=jobs, fn_idx_out=str(index_to_write) if index_to_write else None,
            scoring=list(SHORT_READ_SCORING), **GENOME_WIDE_MAP_OPTIONS,
        )
        if not aligner:
            raise MadcError(f"Could not load or build a minimap2 index from {source}")
        if (aligner.k, aligner.w) != SHORT_READ_INDEX_K_W:
            print(
                f"warning: {source} was indexed with k={aligner.k}, w={aligner.w}; "
                "short-read indexes (minimap2 -x sr) use k=21, w=11",
                file=sys.stderr,
            )
        return cls(aligner, mappy)

    def new_buffer(self):
        return self._mappy.ThreadBuffer()

    def hits(self, sequence: str, buffer) -> List[GenomeHit]:
        return [_genome_hit(hit) for hit in self._aligner.map(sequence, buf=buffer)]

    def fetch(self, contig: str, start: int, end: int) -> str:
        """The genome sequence of an interval; empty when the index holds no sequence there."""
        return self._aligner.seq(contig, start, end) or ""

    def window(self, contig: str, start: int, end: int) -> "WindowIndex":
        return WindowIndex(self._mappy, contig, start, self.fetch(contig, start, end))


class WindowIndex:
    """An index of one genome window with short seeds (``SITE_INDEX_K``,
    ``SITE_INDEX_W``), so that a tag too diverged to seed in the genome-wide
    search still aligns at a site it belongs to.  Each window is an index of
    its own, so windows never compete for the aligner's limit on hits.  The
    mapping options are the short-read preset's, not
    ``GENOME_WIDE_MAP_OPTIONS``.  Hits come back in genome coordinates."""

    def __init__(self, mappy_module, contig: str, start: int, sequence: str) -> None:
        self._contig, self._start = contig, start
        self._aligner = (
            mappy_module.Aligner(
                seq=sequence, preset="sr", k=SITE_INDEX_K, w=SITE_INDEX_W, scoring=list(SHORT_READ_SCORING)
            )
            if sequence
            else None
        )

    def hits(self, sequence: str, buffer=None) -> List[GenomeHit]:
        if self._aligner is None:
            return []
        return [_genome_hit(hit, self._contig, self._start) for hit in self._aligner.map(sequence, buf=buffer)]


def map_sequences(
    aligner, sequences: Sequence[str], jobs: int = 1, max_sites: Optional[int] = MAX_SITES_PER_HAPLOTYPE
) -> Dict[str, GenomePlacement]:
    """Every sequence's distinct genome sites (``distinct_sites``).

    Worker threads share the index; mappy releases Python's interpreter lock
    while it maps, and each batch has a buffer of its own.
    """

    def align_batch(batch: Sequence[str]) -> List[GenomePlacement]:
        buffer = aligner.new_buffer()
        return [distinct_sites(aligner.hits(sequence, buffer), max_sites) for sequence in batch]

    batches = [sequences[start : start + MAPPING_BATCH] for start in range(0, len(sequences), MAPPING_BATCH)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        placements = [placement for found in pool.map(align_batch, batches) for placement in found]
    return dict(zip(sequences, placements))


def clone_position(clone_id: str, contigs: Set[str]) -> Tuple[Optional[Tuple[str, int]], Optional[str]]:
    """The (contig, 1-based position) a CHROM_POS CloneID names, with None; or None, with why it names none.

    Any CloneID that ends in an underscore and digits is read as CHROM_POS
    when what precedes them is a contig of the genome.
    """
    match = CLONE_POSITION_RE.fullmatch(clone_id)
    if not match:
        return None, "no_position"
    if match.group(1) not in contigs:
        return None, "contig_absent"
    return (match.group(1), int(match.group(2))), None


def _site_is_at(site: GenomeHit, contig: str, position: int) -> bool:
    """Whether a 1-based position lies in a site (0-based, half-open) or within ``CLONE_POSITION_SLACK`` of it."""
    return site.contig == contig and site.start - CLONE_POSITION_SLACK < position <= site.end + CLONE_POSITION_SLACK


def place_marker(
    clone_id: str,
    ref_placement: Optional[GenomePlacement],
    contigs: Set[str],
    tag_length: int,
    use_clone_position: bool = True,
    align_ref_tag: Optional[Callable[[str, int, int], Sequence[GenomeHit]]] = None,
) -> MarkerPlacement:
    """Locate a marker's tag from its Ref tag's equal-best genome sites and its CloneID position.

    ``ref_placement`` is the Ref tag's genome sites, or None for a marker
    without a Ref row.  ``tag_length`` sizes the window around a CloneID
    position.

    With ``use_clone_position`` the CloneID position is taken to be on this
    genome.  It picks one of several equal-best sites, and a window around it
    stands in for a Ref tag that does not align.  If two or more equal-best
    sites are at the position (tandem copies), none is picked and every
    equal-best site stays the marker's.  When no equal-best site is at the
    position, ``align_ref_tag(contig, start, end)`` aligns the Ref tag within
    a window around it, and an alignment at the position becomes the marker's
    site even if the tag aligns better elsewhere: the assembly may carry the
    Alt allele at the target and a Ref-like copy elsewhere, or the
    genome-wide search may have missed the site.  With no such alignment the
    marker is unplaced.  Either way the disagreement is reported.  Without
    ``use_clone_position`` (a genome other than the one the CloneIDs refer
    to) the position is never read: an unaligned Ref tag leaves the marker
    unplaced, and equal-best sites stay multi-copy.
    """
    position, no_position_reason = clone_position(clone_id, contigs) if use_clone_position else (None, "ignored")
    equal_best = ref_placement.equal_best if ref_placement is not None else ()
    if position is None:
        if not equal_best:
            return MarkerPlacement((), "unplaced", 0, no_position_reason)
        placed_by = "ref_tag" if len(equal_best) == 1 else "ref_tag_multicopy"
        return MarkerPlacement(equal_best, placed_by, len(equal_best), no_position_reason)
    contig, clone_base = position
    # The genome interval a tag of this length lies in if it contains the CloneID position.
    window_start, window_end = max(0, clone_base - 1 - tag_length), clone_base + tag_length
    if not equal_best:
        return MarkerPlacement((GenomeHit(contig, window_start, window_end, 0, 0),), "clone_id", 0, "fallback")
    at_position = [site for site in equal_best if _site_is_at(site, contig, clone_base)]
    if len(at_position) == 1:
        placed_by = "ref_tag" if len(equal_best) == 1 else "ref_tag+clone_id"
        return MarkerPlacement(tuple(at_position), placed_by, len(equal_best), "agree")
    if at_position:
        return MarkerPlacement(equal_best, "ref_tag_multicopy", len(equal_best), "agree")
    realigned = ()
    if align_ref_tag is not None:
        realigned = align_ref_tag(contig, max(0, window_start - SITE_FLANK), window_end + SITE_FLANK)
    realigned_at_position = [hit for hit in realigned if _site_is_at(hit, contig, clone_base)]
    if realigned_at_position:
        best = max(realigned_at_position, key=lambda hit: hit.score)
        return MarkerPlacement((best,), "clone_id_site", len(equal_best), "disagree")
    return MarkerPlacement((), "unplaced", len(equal_best), "disagree")


def score_at_sites(
    aligner, markers: Mapping[str, MarkerPlacement], haplotypes: Iterable[Tuple[str, str]], jobs: int = 1
) -> Dict[Tuple[str, str], int]:
    """Each (CloneID, sequence) pair's best score within its own marker's sites.

    A site is searched with ``SITE_FLANK`` genome bases either side of it.
    Every placed marker's sites are indexed on their own with short seeds
    (``WindowIndex``), so an allele of the marker that carries several SNPs,
    which the genome-wide search misses or aligns only in part at the
    marker's site, is still scored there.  Markers that share a site do not
    compete for the aligner's limit on hits.  Pairs that do not align at
    their marker's sites, and pairs of unplaced markers, are left out.
    """
    sequences_by_marker: Dict[str, Set[str]] = defaultdict(set)
    for clone_id, sequence in haplotypes:
        if markers[clone_id].sites:
            sequences_by_marker[clone_id].add(sequence)

    def score_marker(clone_id: str) -> Dict[Tuple[str, str], int]:
        windows = [
            aligner.window(site.contig, max(0, site.start - SITE_FLANK), site.end + SITE_FLANK)
            for site in markers[clone_id].sites
        ]
        buffer = aligner.new_buffer()
        best_scores = {}
        for sequence in sequences_by_marker[clone_id]:
            scores = [hit.score for window in windows for hit in window.hits(sequence, buffer)]
            if scores:
                best_scores[(clone_id, sequence)] = max(scores)
        return best_scores

    site_scores: Dict[Tuple[str, str], int] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for best_scores in pool.map(score_marker, list(sequences_by_marker)):
            site_scores.update(best_scores)
    return site_scores


def genome_label(placement: GenomePlacement, marker: MarkerPlacement, site_score: Optional[int] = None) -> GenomeLabel:
    """Label a haplotype by its best score at its marker's own site against its best score anywhere else.

    ``on_target``: the marker's site scores strictly best.  ``paralog``:
    somewhere else scores strictly better (or the haplotype aligns only
    elsewhere).  ``ambiguous``: a tie between the marker's site and another
    site, or a marker whose own Ref tag is multi-copy and unresolved.
    ``unmapped``: no alignment.  ``unplaced``: the marker itself could not be
    located, so the haplotype cannot be judged.

    ``placement`` is the haplotype's sites from the genome-wide search and
    ``marker`` its marker's placement.  ``site_score`` is the haplotype's
    score in the realignment at the marker's site (``score_at_sites``), or
    None when it did not align there.  The target score is the higher of it
    and the genome-wide hit at the site, so the label does not depend on
    which search found the site.
    """
    hits_at_marker = [hit for hit in placement.sites if any(_same_locus(hit, site) for site in marker.sites)]
    hits_elsewhere = [hit for hit in placement.sites if hit not in hits_at_marker]
    off_target_score = hits_elsewhere[0].score if hits_elsewhere else None
    off_target_sites = sum(hit.score == off_target_score for hit in hits_elsewhere)
    off_target_hit = hits_elsewhere[0] if hits_elsewhere else None
    if not marker.sites:
        return GenomeLabel("unplaced", None, ".", off_target_score, off_target_sites, off_target_hit)
    if site_score is not None and (not hits_at_marker or site_score > hits_at_marker[0].score):
        target_score, target_from = site_score, "site"
    elif hits_at_marker:
        target_score, target_from = hits_at_marker[0].score, "genome"
    else:
        target_score, target_from = None, "."
    if target_score is None and off_target_score is None:
        label = "unmapped"
    elif target_score is None or (off_target_score is not None and off_target_score > target_score):
        label = "paralog"
    elif off_target_score == target_score or len(marker.sites) > 1:
        label = "ambiguous"
    else:
        label = "on_target"
    return GenomeLabel(label, target_score, target_from, off_target_score, off_target_sites, off_target_hit)


def _preamble_header_and_rows(handle, path: Path):
    """(preamble rows, header row, reader positioned at the first data row) of an open MADC."""
    reader = csv.reader(handle)
    preamble = []
    for row in reader:
        if is_madc_header(row):
            return preamble, row, reader
        preamble.append(row)
    raise MadcError(f"Could not find an AlleleID header row in MADC file: {path}")


def _haplotype_rows(reader, header: Sequence[str], path: Path):
    """The haplotype rows of a MADC, each as (row, AlleleID, CloneID, sequence in upper case).

    A row with an empty AlleleID, or too short to hold the three identity
    columns, is not a haplotype row and is passed over.
    """
    allele_index, clone_index, sequence_index = madc_identity_columns(header, path)
    width = max(allele_index, clone_index, sequence_index) + 1
    for row in reader:
        if len(row) < width or not row[allele_index].strip():
            continue
        yield row, row[allele_index].strip(), row[clone_index].strip(), row[sequence_index].strip().upper()


def _label_column_index(header: Sequence[str], madc_path: Path) -> int:
    """The column that carries the flag: a GenomeLabel column from an earlier run, or PICRef."""
    stripped = [column.strip() for column in header]
    for name in (GENOME_LABEL_COLUMN, REPLACED_STATISTICS_COLUMN):
        if name in stripped:
            return stripped.index(name)
    raise MadcError(
        f"{madc_path} has no {REPLACED_STATISTICS_COLUMN} column to carry the flag; use drop mode, "
        "or read the labels sidecar"
    )


def check_madc_header(madc_path: Path, flagging: bool) -> None:
    """Stop with an error when the MADC cannot be filtered as asked.

    That is a MADC without a header row or without its identity columns, or,
    when rows are to be flagged, without a column to carry the flag.
    """
    with open_text(madc_path, "rt") as handle:
        _, header, _ = _preamble_header_and_rows(handle, madc_path)
    madc_identity_columns(header, madc_path)
    if flagging:
        _label_column_index(header, madc_path)


def _read_haplotypes(madc_path: Path):
    """Pass 1 over the MADC: its distinct haplotypes, and each marker's Ref tag and tag length.

    Returns the header row and three dictionaries, each in the order of first
    appearance.  ``haplotypes`` maps (CloneID, sequence as written) to the
    sequence that is aligned, which is the sequence without gap characters.
    ``ref_tags`` maps a CloneID to its first Ref row's sequence without gaps,
    and ``tag_lengths`` to the length of its first row as written.
    """
    haplotypes: Dict[Tuple[str, str], str] = {}
    ref_tags: Dict[str, str] = {}
    tag_lengths: Dict[str, int] = {}
    with open_text(madc_path, "rt") as handle:
        _, header, reader = _preamble_header_and_rows(handle, madc_path)
        for _, allele_id, clone_id, sequence in _haplotype_rows(reader, header, madc_path):
            ungapped = sequence.replace("-", "")
            haplotypes[(clone_id, sequence)] = ungapped
            tag_lengths.setdefault(clone_id, len(sequence))
            if clone_id not in ref_tags and allele_type(allele_id) == "Ref":
                ref_tags[clone_id] = ungapped
    return header, haplotypes, ref_tags, tag_lengths


def _place_markers(
    aligner, ref_tags: Mapping[str, str], tag_lengths: Mapping[str, int], use_clone_position: bool, jobs: int
) -> Dict[str, MarkerPlacement]:
    """Every marker's placement (``place_marker``), from the genome sites of its Ref tag."""
    # Ref tags keep every equal-best site, so a capped list cannot hide the one the CloneID names.
    ref_placements = map_sequences(aligner, sorted(set(ref_tags.values())), jobs, max_sites=None)
    markers: Dict[str, MarkerPlacement] = {}
    for clone_id, tag_length in tag_lengths.items():
        ref_tag = ref_tags.get(clone_id)   # None for a marker without a Ref row, which place_marker never aligns

        def align_ref_tag(contig: str, start: int, end: int, ref_tag=ref_tag) -> List[GenomeHit]:
            return aligner.window(contig, start, end).hits(ref_tag)

        markers[clone_id] = place_marker(
            clone_id, ref_placements.get(ref_tag), aligner.contigs, tag_length, use_clone_position, align_ref_tag
        )
    return markers


def _format_site(site: GenomeHit) -> str:
    """A site as contig:start-end(strand), 1-based and inclusive, as the two sidecars show it."""
    strand = {1: "+", -1: "-"}.get(site.strand, ".")
    return f"{site.contig}:{site.start + 1}-{site.end}({strand})"


def _dot_if_none(value) -> object:
    return "." if value is None else value


def _write_markers_sidecar(markers_path: Path, markers: Mapping[str, MarkerPlacement]) -> None:
    markers_path.parent.mkdir(parents=True, exist_ok=True)
    with markers_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(MARKERS_SIDECAR_COLUMNS)
        for clone_id, marker in markers.items():
            sites = ";".join(_format_site(site) for site in marker.sites) or "."
            writer.writerow([clone_id, marker.placed_by, marker.ref_tag_site_count, marker.clone_id_check, sites])


class _RowTallies(NamedTuple):
    """Rows and reads counted while the table is written, in total and by label, outcome and allele category."""

    totals: Dict[str, Dict[str, int]]
    by_label: Dict[str, Dict[str, int]]
    by_outcome: Dict[str, Dict[str, int]]
    by_category: Dict[str, Dict[str, Dict[str, int]]]
    overflow: CountOverflowCheck


def _write_filtered_table(
    madc_path: Path,
    output: Path,
    labels_path: Path,
    labels: Mapping[Tuple[str, str], GenomeLabel],
    actions: Mapping[str, str],
    label_index: Optional[int],
    provenance: str,
) -> _RowTallies:
    """Pass 2 over the MADC: write the rows that are kept, and every row's label to the labels sidecar.

    ``actions`` gives what was asked for ``paralog`` and ``ambiguous`` rows
    (``drop``, ``flag`` or ``keep``); rows with any other label are kept.
    What is done with a row is its outcome, the sidecar's Action column:
    ``kept``, ``flagged``, ``dropped``, or ``protected`` for a canonical Ref
    or Alt row that would have been dropped and is kept.  With ``label_index``
    every row written carries its label in that column, whatever its outcome.
    """

    def tally():
        return defaultdict(lambda: {"rows": 0, "reads": 0})

    by_label, by_outcome = tally(), tally()
    by_category: Dict[str, Dict[str, Dict[str, int]]] = defaultdict(tally)
    overflow = CountOverflowCheck()
    totals = {"input": {"rows": 0, "reads": 0}, "output": {"rows": 0, "reads": 0}}
    output.parent.mkdir(parents=True, exist_ok=True)
    labels_path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(madc_path, "rt") as handle, open_text(output, "wt") as out_handle, labels_path.open(
        "w", newline=""
    ) as label_handle:
        preamble, header, reader = _preamble_header_and_rows(handle, madc_path)
        writer = csv.writer(out_handle, lineterminator="\n")
        label_writer = csv.writer(label_handle, lineterminator="\n")
        label_writer.writerow(LABELS_SIDECAR_COLUMNS)
        if provenance and preamble and preamble[0] and preamble[0][0].startswith("* MADCATS"):
            preamble[0] = [f"{preamble[0][0]} | {provenance}"] + preamble[0][1:]
        writer.writerows(preamble)
        if label_index is not None:
            header[label_index] = GENOME_LABEL_COLUMN
        writer.writerow(header)
        sample_start = madc_sample_start(header)
        for row, allele_id, clone_id, sequence in _haplotype_rows(reader, header, madc_path):
            haplotype_label = labels[(clone_id, sequence)]
            category = allele_type(allele_id)
            requested = actions.get(haplotype_label.label, "keep")
            if requested == "drop":
                outcome = "protected" if category in ("Ref", "Alt") else "dropped"
            else:
                outcome = "flagged" if requested == "flag" else "kept"
            counts = madc_read_counts(row[sample_start:], madc_path)
            overflow.add(clone_id, counts)
            reads = sum(counts)
            label = haplotype_label.label
            for counted in (totals["input"], by_label[label], by_outcome[outcome], by_category[category][label]):
                counted["rows"] += 1
                counted["reads"] += reads
            if outcome != "dropped":
                if label_index is not None:
                    row = row + [""] * max(0, label_index + 1 - len(row))
                    row[label_index] = label
                writer.writerow(row)
                totals["output"]["rows"] += 1
                totals["output"]["reads"] += reads
            off_target_hit = haplotype_label.off_target_hit
            label_writer.writerow([
                clone_id, allele_id, sequence, label, outcome,
                _dot_if_none(haplotype_label.target_score), haplotype_label.target_from,
                _dot_if_none(haplotype_label.off_target_score), haplotype_label.off_target_sites,
                "." if off_target_hit is None else _format_site(off_target_hit),
            ])
    return _RowTallies(totals, by_label, by_outcome, by_category, overflow)


def filter_paralogs(
    madc_path: Path,
    output: Path,
    aligner,
    labels_path: Path,
    markers_path: Path,
    paralog_action: str = "drop",
    ambiguous_action: str = "drop",
    jobs: int = 1,
    provenance: str = "",
    clone_positions: str = "use",
) -> Dict[str, object]:
    """Label every MADC haplotype by where the genome places it; drop or flag paralogous rows.

    Pass 1 reads only row identities (``_read_haplotypes``).  Every distinct
    haplotype is aligned to the genome once; each marker is located from its
    Ref row's tag, all of whose equal-best sites are kept (``place_marker``);
    every haplotype is realigned to its own marker's sites
    (``score_at_sites``) and labelled by ``genome_label``.  Pass 2 streams the
    table again and writes the kept rows unchanged
    (``_write_filtered_table``), so a vendor MADC keeps its preamble,
    statistics and column layout, and a HapApp fixed-allele-ID MADC keeps its
    Code_version column.  Canonical Ref and Alt rows are never dropped, since
    downstream dosage tools need both designed alleles; they are reported as
    ``protected`` instead.  In flag mode the PICRef column becomes
    ``GenomeLabel`` on every row.  The labels sidecar has one line per
    haplotype row, dropped rows included.  ``clone_positions="use"`` takes
    CloneID coordinates to be on this genome; ``"ignore"`` keeps them out of
    marker placement, for a genome other than the one they refer to.

    Returns the summary that goes into the QC JSON.
    """
    actions = {"paralog": paralog_action, "ambiguous": ambiguous_action}
    header, haplotypes, ref_tags, tag_lengths = _read_haplotypes(madc_path)
    label_index = _label_column_index(header, madc_path) if "flag" in actions.values() else None

    sequences = sorted(set(haplotypes.values()))
    print(f"Aligning {len(sequences):,} distinct haplotype sequences on {jobs} thread(s)", file=sys.stderr)
    started = time.time()
    placements = map_sequences(aligner, sequences, jobs)
    markers = _place_markers(aligner, ref_tags, tag_lengths, clone_positions == "use", jobs)
    # Realign every haplotype at its own marker's sites, and label it by the better of the two scores there.
    aligned_pairs = {(clone_id, ungapped) for (clone_id, _), ungapped in haplotypes.items()}
    site_scores = score_at_sites(aligner, markers, aligned_pairs, jobs)
    labels = {
        (clone_id, sequence): genome_label(
            placements[ungapped], markers[clone_id], site_scores.get((clone_id, ungapped))
        )
        for (clone_id, sequence), ungapped in haplotypes.items()
    }
    alignment_seconds = round(time.time() - started, 1)

    _write_markers_sidecar(markers_path, markers)
    try:
        tallies = _write_filtered_table(madc_path, output, labels_path, labels, actions, label_index, provenance)
    except BaseException:
        # A run that fails part of the way through leaves no partial files behind.
        for partial in (output, labels_path, markers_path):
            partial.unlink(missing_ok=True)
        raise

    warning = tallies.overflow.warning(madc_path)
    if warning:
        print(warning, file=sys.stderr)
    sources = Counter(marker.placed_by for marker in markers.values())
    checks = Counter(marker.clone_id_check for marker in markers.values())
    checked = checks["agree"] + checks["disagree"]
    if checked and checks["disagree"] > CLONE_DISAGREEMENT_WARNING_SHARE * checked:
        print(
            f"warning: CloneID positions disagree with the Ref tag's best genome sites for {checks['disagree']:,} "
            f"of {checked:,} markers, so they may refer to another assembly. Those markers were placed at their "
            "CloneID position, or left unplaced where the Ref tag does not align there. If the positions refer "
            "to another assembly, rerun with --clone-positions ignore",
            file=sys.stderr,
        )
    if markers and not sources.keys() - {"unplaced"}:
        print("warning: no marker could be placed on the genome; every row was kept", file=sys.stderr)

    def plain(tally: Mapping[str, Dict[str, int]]) -> Dict[str, Dict[str, int]]:
        return {key: dict(counts) for key, counts in tally.items()}

    return {
        "markers": {
            "total": len(markers),
            "placement": dict(sources),
            "clone_id_check": dict(checks),
            "multicopy_ref_tags": [clone_id for clone_id, marker in markers.items() if marker.ref_tag_site_count > 1],
            "clone_id_disagreements": [
                clone_id for clone_id, marker in markers.items() if marker.clone_id_check == "disagree"
            ],
            "unplaced": [clone_id for clone_id, marker in markers.items() if marker.placed_by == "unplaced"],
        },
        "haplotypes": {
            "distinct": len(haplotypes),
            "distinct_sequences": len(sequences),
            "labels": dict(Counter(result.label for result in labels.values())),
            "realigned_at_site": {
                "haplotypes": sum(bool(markers[clone_id].sites) for clone_id, _ in haplotypes),
                "found_at_site": sum(pair in site_scores for pair in
                                     ((clone_id, ungapped) for (clone_id, _), ungapped in haplotypes.items())),
                "target_from_site": sum(result.target_from == "site" for result in labels.values()),
            },
            "alignment_seconds": alignment_seconds,
        },
        "rows": {
            **tallies.totals,
            "by_label": plain(tallies.by_label),
            "by_outcome": plain(tallies.by_outcome),
            "by_category": {category: plain(by_label) for category, by_label in tallies.by_category.items()},
        },
        "genome_label_column": GENOME_LABEL_COLUMN if label_index is not None else None,
        "count_overflow": tallies.overflow.summary(),
    }


def run_filter_paralogs(args: argparse.Namespace) -> int:
    """Clean a finished MADC with a reference genome (see ``filter_paralogs``)."""
    if not args.genome and not args.genome_index:
        raise MadcError("filter-paralogs needs --genome, --genome-index, or both")
    madc = Path(args.madc)
    output = Path(args.output)
    if not madc.is_file():
        raise MadcError(f"MADC not found: {madc}")
    if output.resolve() == madc.resolve():
        raise MadcError("--output must be a new file, not the input MADC")
    genome = Path(args.genome) if args.genome else None
    index = Path(args.genome_index) if args.genome_index else None
    # A MADC that cannot be filtered as asked is refused before the index is loaded, which can take minutes.
    check_madc_header(madc, flagging="flag" in (args.paralogs, args.ambiguous))
    started = time.time()
    aligner = GenomeAligner.open(genome, index, args.jobs)
    labels_path = output.with_suffix(output.suffix + ".genome_labels.csv")
    markers_path = output.with_suffix(output.suffix + ".genome_markers.csv")
    provenance = (
        f"filter-paralogs v{__version__} paralogs={args.paralogs} ambiguous={args.ambiguous} "
        f"clone_positions={args.clone_positions}"
    )
    summary = filter_paralogs(
        madc, output, aligner, labels_path, markers_path,
        paralog_action=args.paralogs, ambiguous_action=args.ambiguous, jobs=args.jobs, provenance=provenance,
        clone_positions=args.clone_positions,
    )
    qc_path = Path(args.qc_output) if args.qc_output else output.with_suffix(output.suffix + ".qc.json")
    qc = {
        "pipeline": "MADCATS",
        "command": "filter-paralogs",
        "version": __version__,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "inputs": {
            "madc": str(madc.resolve()),
            "genome": str(genome.resolve()) if genome else None,
            "genome_index": str(index.resolve()) if index else None,
        },
        "parameters": {
            "paralogs": args.paralogs, "ambiguous": args.ambiguous, "clone_positions": args.clone_positions,
            "jobs": args.jobs,
        },
        "aligner": aligner.description,
        "outputs": {"madc": str(output.resolve()), "labels": str(labels_path), "markers": str(markers_path)},
        "elapsed_seconds": round(time.time() - started, 1),
        **summary,
    }
    _write_json(qc_path, qc)
    rows = summary["rows"]
    dropped_reads = rows["input"]["reads"] - rows["output"]["reads"]
    dropped_share = dropped_reads / rows["input"]["reads"] if rows["input"]["reads"] else 0.0
    print(
        f"Kept {rows['output']['rows']:,} of {rows['input']['rows']:,} rows; "
        f"dropped {dropped_reads:,} reads ({dropped_share:.1%})",
        file=sys.stderr,
    )
    print(
        f"Wrote {output}, its labels {labels_path}, its markers {markers_path}, and QC to {qc_path}", file=sys.stderr
    )
    return 0


def add_filter_paralogs_command(subparsers) -> None:
    filtering = subparsers.add_parser(
        "filter-paralogs",
        help="drop or flag MADC rows that a reference genome places at another locus (needs mappy)",
        description=(
            "Drop or flag the rows of a MADC that a reference genome places at another locus, for panels where "
            "paralogous copies of a target co-amplify with it. Each marker is placed on the genome by its Ref "
            "tag and, for CHROM_POS CloneIDs, by that position. Every row's sequence is then aligned and "
            "labelled on_target, paralog (it aligns better elsewhere) or ambiguous (it aligns equally well at "
            "the marker and elsewhere). Needs the mappy aligner: pip install 'madcats[genome]'."
        ),
    )
    filtering.add_argument(
        "--madc",
        required=True,
        help=(
            "MADC to clean, plain or gzipped: one written by generate, a vendor DArT MADC (with or without its "
            "starred preamble) or a HapApp fixed-allele-ID MADC"
        ),
    )
    filtering.add_argument(
        "--genome",
        help=(
            "reference genome FASTA, plain or gzipped. It should hold each locus once: of a haplotype-resolved "
            "assembly, give one set of chromosomes. Not needed when --genome-index exists"
        ),
    )
    filtering.add_argument(
        "--genome-index",
        help=(
            "minimap2 index (.mmi) built with -x sr, in one part; loaded when it exists, otherwise built from "
            "--genome and saved here for the next run"
        ),
    )
    filtering.add_argument(
        "--clone-positions",
        choices=("use", "ignore"),
        default="use",
        help=(
            "CloneIDs of the form CHROM_POS: use takes them as positions on this genome and places each marker "
            "there, leaving it unplaced when its Ref tag aligns only elsewhere; ignore places markers by "
            "their Ref tag alone, for CloneIDs that refer to another assembly (default: use)"
        ),
    )
    filtering.add_argument(
        "--output",
        required=True,
        help=(
            "cleaned MADC to write. Written beside it: OUTPUT.genome_labels.csv (every input row's label and "
            "scores, and what was done with it) and OUTPUT.genome_markers.csv (where each marker was placed)"
        ),
    )
    filtering.add_argument(
        "--paralogs",
        choices=("drop", "flag"),
        default="drop",
        help=(
            "paralog rows: drop them, or flag them (keep them, and write every row's label in a GenomeLabel "
            "column that takes the place of PICRef). Ref and Alt rows are never dropped, nor are rows that do "
            "not align or whose marker could not be placed (default: drop)"
        ),
    )
    filtering.add_argument(
        "--ambiguous",
        choices=("drop", "flag", "keep"),
        default="drop",
        help=(
            "ambiguous rows, which include the rows of a marker whose Ref tag has several equal-best sites: "
            "drop them, flag them as --paralogs flag does, or keep them (default: drop)"
        ),
    )
    filtering.add_argument("--qc-output", help="QC JSON (default: OUTPUT.qc.json)")
    filtering.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(), help="alignment threads (default: up to 4)"
    )
    filtering.set_defaults(func=run_filter_paralogs)
