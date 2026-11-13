"""``madcats filter-paralogs``: clean a finished MADC with a reference genome.

In panels where paralogous copies co-amplify with a target, a paralog read
that differs from the Ref tag only at non-target bases becomes a RefMatch or
AltMatch row of that marker.  This module places each marker on the genome
by its Ref tag, aligns every haplotype, and drops or flags the rows whose
best alignment lies at another locus (paralogs) or ties between loci
(ambiguous).

It is kept apart from ``madcats.py`` because it is the only command with a
dependency outside the standard library: the ``mappy`` aligner (minimap2),
imported on first use (``pip install 'madcats[genome]'``), so ``madcats``
runs without it.  The MADC reading helpers it shares with ``validate``
(header detection, sample columns, allele categories, the count-overflow
check) stay in ``madcats``.

Reading map: scoring and alignment (``hit_score``, ``GenomeAligner``,
``WindowIndex``, ``map_sequences``); marker placement and labels
(``place_marker``, ``score_at_sites``, ``genome_label``); the two-pass
table filter (``filter_paralogs``); the command line
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

GENOME_LABEL_COLUMN = "PICRef"     # the vendor statistics column the GenomeLabel flag replaces
GENOME_LABEL_HEADER = "GenomeLabel"
# minimap2 -x sr scoring: match, mismatch, gap open/extend, long-gap open/extend, ambiguous base.
SR_SCORING = (2, 8, 12, 2, 24, 1, 1)
# Map-time sensitivity on top of -x sr (the index is unchanged): single-seed chains and up to 50
# hits, so an 81-bp allele with several SNPs still finds its copies.  On lettuce this lifted
# drop/keep agreement with bwa mem from 97.3% to 98.2%.
SR_MAP_OPTIONS = {"min_cnt": 1, "min_chain_score": 20, "best_n": 50}
GENOME_SITES_KEPT = 8              # distinct genome sites remembered per haplotype, best first (Ref tags keep all)
CLONE_POSITION_SLACK = 10          # a CloneID position this far outside the placed Ref tag still agrees
CLONE_POSITION_RE = re.compile(r"(.+)_0*(\d+)")
SITE_PAD = 20                      # genome bases either side of a site in the window realignment
MAPPING_BATCH = 2000               # sequences handed to a worker thread at a time
SITE_INDEX_K, SITE_INDEX_W = 11, 3 # the site realignment's seeds: short, so diverged alleles still seed
CLONE_DISAGREEMENT_WARNING_SHARE = 0.1   # warn when more than this share of CloneID positions disagree with the Ref tags
GENOME_LABEL_SIDECAR_COLUMNS = [
    "CloneID", "AlleleID", "AlleleSequence", "Label", "Action",
    "TargetScore", "TargetFrom", "OffTargetScore", "OffTargetSites", "OffTargetHit",
]
GENOME_MARKER_SIDECAR_COLUMNS = ["CloneID", "Placement", "RefTagSites", "CloneIDCheck", "Sites"]


class GenomeHit(NamedTuple):
    """One alignment to the genome: a 0-based, half-open reference interval and its score."""

    contig: str
    start: int
    end: int
    strand: int
    score: int


class GenomePlacement(NamedTuple):
    """A sequence's distinct genome sites, best score first (hits overlapping a better one are merged into it)."""

    sites: Tuple[GenomeHit, ...]

    @property
    def best(self) -> Tuple[GenomeHit, ...]:
        """The equal-best sites."""
        return tuple(site for site in self.sites if site.score == self.sites[0].score)


class MarkerPlacement(NamedTuple):
    """Where a marker's own tag lies in the genome.

    ``source``: ``ref_tag`` (the Ref tag has one best site), ``ref_tag+clone_id``
    (several equal-best sites; the CloneID position picks one),
    ``ref_tag_multicopy`` (several, none picked: every site is the marker's),
    ``clone_id_site`` (no equal-best site covers the CloneID position; the Ref
    tag's alignment within a window around it), ``clone_id`` (the Ref tag did
    not align; the window around the CloneID position stands in) or
    ``unplaced``.  ``clone_check`` compares the CloneID position with the
    Ref-tag sites: ``agree``, ``disagree``, ``no_position`` (the CloneID is not
    CHROM_POS), ``contig_absent``, ``fallback``, or ``ignored``
    (``--clone-positions ignore``).
    """

    sites: Tuple[GenomeHit, ...]
    source: str
    ref_tag_sites: int
    clone_check: str


class GenomeLabel(NamedTuple):
    """A haplotype's label and the two scores that decided it."""

    label: str
    target_score: Optional[int]       # best alignment at the marker's own site
    target_from: str                  # genome (the genome-wide hit there), site (the realignment there scored higher) or '.'
    off_target_score: Optional[int]   # best alignment anywhere else
    off_target_sites: int             # distinct sites elsewhere at that score
    off_target_hit: Optional[GenomeHit]


def hit_score(cigar: Sequence[Sequence[int]], edit_distance: int, matches: int, block_length: int) -> int:
    """An alignment's DP score recomputed from what mappy reports (it gives neither AS nor XS).

    minimap2's block length counts aligned columns and gap bases but not
    ambiguous columns, and NM counts all three kinds of difference, so
    mismatches = block - matches - gap bases and ambiguous = NM - (block -
    matches).  Gaps cost the cheaper of minimap2's two affine pieces.
    """
    match, mismatch, gap_open, gap_extend, long_open, long_extend, ambiguous = SR_SCORING
    gap_bases = gap_cost = 0
    for length, operation in cigar:
        if operation in (1, 2):  # I, D
            gap_bases += length
            gap_cost += min(gap_open + gap_extend * length, long_open + long_extend * length)
    differences = block_length - matches
    mismatches = max(0, differences - gap_bases)
    ambiguous_bases = max(0, edit_distance - differences)
    return match * matches - mismatch * mismatches - ambiguous * ambiguous_bases - gap_cost


def _sites_overlap(hit: GenomeHit, site: GenomeHit) -> bool:
    """Same locus: the intervals share at least half of the shorter one.  Adjacent
    tandem copies, which share a few bases, stay distinct."""
    if hit.contig != site.contig:
        return False
    overlap = min(hit.end, site.end) - max(hit.start, site.start)
    return overlap > 0 and 2 * overlap >= min(hit.end - hit.start, site.end - site.start)


def summarize_hits(hits: Iterable[GenomeHit], keep: Optional[int] = GENOME_SITES_KEPT) -> GenomePlacement:
    """A sequence's distinct sites, best first, at most ``keep`` (None: all).  A hit at
    the same locus as a better-scoring one (``_sites_overlap``) is merged into it."""
    sites: List[GenomeHit] = []
    for hit in sorted(hits, key=lambda hit: (-hit.score, hit.contig, hit.start)):
        if any(_sites_overlap(hit, site) for site in sites):
            continue
        sites.append(hit)
        if len(sites) == keep:
            break
    return GenomePlacement(tuple(sites))


def _index_parts(path: Path) -> Optional[List[List[str]]]:
    """The contig names in each part of a minimap2 index, read from its layout
    (``mm_idx_dump``), or None when the file is not one.  minimap2 writes one
    part per -I bases of genome; mappy loads only the first."""
    size = path.stat().st_size
    with path.open("rb", buffering=0) as handle:

        def read(count: int) -> bytes:
            data = handle.read(count)
            if len(data) != count:
                raise MadcError(f"{path} ends inside a minimap2 index part; rebuild the index")
            return data

        parts: List[List[str]] = []
        while handle.tell() < size:
            if read(4) != b"MMI\x02":
                if not parts:
                    return None
                raise MadcError(f"{path} has unreadable data after its first index part; rebuild the index")
            _, _, bucket_bits, sequences, flag = struct.unpack("<5I", read(20))
            names, bases = [], 0
            for _ in range(sequences):
                names.append(read(read(1)[0]).decode())
                bases += struct.unpack("<I", read(4))[0]
            for _ in range(1 << bucket_bits):
                handle.seek(8 * struct.unpack("<I", read(4))[0], 1)
                handle.seek(16 * struct.unpack("<I", read(4))[0], 1)
            if not flag & 2:  # MM_I_NO_SEQ: the packed sequence is left out
                handle.seek(4 * ((bases + 7) // 8), 1)
            parts.append(names)
        if handle.tell() != size:
            raise MadcError(f"{path} ends inside a minimap2 index part; rebuild the index")
        return parts


def _genome_hit(hit, contig: Optional[str] = None, offset: int = 0) -> GenomeHit:
    """A mappy alignment as a GenomeHit; ``contig`` and ``offset`` place a hit within a window on the genome."""
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
    """minimap2's short-read preset, through mappy, with the scoring pinned to SR_SCORING."""

    def __init__(self, aligner, mappy_module) -> None:
        self._aligner = aligner
        self._mappy = mappy_module
        self.contigs = frozenset(aligner.seq_names)
        self.description = {
            "mappy": getattr(mappy_module, "__version__", None),
            "preset": "sr",
            "k": aligner.k,
            "w": aligner.w,
            "scoring": list(SR_SCORING),
            "map_options": dict(SR_MAP_OPTIONS),
            "site_realignment": {"k": SITE_INDEX_K, "w": SITE_INDEX_W, "pad": SITE_PAD},
        }

    @classmethod
    def open(cls, genome: Optional[Path], index: Optional[Path], threads: int) -> "GenomeAligner":
        """Load ``index`` when it exists; otherwise index ``genome`` (saving to ``index`` when given)."""
        mappy = _load_mappy()
        if index is not None and index.exists():
            source, save = index, None
        elif genome is not None:
            source, save = genome, index
        else:
            raise MadcError(f"Genome index not found: {index} (pass --genome to build it)")
        if not source.exists():
            raise MadcError(f"Genome file not found: {source}")
        # The file may be an index whichever option named it; an index is loaded, never rebuilt.
        parts = _index_parts(source)
        if parts is not None:
            save = None
        if save is not None:
            save.parent.mkdir(parents=True, exist_ok=True)
        if parts is not None and len(parts) > 1:
            raise MadcError(
                f"{source} is a minimap2 index in {len(parts)} parts ({sum(map(len, parts)):,} contigs), and mappy "
                f"searches only the first ({len(parts[0]):,} contigs), so the rest of the genome would be missed. "
                "Rebuild it as one part: give --genome and a new --genome-index path, or run minimap2 -x sr "
                "with -I above the genome size"
            )
        verb = "Loading index" if parts is not None else "Indexing genome"
        print(f"{verb} {source}" + (f" (saving {save})" if save else ""), file=sys.stderr)
        aligner = mappy.Aligner(
            str(source), preset="sr", n_threads=threads, fn_idx_out=str(save) if save else None,
            scoring=list(SR_SCORING), **SR_MAP_OPTIONS,
        )
        if not aligner:
            raise MadcError(f"Could not load or build a minimap2 index from {source}")
        if (aligner.k, aligner.w) != (21, 11):
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
        return self._aligner.seq(contig, start, end) or ""

    def window(self, contig: str, start: int, end: int) -> "WindowIndex":
        return WindowIndex(self._mappy, contig, start, self.fetch(contig, start, end))


class WindowIndex:
    """A sensitive index (short seeds) of one genome window, so a tag too
    diverged to seed genome-wide (k=21) still aligns at a site it belongs to.
    Each window is its own index, so windows never compete for hits; hits
    come back in genome coordinates."""

    def __init__(self, mappy_module, contig: str, start: int, sequence: str) -> None:
        self._contig, self._start = contig, start
        self._aligner = (
            mappy_module.Aligner(seq=sequence, preset="sr", k=SITE_INDEX_K, w=SITE_INDEX_W, scoring=list(SR_SCORING))
            if sequence
            else None
        )

    def hits(self, sequence: str, buffer=None) -> List[GenomeHit]:
        if self._aligner is None:
            return []
        return [_genome_hit(hit, self._contig, self._start) for hit in self._aligner.map(sequence, buf=buffer)]


def map_sequences(
    aligner, sequences: Sequence[str], jobs: int = 1, keep: Optional[int] = GENOME_SITES_KEPT
) -> Dict[str, GenomePlacement]:
    """Every sequence's distinct genome sites (``summarize_hits``).  Worker threads share
    the index; mappy releases the GIL while it maps, and each batch has its own buffer."""

    def place(batch: Sequence[str]) -> List[GenomePlacement]:
        buffer = aligner.new_buffer()
        return [summarize_hits(aligner.hits(sequence, buffer), keep) for sequence in batch]

    batches = [sequences[start : start + MAPPING_BATCH] for start in range(0, len(sequences), MAPPING_BATCH)]
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        placements = [placement for found in pool.map(place, batches) for placement in found]
    return dict(zip(sequences, placements))


def clone_position(clone_id: str, contigs: Set[str]) -> Tuple[Optional[Tuple[str, int]], Optional[str]]:
    """The CHROM_POS a CloneID names (1-based; zero padding allowed) and None, or None and why it names none."""
    match = CLONE_POSITION_RE.fullmatch(clone_id)
    if not match:
        return None, "no_position"
    if match.group(1) not in contigs:
        return None, "contig_absent"
    return (match.group(1), int(match.group(2))), None


def _site_covers(site: GenomeHit, contig: str, position: int) -> bool:
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

    With ``use_clone_position`` the CloneID position is taken to be on this
    genome.  It picks one of several equal-best sites, and a window around it
    stands in for a Ref tag that does not align.  When no equal-best site
    covers it, ``align_ref_tag(contig, start, end)`` aligns the Ref tag within
    a window around the position, and an alignment covering the position
    becomes the marker's site even if the tag aligns better elsewhere: the
    assembly may carry the Alt allele at the target and a Ref-like copy
    elsewhere, or the genome-wide search may have missed the site.  With no
    such alignment the marker is unplaced.  Either way the disagreement is
    reported.  Without ``use_clone_position`` (a genome other than the one the
    CloneIDs refer to) the position is never read: an unaligned Ref tag leaves
    the marker unplaced, and equal-best sites stay multi-copy.
    """
    position, no_position_reason = clone_position(clone_id, contigs) if use_clone_position else (None, "ignored")
    equal_best = ref_placement.best if ref_placement is not None else ()
    if position is None:
        if not equal_best:
            return MarkerPlacement((), "unplaced", 0, no_position_reason)
        source = "ref_tag" if len(equal_best) == 1 else "ref_tag_multicopy"
        return MarkerPlacement(equal_best, source, len(equal_best), no_position_reason)
    contig, base = position
    # The genome interval a tag of this length lies in if it contains the CloneID position.
    start, end = max(0, base - 1 - tag_length), base + tag_length
    if not equal_best:
        return MarkerPlacement((GenomeHit(contig, start, end, 0, 0),), "clone_id", 0, "fallback")
    at_position = [site for site in equal_best if _site_covers(site, contig, base)]
    if len(at_position) == 1:
        source = "ref_tag" if len(equal_best) == 1 else "ref_tag+clone_id"
        return MarkerPlacement(tuple(at_position), source, len(equal_best), "agree")
    if at_position:
        return MarkerPlacement(equal_best, "ref_tag_multicopy", len(equal_best), "agree")
    realigned = align_ref_tag(contig, max(0, start - SITE_PAD), end + SITE_PAD) if align_ref_tag else ()
    covering = [hit for hit in realigned if _site_covers(hit, contig, base)]
    if covering:
        return MarkerPlacement((max(covering, key=lambda hit: hit.score),), "clone_id_site", len(equal_best), "disagree")
    return MarkerPlacement((), "unplaced", len(equal_best), "disagree")


def score_at_sites(
    aligner, markers: Mapping[str, MarkerPlacement], haplotypes: Iterable[Tuple[str, str]], jobs: int = 1
) -> Dict[Tuple[str, str], int]:
    """Each (CloneID, sequence) pair's best score within its own marker's sites (each padded by SITE_PAD).

    Every placed marker's sites are indexed on their own with short seeds,
    so a real allele carrying several SNPs, which the genome-wide search
    misses or aligns only partly at home, is still scored there, whatever
    other markers the panel has.
    """
    by_marker: Dict[str, Set[str]] = defaultdict(set)
    for clone, sequence in haplotypes:
        if markers[clone].sites:
            by_marker[clone].add(sequence)

    def score(clone: str) -> Dict[Tuple[str, str], int]:
        windows = [
            aligner.window(site.contig, max(0, site.start - SITE_PAD), site.end + SITE_PAD) for site in markers[clone].sites
        ]
        buffer = aligner.new_buffer()
        found = {}
        for sequence in by_marker[clone]:
            scores = [hit.score for window in windows for hit in window.hits(sequence, buffer)]
            if scores:
                found[(clone, sequence)] = max(scores)
        return found

    results: Dict[Tuple[str, str], int] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
        for found in pool.map(score, list(by_marker)):
            results.update(found)
    return results


def genome_label(placement: GenomePlacement, marker: MarkerPlacement, site_score: Optional[int] = None) -> GenomeLabel:
    """Label a haplotype by its best score at its marker's own site against its best score anywhere else.

    ``on_target``: the marker's site scores strictly best.  ``paralog``:
    somewhere else scores strictly better (or the haplotype aligns only
    elsewhere).  ``ambiguous``: a tie between the marker's site and another
    site, or a marker whose own Ref tag is multi-copy and unresolved.
    ``unmapped``: no alignment.  ``unplaced``: the marker itself could not be
    located, so the haplotype cannot be judged.  ``site_score`` is the
    sensitive realignment at the marker's site (``score_at_sites``).  The
    target score is the higher of it and the genome-wide hit at the site, so
    the label does not depend on which search found the site.
    """
    at_marker = [hit for hit in placement.sites if any(_sites_overlap(hit, site) for site in marker.sites)]
    elsewhere = [hit for hit in placement.sites if hit not in at_marker]
    off_score = elsewhere[0].score if elsewhere else None
    off_sites = sum(hit.score == off_score for hit in elsewhere)
    off_hit = elsewhere[0] if elsewhere else None
    if not marker.sites:
        return GenomeLabel("unplaced", None, ".", off_score, off_sites, off_hit)
    if site_score is not None and (not at_marker or site_score > at_marker[0].score):
        target_score, target_from = site_score, "site"
    elif at_marker:
        target_score, target_from = at_marker[0].score, "genome"
    else:
        target_score, target_from = None, "."
    if target_score is None and off_score is None:
        label = "unmapped"
    elif target_score is None or (off_score is not None and off_score > target_score):
        label = "paralog"
    elif off_score == target_score or len(marker.sites) > 1:
        label = "ambiguous"
    else:
        label = "on_target"
    return GenomeLabel(label, target_score, target_from, off_score, off_sites, off_hit)


def _madc_stream(handle, path: Path):
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
    for name in (GENOME_LABEL_HEADER, GENOME_LABEL_COLUMN):
        if name in stripped:
            return stripped.index(name)
    raise MadcError(
        f"{madc_path} has no {GENOME_LABEL_COLUMN} column to carry the flag; use drop mode, "
        "or read the labels sidecar"
    )


def check_madc_header(madc_path: Path, flagging: bool) -> None:
    """Stop with an error when the MADC has no header row, no identity columns, or (to flag rows) no column for the flag."""
    with open_text(madc_path, "rt") as handle:
        _, header, _ = _madc_stream(handle, madc_path)
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
        _, header, reader = _madc_stream(handle, madc_path)
        for _, allele_id, clone, sequence in _haplotype_rows(reader, header, madc_path):
            ungapped = sequence.replace("-", "")
            haplotypes[(clone, sequence)] = ungapped
            tag_lengths.setdefault(clone, len(sequence))
            if clone not in ref_tags and allele_type(allele_id) == "Ref":
                ref_tags[clone] = ungapped
    return header, haplotypes, ref_tags, tag_lengths


def _place_markers(
    aligner, ref_tags: Mapping[str, str], tag_lengths: Mapping[str, int], use_clone_position: bool, jobs: int
) -> Dict[str, MarkerPlacement]:
    """Every marker's placement (``place_marker``), from the genome sites of its Ref tag."""
    # Ref tags keep every equal-best site, so a capped list cannot hide the one the CloneID names.
    ref_placements = map_sequences(aligner, sorted(set(ref_tags.values())), jobs, keep=None)
    markers: Dict[str, MarkerPlacement] = {}
    for clone, tag_length in tag_lengths.items():
        ref_tag = ref_tags.get(clone)   # None for a marker without a Ref row, which place_marker never aligns

        def align_ref_tag(contig: str, start: int, end: int, ref_tag=ref_tag) -> List[GenomeHit]:
            return aligner.window(contig, start, end).hits(ref_tag)

        markers[clone] = place_marker(
            clone, ref_placements.get(ref_tag), aligner.contigs, tag_length, use_clone_position, align_ref_tag
        )
    return markers


def _format_site(site: GenomeHit) -> str:
    strand = {1: "+", -1: "-"}.get(site.strand, ".")
    return f"{site.contig}:{site.start + 1}-{site.end}({strand})"


def _dot_if_none(value) -> object:
    return "." if value is None else value


def _write_markers_sidecar(markers_path: Path, markers: Mapping[str, MarkerPlacement]) -> None:
    markers_path.parent.mkdir(parents=True, exist_ok=True)
    with markers_path.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(GENOME_MARKER_SIDECAR_COLUMNS)
        for clone, marker in markers.items():
            sites = ";".join(_format_site(site) for site in marker.sites) or "."
            writer.writerow([clone, marker.source, marker.ref_tag_sites, marker.clone_check, sites])


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

    ``actions`` gives what to do with ``paralog`` and ``ambiguous`` rows
    (``drop``, ``flag`` or ``keep``); rows with any other label are kept.  A
    canonical Ref or Alt row that would be dropped is kept and reported as
    ``protected``.  With ``label_index`` every row written carries its label
    in that column.
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
        preamble, header, reader = _madc_stream(handle, madc_path)
        writer = csv.writer(out_handle, lineterminator="\n")
        label_writer = csv.writer(label_handle, lineterminator="\n")
        label_writer.writerow(GENOME_LABEL_SIDECAR_COLUMNS)
        if provenance and preamble and preamble[0] and preamble[0][0].startswith("* MADCATS"):
            preamble[0] = [f"{preamble[0][0]} | {provenance}"] + preamble[0][1:]
        writer.writerows(preamble)
        if label_index is not None:
            header[label_index] = GENOME_LABEL_HEADER
        writer.writerow(header)
        sample_start = madc_sample_start(header)
        for row, allele_id, clone, sequence in _haplotype_rows(reader, header, madc_path):
            result = labels[(clone, sequence)]
            category = allele_type(allele_id)
            action = actions.get(result.label, "keep")
            if action == "drop":
                outcome = "protected" if category in ("Ref", "Alt") else "dropped"
            else:
                outcome = "flagged" if action == "flag" else "kept"
            counts = madc_read_counts(row[sample_start:], madc_path)
            overflow.add(clone, counts)
            reads = sum(counts)
            for bucket in (totals["input"], by_label[result.label], by_outcome[outcome], by_category[category][result.label]):
                bucket["rows"] += 1
                bucket["reads"] += reads
            if outcome != "dropped":
                if label_index is not None:
                    row = row + [""] * max(0, label_index + 1 - len(row))
                    row[label_index] = result.label
                writer.writerow(row)
                totals["output"]["rows"] += 1
                totals["output"]["reads"] += reads
            off_hit = _format_site(result.off_target_hit) if result.off_target_hit is not None else "."
            label_writer.writerow([
                clone, allele_id, sequence, result.label, outcome, _dot_if_none(result.target_score), result.target_from,
                _dot_if_none(result.off_target_score), result.off_target_sites, off_hit,
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
    site_scores = score_at_sites(aligner, markers, {(clone, ungapped) for (clone, _), ungapped in haplotypes.items()}, jobs)
    labels = {
        (clone, sequence): genome_label(placements[ungapped], markers[clone], site_scores.get((clone, ungapped)))
        for (clone, sequence), ungapped in haplotypes.items()
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
    sources = Counter(marker.source for marker in markers.values())
    checks = Counter(marker.clone_check for marker in markers.values())
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
            "multicopy_ref_tags": [clone for clone, marker in markers.items() if marker.ref_tag_sites > 1],
            "clone_id_disagreements": [clone for clone, marker in markers.items() if marker.clone_check == "disagree"],
            "unplaced": [clone for clone, marker in markers.items() if marker.source == "unplaced"],
        },
        "haplotypes": {
            "distinct": len(haplotypes),
            "distinct_sequences": len(sequences),
            "labels": dict(Counter(result.label for result in labels.values())),
            "realigned_at_site": {
                "haplotypes": sum(bool(markers[clone].sites) for clone, _ in haplotypes),
                "found_at_site": sum((clone, ungapped) in site_scores for (clone, _), ungapped in haplotypes.items()),
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
        "genome_label_column": GENOME_LABEL_HEADER if label_index is not None else None,
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
    dropped = rows["input"]["reads"] - rows["output"]["reads"]
    share = dropped / rows["input"]["reads"] if rows["input"]["reads"] else 0.0
    print(
        f"Kept {rows['output']['rows']:,} of {rows['input']['rows']:,} rows; dropped {dropped:,} reads ({share:.1%})",
        file=sys.stderr,
    )
    print(f"Wrote {output}, its labels {labels_path}, and QC to {qc_path}", file=sys.stderr)
    return 0


def add_filter_paralogs_command(subparsers) -> None:
    filtering = subparsers.add_parser(
        "filter-paralogs",
        help="drop or flag MADC rows that a reference genome places at another locus (needs mappy)",
        description=(
            "Clean a MADC with a reference genome, for panels where paralogous copies of a target co-amplify "
            "with it. Each marker is placed by aligning its Ref tag, with the CloneID position as a "
            "cross-check. Every haplotype is then aligned and labelled on_target, paralog (it aligns better "
            "elsewhere) or ambiguous (it aligns equally well at the marker and elsewhere). Paralog and "
            "ambiguous rows are dropped or flagged; canonical Ref and Alt rows are always kept. Reads MADCs "
            "from generate, vendor DArT MADCs and HapApp fixed-allele-ID MADCs. Needs the mappy aligner: "
            "pip install 'madcats[genome]'."
        ),
    )
    filtering.add_argument(
        "--madc",
        required=True,
        help=(
            "MADC to clean: from generate, a vendor DArT MADC (with or without its starred preamble), or a "
            "HapApp fixed-allele-ID MADC"
        ),
    )
    filtering.add_argument("--genome", help="reference genome FASTA (optionally gzipped)")
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
            "CHROM_POS CloneIDs: use them as positions on the genome given, to choose between equal-best Ref-tag "
            "sites, stand in for Ref tags that do not align, and place a marker whose Ref tag aligns better "
            "elsewhere at its CloneID position (or leave it unplaced if the tag does not align there); or ignore "
            "them when they refer to another assembly (default: use)"
        ),
    )
    filtering.add_argument("--output", required=True, help="cleaned MADC to write")
    filtering.add_argument(
        "--paralogs",
        choices=("drop", "flag"),
        default="drop",
        help="paralog rows: drop them, or keep them labelled in the GenomeLabel column (default: drop)",
    )
    filtering.add_argument(
        "--ambiguous",
        choices=("drop", "flag", "keep"),
        default="drop",
        help="rows with equal-best hits at the marker and elsewhere: drop, flag, or keep unlabelled (default: drop)",
    )
    filtering.add_argument("--qc-output", help="QC JSON (default: OUTPUT.qc.json)")
    filtering.add_argument(
        "--jobs", type=positive_int, default=_default_jobs(), help="alignment threads (default: up to 4)"
    )
    filtering.set_defaults(func=run_filter_paralogs)
