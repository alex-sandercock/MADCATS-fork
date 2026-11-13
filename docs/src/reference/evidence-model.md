# Evidence Model

This page explains how `madcats generate` decides which reads count, which sequences become rows,
and how sequencing errors are removed. The design rationale and the measurements behind each default
are in the [development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md) on GitHub.

## Read Assignment

For each marker, the panel holds a fixed-length Ref tag and Alt tag that both start with an
assignment anchor (the left primer, or a pseudo-primer prefix).

1. A read is assigned to a marker by its prefix (the anchor). One mismatch is allowed by default
   (`--primer-mismatches 1`). Anchors that are prefixes of longer anchors are supported, and every
   anchor within the mismatch limit is considered. A read that matches several markers goes to the
   one whose full tag it matches best; if that is a tie, the read is dropped as ambiguous. The one
   exception is a tie within a family of [tiled markers](#tiled-markers).
2. The read's first tag-length bases are kept if they have at least `--min-identity` (80%) identity
   to the marker's Ref or Alt tag, counting substitutions only.
3. If the read fails for every candidate marker, it gets one gapped alignment under the same 20%
   edit limit, with gaps of up to `--max-indel-length` (10) bases. The alignment is kept only if it
   looks like an indel allele: at most 6 substitutions and 6 gap events. Reads from co-amplified
   paralogs need many substitutions to fit, so they are rejected. A fast bit-vector edit-distance
   bound skips tags that cannot match before any alignment runs.
4. A read shorter than the tag is added to the one full-length allele whose start it matches
   exactly, but only if that allele has already passed the evidence rules and the read is longer
   than every prefix shared between two alleles of its locus. A partial read can therefore add to an
   allele's count, but it cannot create an allele or change an allele ratio. Reads shorter than
   `--min-partial-length` (tag length − 9) are never counted.

### Tiled markers

Tiles are markers that read the same tag window of one amplicon: DArTag panels place several
target SNPs in one amplicon, under one shared (or nested) primer, and write each tile's tag with
its neighbours' SNPs as IUPAC codes. A read from the amplicon fits every tile equally well, so the
tie rule above would drop all of them and leave the tiles with zero reads. `--tile-groups` (default
`auto`) links tiles into families, either declared with `LocusGroup` or detected; see
[Tiled Markers](../getting-started/input-formats.md#tiled-markers) for the same-window test.

When every marker tied for the best fit belongs to one family, the read came from that family's
amplicon. It is credited to every member it fits within the identity limit, by substitution-only
identity or through one gapped alignment when gapped rescue is on. Each member
records the read's exact tag as its own observation. This holds whether or not the read tied; a
family member it fits less well still gets it. So each tile gets the reads it would get in a panel
of its own, as in DArT's pipeline. Partial reads follow the same rule. A tie between markers of different families, or between a family member and an unlinked
marker, is still dropped as ambiguous.

Tiles share their reads, so their counts are not independent: the same molecules support every
member's row. Use the tiles' calls as one locus, or keep one tile per family, when an analysis
assumes independent markers.

## Which Sequences Become Rows

Canonical Ref and Alt rows are always written, even with zero counts.

Any other sequence becomes a row when it reaches `--min-global-count` reads project-wide **and**
passes at least one enabled evidence route:

| Route | Default | Rule |
|-------|---------|------|
| **Reproducibility** (primary) | on | At least `--reproducible-min-reads` (10) reads in each of at least `--reproducible-min-samples` (2) samples |
| **Depth** | on | At least `--depth-min-reads` (30) reads within a single sample. Keeps private alleles discoverable in unreplicated designs |
| **Quality** | off | At least `--quality-min-clean-reads` reads whose every tag base reaches `--discovery-min-quality` (Q20) |

Routes combine by OR. Identical tag sequences that recur independently across samples are strong
evidence, even when one base position is systematically low quality, which is why reproducibility
is primary and the whole-tag quality route is off by default.

### Supporting-frequency floor

`--min-supporting-freq` (default 0, off) additionally requires a sample's reads for the sequence to
reach that share of the sample's locus total before the sample can support the reproducibility or
depth route.

- At very high depth (above roughly 500x per sample per marker, e.g. AgriSeq), use `0.10`. At these
  depths sequencing errors pass absolute thresholds, and per-marker depth varies tenfold within a
  run, so a share-based floor works better than raising the absolute thresholds.
- For unreplicated designs that rely on the depth route, around `0.01` is a reasonable safeguard.
- Otherwise, leave it off. Co-amplified paralogs add reads to a locus, so real minor haplotypes can
  fall below 1% of its reads.

## Denoising

Canonical sequences and sequences that pass the evidence rules are protected: denoising never merges
them into anything else.

1. Per-sample denoising (UNOISE-style abundance skew). Within each sample, a rare sequence that differs at
   no more than `--denoise-distance` (1) bases from a more abundant parent is merged into the parent when the
   parent is at least `2^(alpha·d + 1)` times more abundant (`--denoise-alpha` 2.0). The closest
   parent is chosen first, then the most abundant.
2. Pooled same-allele merge (`--global-denoise-distance`, 1). A sequence that failed the evidence
   rules is merged into the protected parent chosen from project-wide abundances, but only when both
   carry the same target allele, so no merge can change a sample's Ref/Alt dosage.

Reads in sequences that pass no route and have no parent are dropped, and the QC read accounting
reports them. Setting every route threshold to 0 writes every sequence with at least
`--min-global-count` reads (exhaustive discovery).

The same evidence rule decides discovery, denoising protection and partial-read crediting, so the
three stay consistent.

## Classification

Each row is classified against its marker's canonical tags:

| Category | Meaning |
|----------|---------|
| `Ref` | The canonical Ref tag |
| `Alt` | The canonical Alt tag |
| `RefMatch` | Carries the Ref target allele plus other differences |
| `AltMatch` | Carries the Alt target allele plus other differences |
| `Other` | Neither target allele: the target base is a third base, missing, or deleted |

The target base is normally read at its fixed offset in the tag. When a gapped alignment explains
the row with fewer edits, stays within the substitution and gap limits, and is confirmed by at least
12 bases after its last gap, the target base is read from that alignment instead. A target inside a
deletion, or pushed past the end of the read, is `Other`.

## IUPAC Codes and N Bases

Classification respects ambiguity codes in the design: if a design allele carries an ambiguity
code, an observed base compatible with it classifies as RefMatch/AltMatch.

Row identity and counts use exact sequence equality, so an ambiguity-coded canonical tag reports
zero exact counts. MADCATS does not guess which concrete variants a curated IUPAC consensus was
meant to cover, because such codes can stand for real secondary SNPs, and collapsing them would
merge distinct microhaplotypes.

A read containing N never becomes a row of its own. It is added to the one N-free sequence at its
locus that it is compatible with, or dropped (`reads_n_mapped` / `reads_n_dropped` in the QC).

## Indels

Designed indel markers are fully supported. Vendors report indel alleles as fixed-length windows
whose suffix shifts at the indel, so both tags are ordinary same-length sequences.

Undesigned indels are recovered by the gapped alignment step above and written as fixed-length
rows in the same shifted-suffix convention. Indels longer than `--max-indel-length` are not
recovered, and neither are alignments that need many substitutions, which usually come from
paralogs.

## Read Accounting

The QC JSON accounts for every read: FASTQ records, reads assigned at full length, reads assigned
through a gapped alignment (`reads_assigned_gapped`), partial reads credited, ambiguous or unmapped,
reads moved by each denoising pass, reads dropped with low-evidence sequences, and per-sample
`reads_in_emitted_rows`.

A read credited to several tiles counts once in `reads_assigned_full_length` (per sample:
`reads_assigned`). The further credits are counted separately, as `tile_extra_credits_full_length`
and `tile_extra_credits_partial` (per sample: `tile_extra_credits`, `partial_tile_extra_credits`).
From assignment on, the ledger counts (marker, read) credits: `reads_n_mapped`, `reads_n_dropped`,
`partial_reads_added_to_counts`, `reads_in_final_count_table`, the denoising moves and
`reads_in_emitted_rows`. For a panel without tiles, credits and reads are the same. The ledger
balances as

```
reads_in_final_count_table = reads_assigned_full_length + tile_extra_credits_full_length
                             − reads_n_dropped + partial_reads_added_to_counts
```
