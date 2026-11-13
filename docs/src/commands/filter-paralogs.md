# madcats filter-paralogs

Drop or flag MADC haplotypes that a reference genome places at another locus. Optional; needs the
`mappy` aligner.

## Synopsis

```bash
madcats filter-paralogs --madc <FILE> --genome <FASTA> [--genome-index <MMI>] --output <FILE> [OPTIONS]
madcats filter-paralogs --madc <FILE> --genome-index <MMI> --output <FILE> [OPTIONS]
```

## Description

In some panels, paralogous copies of a target co-amplify with it, and their reads land in the
marker's rows. A paralog that differs from the Ref tag only at non-target bases becomes a RefMatch
or AltMatch row. When counts are summed per target SNP (as `BIGr::madc2vcf_all` does), those rows
push homozygous samples towards false heterozygous calls. `filter-paralogs` uses a reference genome
to decide where each haplotype comes from:

1. **Place each marker.** Its Ref tag is aligned to the genome. A CloneID of the form `CHROM_POS`
   (zero padding allowed, e.g. `Chr1_000333995`) is read as the marker's position on that genome.
   When the Ref tag has several equal-best sites, the CloneID position picks one; if it cannot, the
   marker is reported as multi-copy. When the Ref tag does not align, a window around the CloneID
   position is used. When the Ref tag's best sites are all away from the CloneID position, the tag
   is aligned again within that window (see [CloneID Positions](#cloneid-positions)). All of these
   uses assume the CloneID refers to the genome given; see
   [Using Another Assembly](#using-another-assembly).
2. **Align every distinct haplotype.** Each haplotype gets two scores: its best alignment at its
   marker's site, and its best alignment anywhere else. Every haplotype is also realigned to its own
   marker's site alone with short seeds, and the better of that score and the genome-wide one at
   the site counts. A real allele carrying several SNPs, which the genome-wide search misses or
   aligns only partly, is still scored at home.
3. **Label it**:

   | Label | Meaning | Default action |
   |-------|---------|----------------|
   | `on_target` | The marker's site scores strictly best | kept |
   | `paralog` | Another site scores strictly better, or the haplotype aligns only elsewhere | dropped (`--paralogs`) |
   | `ambiguous` | The marker's site ties with another site, or the marker's own Ref tag is multi-copy and unresolved | dropped (`--ambiguous`) |
   | `unmapped` | No alignment anywhere | kept |
   | `unplaced` | The marker itself could not be placed | kept |

4. **Write the cleaned MADC.** Canonical Ref and Alt rows are never dropped, because downstream
   dosage tools need both designed alleles. A canonical row whose label would have dropped it is
   counted as `protected`. Kept rows are copied unchanged, and the output keeps the input's layout:
   the vendor preamble, the statistics columns, HapApp's `Code_version` column and the column order.

The filter is opt-in, meant for panels known to co-amplify paralogs. Check the share of reads in
`paralog` rows in the QC JSON: in the lettuce 3K panel it was 63% of all reads in the vendor MADC.
On panels without the problem the filter removes little and is not needed.

It reads three layouts: a MADC from `generate`, a vendor DArT MADC, and the fixed-allele-ID MADC
that HapApp writes (`Code_version`, `AlleleID`, `CloneID`, `AlleleSequence`, then the samples, with
IDs such as `|Ref_0001`). It can therefore run on the raw MADC before HapApp, or on HapApp's output.

Scores use minimap2's short-read scoring (match +2, mismatch −8, gap open −12 and extension −2,
ambiguous base −1), recomputed from each alignment, so ties are exact.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--madc <FILE>` | MADC to clean: from `generate`, a vendor DArT MADC (with or without its starred preamble), or a HapApp fixed-allele-ID MADC. Gzipped input is read directly |
| `--output <FILE>` | Cleaned MADC to write (must differ from `--madc`) |

Give `--genome`, `--genome-index`, or both.

## Genome

| Option | Default | Description |
|--------|---------|-------------|
| `--genome <FASTA>` | none | Reference genome FASTA, optionally gzipped |
| `--genome-index <MMI>` | none | minimap2 index. Loaded when it exists; otherwise built from `--genome` and saved here for the next run |
| `--clone-positions {use,ignore}` | `use` | `CHROM_POS` CloneIDs: `use` them as positions on the genome given, to choose between equal-best Ref-tag sites, stand in for Ref tags that do not align and place markers whose Ref tag aligns better elsewhere (see [CloneID Positions](#cloneid-positions)), or `ignore` them when they refer to another assembly |

Indexing a large genome is slow and memory-heavy, so build the index once and reuse it. The first
run with both options writes it:

```bash
madcats filter-paralogs --madc MADC.csv --genome genome.fa --genome-index genome.sr.mmi --output MADC_filtered.csv
```

Later runs need only `--genome-index`. An index from the minimap2 command line works too, built with
the short-read preset in one part, with `-I` above the genome size:
`minimap2 -x sr -I 20G -d genome.sr.mmi genome.fa`. minimap2 splits the index of a genome larger
than `-I` (8 Gb by default) into parts, and mappy searches only the first part, so the command
refuses a multi-part index rather than miss the rest of the genome. An index the command builds
itself is always one part. On the 2.4 Gb lettuce genome indexing took 50 s and 10.6 GB of memory,
and the index file is 5.7 GB.

## Options

| Option | Default | Description |
|--------|---------|-------------|
| `--paralogs {drop,flag}` | `drop` | `paralog` rows: drop them, or keep them labelled in the `GenomeLabel` column |
| `--ambiguous {drop,flag,keep}` | `drop` | `ambiguous` rows: drop, flag, or keep. A kept row still carries its label when `--paralogs flag` adds the `GenomeLabel` column |
| `--qc-output <FILE>` | `OUTPUT.qc.json` | QC JSON |
| `--jobs <INT>` | up to 4 | Alignment threads |

Dropping ambiguous rows is the default because, on lettuce, keeping them left false heterozygous
calls at 3.4% of inbred genotypes, against 0.87% when they were dropped.

## Flag Mode

With either option set to `flag`, the `PICRef` column is renamed `GenomeLabel` and every row
carries its label there: `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced`. Flagged
rows are kept. Column positions do not change, so positional consumers still see the 16-column
layout. In a vendor MADC this overwrites the vendor's `PICRef` statistic. Flag mode needs the
`PICRef` (or an earlier `GenomeLabel`) column. HapApp's fixed-allele-ID MADC has none, so on that
layout use drop mode and read the labels sidecar.

A flagged file is for review. Tools that check the vendor's column names, such as HapApp-online's
upload check, reject `GenomeLabel` in place of `PICRef`; give them the drop-mode output.

## Outputs

| File | Contents |
|------|----------|
| `OUTPUT` | The cleaned MADC |
| `OUTPUT.genome_labels.csv` | One line per haplotype row, dropped rows included. A row with an empty `AlleleID` is not a haplotype row; it is passed over and not written to `OUTPUT` either |
| `OUTPUT.genome_markers.csv` | One line per marker: where it was placed and how |
| `OUTPUT.qc.json` | Parameters, aligner settings, marker placement summary, the number of distinct haplotypes and sequences, haplotypes per label and how many were found again in the realignment at their marker's site, alignment time, rows and reads per label, per outcome and per allele category, and `count_overflow` (see below) |

Read totals are summed as stored. If the input has negative counts, the command warns and
`count_overflow` records how many cells and markers carry them. A negative count with a maximum of
exactly 127 points to signed 8-bit overflow in a vendor file; see the
[FAQ](../reference/faq.md#my-vendor-madc-has-negative-read-counts).

When the input carries a MADCATS provenance cell (`generate --dart-preamble`), the filter's
settings are appended to it.

### Labels sidecar

| Column | Contents |
|--------|----------|
| `CloneID`, `AlleleID`, `AlleleSequence` | Row identity |
| `Label` | `on_target`, `paralog`, `ambiguous`, `unmapped` or `unplaced` |
| `Action` | `kept`, `flagged`, `dropped`, or `protected` (a canonical row kept despite its label) |
| `TargetScore` | Best alignment score at the marker's site (`.` if none) |
| `TargetFrom` | `genome` (the genome-wide search's hit at the site), `site` (the realignment to the marker's site, when it scored higher or the genome-wide search found nothing there) or `.` |
| `OffTargetScore` | Best alignment score anywhere else (`.` if none) |
| `OffTargetSites` | Distinct sites elsewhere at that score, among the 8 sites kept for a haplotype |
| `OffTargetHit` | The best site elsewhere, 1-based: `Chr4:61708465-61708545(+)` |

A perfect 81-base match scores 162 and each mismatch costs 10, so `TargetScore` and
`OffTargetScore` show how many mismatches separate a haplotype from each copy.

### Markers sidecar

| Column | Contents |
|--------|----------|
| `CloneID` | Marker |
| `Placement` | `ref_tag`, `ref_tag+clone_id` (multi-copy Ref tag resolved by the CloneID position), `ref_tag_multicopy` (unresolved), `clone_id_site` (the Ref tag's alignment at the CloneID position, though it aligns better elsewhere), `clone_id` (Ref tag did not align) or `unplaced` |
| `RefTagSites` | Equal-best genome sites of the Ref tag (0 when it did not align). The genome-wide search returns at most 51 hits, so this counts the sites found, not necessarily every copy |
| `CloneIDCheck` | `agree`, `disagree`, `no_position` (CloneID is not `CHROM_POS`), `contig_absent`, `fallback` or `ignored` (`--clone-positions ignore`) |
| `Sites` | The marker's site(s), 1-based, `;`-separated |

## CloneID Positions

With `--clone-positions use` (the default), a `CHROM_POS` CloneID is read as the marker's position
on the genome given, where the panel was designed. When none of the Ref tag's equal-best sites is
at that position (`CloneIDCheck` `disagree`), the Ref tag is aligned again, alone, within a window
around the position. If it aligns there, that alignment is the marker's site (`clone_id_site`),
even when the tag scores better elsewhere: the genome-wide search can miss the design site, and
the assembly may carry the Alt allele at the target while a paralog carries Ref. If the tag does
not align at the position, the marker is `unplaced` and keeps every row.

The command warns when more than 10% of markers disagree. Either the genome is not the one you
meant to use, or the CloneIDs refer to another assembly, in which case rerun with
`--clone-positions ignore`.

## Using Another Assembly

Any genome that contains the panel's target loci can be used, not only the one the panel was
designed on, because markers are placed by aligning their Ref tags. Labels against a second
assembly, such as another cultivar's or a newer version, show where the two genomes support a
haplotype differently. A haplotype that is `on_target` against one genome and `paralog` against
another aligns better to some other locus in the second genome. That locus may be a copy the first
genome lacks, but a copy present in both with a different sequence changes the label in the same
way, so confirm a copy difference in the assemblies themselves.

CloneID positions refer to the design assembly, so against another assembly use
`--clone-positions ignore`. Otherwise, where the chromosome names match, a CloneID position is read
in the new genome's coordinates. A marker whose Ref tag does not align is then judged at the wrong
place, and its non-canonical rows are dropped as paralogs. A marker whose Ref tag aligns away from
its CloneID position is left `unplaced`, unless the tag also aligns at that position. With
`ignore`, a marker whose Ref tag does not align is `unplaced` and keeps every row, and a Ref tag
with several equal-best sites stays multi-copy, so its rows are `ambiguous`. Where the chromosome
names differ, CloneID positions go unused anyway (`contig_absent`).

Keep the default on the design assembly. On the lettuce vendor MADC, CloneID positions decide the
site of 11 of 2,987 markers: 9 whose Ref tag does not align, 1 tie, and `Chr2_001081815`, whose Ref
tag carries two IUPAC ambiguity codes. The genome-wide search placed that tag on Chr4 (score 96),
but it aligns at its CloneID position on Chr2 with 156. With `ignore`, the 9 markers are
`unplaced`, and 1,032 rows that the design assembly labels `paralog` are kept, while
`Chr2_001081815` is judged on Chr4, which drops 11 of its rows that the default keeps.

Each genome needs its own index:

```bash
madcats filter-paralogs --madc MADC.csv --genome-index v8.sr.mmi --output MADC_v8.csv
madcats filter-paralogs --madc MADC.csv --genome other.fa --genome-index other.sr.mmi \
  --clone-positions ignore --output MADC_other.csv
```

To compare the two runs, join their labels sidecars on `CloneID` and `AlleleSequence`:

```python
import pandas as pd

key = ["CloneID", "AlleleSequence"]
design = pd.read_csv("MADC_v8.csv.genome_labels.csv")
other = pd.read_csv("MADC_other.csv.genome_labels.csv")
both = design.merge(other, on=key, suffixes=("_v8", "_other"))
print(pd.crosstab(both.Label_v8, both.Label_other))
```

Rows that are `on_target` against the design assembly and `paralog` against the other align
better to another locus in the other assembly. They point to loci worth checking for a copy that
only the other assembly has, not to such a copy by themselves.

## Example

```bash
madcats filter-paralogs \
  --madc DLet22-7611_MADC.csv \
  --genome-index ref/Lsativa_467_v8_Chr.sr.mmi \
  --output DLet22-7611_MADC_filtered.csv \
  --jobs 8
```

On the lettuce 3K DArT MADC (732,534 rows, 376 samples) this took 2 min 15 s and 9.4 GB of memory.
It kept 113,667 rows, with 35.7% of the reads. See the
[paralog filter tutorial](../tutorials/paralog-filter.md) for the whole workflow.

## Limitations

- The genome must contain the panel's target loci: the assembly the panel was designed on, or one
  close to it. Markers whose region is missing or collapsed in the assembly cannot be judged.
- A copy that exists in your germplasm but not in the reference assembly cannot be seen; its
  haplotypes look like diverged alleles of the target. A second assembly that has the copy can
  show it (see [Using Another Assembly](#using-another-assembly)).
- minimap2 finds fewer distant copies than bwa mem, so a haplotype far from every genome copy is
  judged on the hits it has. The rows kept or dropped agreed with a bwa-mem classification for
  98.2% of rows of a MADCATS MADC and 94.9% of a vendor MADC, whose rows are more diverged. Those
  figures count the canonical Ref and Alt rows, which are always kept; among the rows the filter
  decides, they are 97.3% and 94.8%.
- Only the tag is aligned (81 bases by default), not the full read.
- An allele of the marker that is as near a copy as it is to its own site is `ambiguous`. With
  copies planted one or two bases from the Ref tag, 8 of 349 rows of the markers themselves were
  lost this way; with copies that differ more, none of 515.
- In a polyploid whose reference carries every subgenome or haplotype, a tag has a copy on each,
  and those copies are treated as any other: reads of a diverged copy are labelled `paralog`, and
  an identical copy makes the marker's rows `ambiguous`. This has not been measured on a
  polyploid panel.

## How It Was Checked

- **Planted copies.** In the delta smelt genome (437 Mb) copies of 105 markers' tags were added
  as new contigs, with the rows a co-amplified copy would give added to a real MADC. Every row
  from a copy 1 to 8 bases, or a base and a 2-base deletion, from the Ref tag was labelled
  `paralog` (450 of 450 over three draws), and no row of the 263 untouched markers was mislabelled
  (3,348 rows). A copy left out of the genome was not seen.
- **Placement.** On the same data all 368 markers were placed where their two primers lie in the
  genome, with CloneID positions or without. Markers given a wrong position were left unplaced,
  not placed elsewhere.
- **Against bwa mem.** The lettuce figures above.
- **Alfalfa, against bwa mem and blastn.** Two tetraploid F1 families of the alfalfa 3K panel
  (7,517 rows of 54-base tags and 8,731 rows of 81-base tags) and the panel's haplotype database
  (23,469 haplotypes) were run on a monoploid genome of another cultivar than the one the panel
  was designed on. Every row was also aligned with bwa mem and with blastn and labelled from
  each aligner's own score. The filter kept or dropped the same rows as bwa mem for 99.4 to
  99.9% of rows and as blastn for 99.3 to 99.7%; the two aligners agree with each other on 99.4
  to 99.8%.
- **Alfalfa, by inheritance.** In those families a row can be mapped without aligning it, by the
  rows of other markers it is inherited with. Of the rows the filter kept, 97.6% map to their own
  marker in both families. Of the rows labelled `paralog` whose better site is far from the
  marker, 84 to 95% map away from it, most of them to the region the filter named. A row
  whose better site is a nearby copy on the same chromosome is inherited with the marker, and
  this test cannot judge it.
- **Alfalfa, placement.** With `--clone-positions ignore` on that genome, 97.8% of the markers
  with 81-base tags and 91.0% of those with 54-base tags were placed where their 361-base design
  sequence lies, 0.5 to 0.9% in another copy, and the rest left unplaced. With CloneIDs that are
  positions on the genome given, all but 7 of 2,734 were placed there.
- Output does not depend on `--jobs`.

## See Also

- [Paralog filter tutorial](../tutorials/paralog-filter.md)
- [madcats generate](generate.md)
- [Output Formats](../reference/output-formats.md)
