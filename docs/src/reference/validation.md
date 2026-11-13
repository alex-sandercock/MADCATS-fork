# Validation & Accuracy

MADCATS has been validated on five platforms, each against the best independent comparison available
for that dataset. This page summarizes the results. The full record, including every
version-to-version comparison, is in the [development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md) on GitHub.

## Performance at a Glance

"Usable windows" is the fraction of raw reads that end up counted in a panel window. Some of the
loss is off-target and primer-dimer reads that should be rejected, and other pipelines calculate
their percentages differently.

"Genotype agreement" means the MADCATS table was run through HapApp, BIGr and updog to call dosages,
and those calls were compared with the reference genotypes.

| Platform (dataset) | Reads → usable windows | Accuracy | Validated against |
|---|---|---|---|
| DArTag (sweetpotato, red clover, blackberry, alfalfa, pecan) | 56.5% of delivered records, vs the vendor's 50.9% on the same FASTQs | Sweetpotato count r = 0.9997 (1.20M cells); red clover dosage agreement 96.98%; blackberry 95.39%; alfalfa 98.40% / 97.40%; pecan 98.78% | Vendor MADC, vendor/updog dosage, families, skim-seq VCF |
| AgriSeq, Ion Torrent (2 honeybee chips; 372-sample oyster) | 88.4% | 94.9% vs Torrent VariantCaller (152,698 cells); honeybee replicate concordance 98.0% | Vendor caller + replicates |
| GT-seq (96-sample delta smelt) | 84.5% merged | 99.49% R1-only; 98.42% genotype agreement; `infer-design` 368/368 | Independent BWA/bcftools truth VCF |
| Flex-Seq capture (192-sample cranberry) | 81.9% | 95.13% genotype agreement (295,068 cells); trio incompatibility 3.6% vs 6.2% for the comparison VCF's source | Published population VCF + 20 pedigree trios |
| rhAmpSeq (608 runs, 8 grape families) | 86.6% | Family transmission 99.92–99.99%; 93.0% exact agreement with the platform's own pipeline on identical reads, with +19% / +79% more callable cells at 10x / 30x | Mendelian constraints + same-reads head-to-head |
| Synthetic truth suite | — | Dosage 100% on clean, linked, paralog and repeat classes; paired-end merge 100% with 98.6% molecular recovery | Exact simulated truth |

Accuracy depends more on the quality of the comparison data than on the platform. Where the
reference genotypes are reliable (GT-seq, DArTag), MADCATS agrees at 98–99% or better. Against
noisier references (Torrent VariantCaller, raw mpileup), agreement is 94–95%, and replicate and trio
checks show the MADCATS calls to be the more consistent of the two.

## DArTag

### Sweetpotato (hexaploid, 81-bp tags)

The development project: 80 samples, 3,120 markers, 154.6M reads. Against the vendor MADC (19,319
rows), every shared row had identical classification.

| Configuration | Rows | Precision | Recall |
|---|---|---|---|
| Default (reproducibility 10×2 OR depth 30×1) | 17,844 | 84.0% | 77.6% |
| `--reproducible-min-reads 20` | 16,517 | 89.8% | 76.8% |
| `--min-supporting-freq 0.01` | 17,319 | 86.2% | 77.2% |

Row recall against a vendor table has a ceiling: 27.3% of the vendor's positive-count rows have
fewer than 10 reads in fewer than 2 samples of the vendor's own counts, below any within-project
evidence threshold.

### Red clover (diploid, 109-bp tags, ~60x depth, no primer file)

- `infer-design` recovered 2,968 loci; 2,902 exactly match the vendor canonicals and 56 differ only
  by vendor IUPAC curation.
- Recall against vendor rows with enough read support to pass an evidence route: 94.1%.
- Updog dosage concordance with the vendor-derived calls: 99.73% at DP ≥ 10.
- In a 152-progeny mapping family, the Mendelian-violation rate was 1.30% under MADCATS against
  1.59% under the vendor MADC. The heterozygotes created by recovered indel alleles are inherited.

### Blackberry (mixed ploidy, 109-bp tags, indels, degenerate primers)

The first project validated through all three panel routes:
- Design files: 2,881 / 3,015 exact; all 134 others differ only by vendor IUPAC codes; all 10
  indel loci exact.
- Vendor-MADC bootstrap: 92.9% precision, and 94.7% recall against vendor rows with enough read
  support to pass an evidence route.
- Primer-free inference: 2,943 inferred; 18 wrong windows (0.6%), all at read-starved or
  short-amplicon loci.

### Alfalfa and pecan (raw deliveries)

- Alfalfa (two 2021 orders, 54-bp design): category agreement 99.60% / 99.83%, dosage agreement
  98.40% / 97.40%; F1 and BC1 progeny haplotypes within the parental sets at 98.6–98.7% (vendor
  98.9%).
- Pecan 3K (81-bp design): category agreement 98.89%, dosage agreement 98.78% on 1.02M cells; F1
  progeny within the parental set at 99.82% (vendor 99.86%); 93.3% agreement with skim-sequencing
  genotypes on 11 shared cultivars (vendor 93.3%).

## AgriSeq (Ion Torrent)

- Honeybee: two chips on an 81-marker panel. The curated frozen 41-bp panel scored 93.4% /
  93.9% vendor-caller concordance and 98.0% concordance across 30 technical replicate pairs.
  Freezing windows raised cross-chip haplotype sharing from 42% to 88%.
- Pacific oyster (PRJNA1118744; 372 samples, 591 amplicons): with the honeybee defaults
  unmodified, 88.4% of reads normalized and genotypes agreed with Torrent VariantCaller at 94.9%,
  on a different lab, species, chip generation and target length.

## GT-seq (delta smelt, Zenodo 10.5281/zenodo.19069550)

96 samples, 410 loci, and the panel's real primer file. Genotypes agreed with the archive's
independent BWA/bcftools VCF at 99.49%. Merged reads with frozen tags plus tiles kept that accuracy
on the original markers while adding 54% more truth-validated calls. `infer-design` with the primer
file withheld recovered 368/368 loci exactly.

## Flex-Seq Capture (cranberry)

192 samples with a K=4 tiled 110-bp panel and one cohort-level `generate`. Genotype agreement with
the published population VCF was 95.13%. Trio incompatibility was 3.6%, against
6.2% for the raw calls behind the comparison VCF.

## rhAmpSeq (grape, PRJNA281110)

Eight biparental families (763 runs) processed against one frozen 2,951-window panel. Across
~14,900 locus-family tests, 99.92–99.99% of haplotype calls in every family are alleles of that
locus's four inferred parental haplotypes. Against the platform's own `analyze_amplicon` pipeline
on identical reads, MADCATS was the more Mendelian-consistent of the pair (99.947% vs 99.850%).

## Adapters

Across 127,667 emitted rows in six MADCs spanning five platforms, zero rows contain TruSeq or Ion P1
adapter at a strict 13-mer probe. Trimming adapters before a run changed dosage concordance by 0.05
points on the capture benchmark while costing 2 points of extraction.

## Genome Paralog Filter (`filter-paralogs`)

- **Lettuce 3K DArTag** (the panel the filter was built on): the rows kept or dropped agreed with a
  bwa-mem classification for 98.2% of the rows of a MADCATS MADC and 94.9% of a vendor MADC; among
  the rows the filter decides, 97.3% and 94.8%.
- **Delta smelt GT-seq, placement:** all 368 markers were placed where their two primers lie in the
  437 Mb genome, and no row of this paralog-free panel was dropped.
- **Delta smelt GT-seq, planted copies:** with copies of 105 markers' tags added to the genome and
  the matching rows added to the MADC, every row from a copy 1 to 8 bases from the Ref tag, or
  with a 2-base deletion, was labelled `paralog` (450 of 450), and no row of an untouched marker
  was mislabelled (3,348 rows). A copy that is not in the genome is not seen.

See [`filter-paralogs`](../commands/filter-paralogs.md#how-it-was-checked).
