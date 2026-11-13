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

## Genotype Calling (`call`)

`call` adapts polyRAD's model and parts of updog's, so it was checked against both packages and
against data with known or independent genotypes. The tables are in the development record.

- **Same answers as polyRAD.** With polyRAD's settings `call` gives polyRAD's genotypes: on three
  alfalfa MADCs all but 0.01% of 6.8 million target genotype values are identical, and on the polyRAD
  paper's own benchmark (five diversity panels) the error of `call --prior hwe` equals the
  published error of polyRAD's HWE model to the third decimal.
- **Potato against a SNP array.** On the deposit of Endelman et al. (2024), the study's stored
  genotypes scored here reproduce the paper's own per-marker error table. With the bias fitted,
  `call` has the array's dosage for 87.5% of 550,052 genotypes of 298 tetraploid clones, where the
  study's updog genotypes have 87.2%; without it 82.8%, level with polyRAD at the same
  overdispersion and 5 points above polyRAD's defaults.
- **A vendor MADC of mixed ploidy.** The same study's DArT MADC (375 samples, 52 of them diploid)
  was called as delivered: 91.9% of 202,143 genotypes have the array's dosage with the bias
  fitted (the study's genotypes: 91.8%) and 89.1% without.
- **Real diploid reads.** On 96 delta smelt (GT-seq), 99.6% of target genotypes and 99.8% of
  off-target genotypes that pass the default filters agree with genotypes called from alignments
  of the same reads; polyRAD and updog on the same counts agree less often.
- **Polyploids with known dosages from real reads.** Tetraploids and hexaploids pooled from the
  smelt diploids, at full depth and thinned to 50, 25 and 10 reads: `call` is level with updog at
  full depth (98.3% of tetraploid and 97.0% of hexaploid target genotypes right) and ahead of it
  once reads are thinned, by 1 to 4 points for tetraploids and 6 to 10 for hexaploids.
- **Hexaploid replicates.** On ten samples of each of eight sweetpotato clones, with the read
  counts the study's updog run used, two samples of a clone get the same dosage in 86.7% of
  pairs with `call` (norm prior, bias fitted) and 87.3% with updog.
- **Simulated tables with known dosages**, from updog's and polyRAD's own simulators: diploid,
  tetraploid and hexaploid, populations, groups and families, 8 to 150 reads, with and without
  allelic bias.

Limits found: markers on one subgenome of an allopolyploid are outside the model; linkage between
markers is not used, and polyRAD's models that use it are better on mapping populations; fitting
the allelic bias helps where one allele gives fewer reads (potato) and costs accuracy where it does
not, so it is off by default.

See [`call`](../commands/call.md).

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

- **Alfalfa 3K DArTag, two aligners:** on two tetraploid F1 families and the panel's haplotype
  database, with a monoploid genome of another cultivar, the filter kept or dropped the same rows
  as bwa mem for 99.4 to 99.9% and as blastn for 99.3 to 99.7%; bwa mem and blastn agree with each
  other on 99.4 to 99.8%.
- **Alfalfa 3K DArTag, inheritance:** mapped genetically in the two families, 97.6% of the rows the
  filter kept lie at their own marker, and 84 to 95% of the `paralog` rows that the test can
  judge lie somewhere else. The filter removes 4% of the reads of this panel, against 63% in
  lettuce.

See [`filter-paralogs`](../commands/filter-paralogs.md#how-it-was-checked).
