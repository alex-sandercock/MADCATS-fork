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
  full depth (98.3% of tetraploid and 96.9% of hexaploid target genotypes right) and ahead of it
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

### Priors

**`norm` in families.** It helped in simulations and at low depth, and not on real families at
high depth:

- Simulated tetraploid F1 families at 20 to 60 reads, five draws: 91.0% of dosages right against
  90.1% with `hwe`, and ahead in every draw (0.6 to 1.5 points). On a biased panel it was 2.5
  points ahead with `--allele-bias auto` and 4.4 with `fit`; on updog's simulator at 150 reads,
  7.6 and 10.2 (91.5% and 91.6% against 84.0% and 81.4%).
- With `norm`, `auto` kept no bias on F1 families without one; with `hwe` it kept a bias at 2
  markers in 100 and lost 0.7 points.
- The two F1 populations of the polyRAD paper's benchmark, at about 10 reads: error 0.186 against
  0.200 (*Miscanthus*, diploid) and 0.385 against 0.409 (potato, tetraploid), level with polyRAD's
  mapping prior, which is told the parents.
- Two tetraploid alfalfa F1 populations against published dosages: 90.1% and 89.8% with `norm`,
  90.9% and 90.0% with `hwe`. With the bias fitted, 92.3% and 91.9% against 92.2% and 91.8%.
- 85 potato clones of ten families against a SNP array: 92.1% against 92.3%, and 96.57% against
  96.60% with the bias fitted.
- In one randomly mating population (five draws) it was level with `hwe` where the data have no
  bias, and 0.5 to 0.7 points below where they have one and it is not fitted or `auto` is used.
- On simulated inbred lines without a selfing rate it was 1.1 points ahead of `hwe`.

**`structure`.** It matters at low read depth:

- On simulated tetraploids with 8 reads per sample and marker it called 0.8 to 5.5 points more
  genotypes correctly than `hwe`, depending on how strong the structure was; with 25 reads 0.5 to
  2.5 points, and with 150 reads 0.8.
- On 85 potato clones from ten families, thinned to 10 and 25 reads, the gain against SNP-array
  genotypes was 0.9 and 0.7 points; at their full depth of about 300 reads there was none.
- Compared with polyRAD's `IteratePopStruct` on the same thirteen tables it called the same share
  of genotypes correctly or more, by up to 0.7 points. On a table without structure polyRAD's
  prior lost 3.4 points and `call` chose 0 neighbours.
- polyRAD models structure with principal components; that approach and a relationship-matrix
  predictor were tried in place of the neighbours and were no better (see the development
  record).

**The selfing rate.** On the inbred soybean set of the polyRAD benchmark the error was 0.037 with
`--selfing-rate 0.95`, as polyRAD has, and 0.157 without. On simulated inbred tetraploids 99.6%
and 99.9% of dosages were right at 0.95 and 0.99, against 97.4% without. On simulated half-inbred
tetraploids (four draws) 93.9% were right with `--selfing-rate 0.5`, against 93.5% without and
94.1% with `--prior norm`.

**Haplotype records.** On simulated microhaplotype tables (3 to 6 alleles per marker, diploid and
tetraploid; one population, groups and families):

- `hwe` gave the same share of correct genotypes as polyRAD's HWE prior, to within 0.4 points.
- `structure` was never below `hwe` and gained 1 to 4 points in tetraploids with groups or
  families.
- `norm` was never below `hwe`, gained at most 1.7 points, and took two to five times as long.
- `uniform` was the lowest throughout, by up to 6 points.
- With the default `--min-marginal-posterior`, 94 to 99% of the genotypes kept were right under
  `hwe`, `norm` and `structure` alike.

These are simulations with one simulator; there is no published accuracy baseline for haplotype
calls to compare with, and no real haplotype data with known genotypes has been tested.

**Off-target records.** On simulated microhaplotype tables (diploid and tetraploid, 8 to 150
reads) off-target records were called as accurately as target records, whether a substitution was
on one haplotype or shared by two, and 96 to 99.8% of the off-target genotypes that passed the
default filters were right. On real reads of 96 diploid delta smelt, 99.8% of the off-target
genotypes that passed the default filters agreed with genotypes called from alignments of the
same reads, against 99.6% at the targets.

### Allelic Bias

Share of dosages right with `--prior hwe`, with the rules of `auto` as released. Real data:

| Data | `off` | `fit` | `auto` | Biases `auto` used |
|------|-------|-------|--------|--------------------|
| Potato, 85 tetraploid clones, against a SNP array | 92.3% | 96.6% | 96.4% | 785 of 1,544 |
| Potato vendor MADC, tetraploids | 86.4% | 91.3% | 91.1% | 2,566 of 3,904, with the diploids |
| Potato vendor MADC, diploids | 92.0% | 92.6% | 92.5% | |
| Delta smelt (diploid), against alignment-based genotypes | 99.4% | 99.3% | 99.3% | 127 of 946 |
| Tetraploids pooled from real diploid reads, 50 reads, three draws | 90.4% | 89.5% | 90.4% | 15 to 29 of 946 |
| Hexaploids pooled from real diploid reads, 50 reads, three draws | 77.8% | 71.0% | 77.8% | 3 to 12 of 944 |

Simulated tables, as the mean over five draws (four for updog's simulator):

| Scenario | `off` | `fit` | `auto` | Biases `auto` used |
|----------|-------|-------|--------|--------------------|
| Tetraploids, random mating, no bias, 20 to 60 reads | 89.7% | 89.2% | 89.7% | none |
| Tetraploids, random mating, biased, 20 to 60 reads | 83.9% | 88.7% | 87.0% | 24% |
| Hexaploids, random mating, no bias, 20 to 60 reads | 75.0% | 73.3% | 75.0% | none |
| Hexaploids, random mating, biased, 20 to 60 reads | 66.7% | 72.6% | 67.8% | 7% |
| Tetraploid inbred lines, no bias | 97.3% | 73.7% | 97.3% | none |
| updog's simulator: tetraploids, biased, 150 reads | 84.8% | 96.5% | 96.3% | 72% |
| updog's simulator: tetraploids, biased, 25 reads | 76.1% | 84.7% | 82.8% | 35% |
| updog's simulator: hexaploids, biased, 150 reads | 66.7% | 85.8% | 82.6% | 52% |
| updog's simulator: diploids, biased, 150 reads | 99.6% | 100.0% | 100.0% | 75% |

- **Where a bias is there, fitting it is worth 4 to 5 points in tetraploid potato** and up to 19
  in simulated hexaploids. It brought `call` level with updog on biased polyploid data (within 0.7
  points on seven tables from updog's simulator), where calling without it was 8 to 20 points
  behind. Against the potato array, at the markers updog found biased the share rose from 80% to
  95%. On the whole deposit of that study (703 clones, 1,862 markers, none selected) `fit` raised
  the share from 82.8% to 87.5%, where the study's own updog genotypes have 87.2%.
- **Where no bias is there, `auto` leaves the calls alone.** In one population and in inbred
  lines it kept no bias in any draw, and its genotypes are those of `off`. `fit` lost 0.1 to 1.8
  points in one population, 7 on hexaploids pooled from real reads and 24 on inbred lines.
- **`auto` needs depth.** Of what `fit` gains it kept 97 to 100% at 150 reads in diploids and
  tetraploids and 83% in hexaploids, 78% at 25 reads, and at 20 to 60 reads with smaller biases
  54 to 68% in tetraploids and 19% in hexaploids.
- **Inbred lines.** A bias changes only what a heterozygote shows. Under a prior that expects
  heterozygotes, every fit ran to a limit of the search and homozygous samples passed for
  heterozygous ones. `auto` rejects a bias at a limit; `fit` does not.
- **Families, subpopulations and partly inbred samples without bias.** The plain `hwe` prior
  does not fit them, and a bias can make up for it: `auto` kept a bias at 2, 4 and 12 markers in
  100 and lost 0.7, 0.6 and 0.8 points. With `--prior norm` it kept none in F1 families and lost
  0.1 points in subpopulations and in half-inbred samples; with `--selfing-rate 0.5` it lost 0.2
  in half-inbred samples. None of these biases ends at a search limit.
- **Mixed ploidy.** With diploids, tetraploids and hexaploids in one table (polyRAD's simulator,
  no bias), `auto` was within 0.1 points of `off` in every draw.
- **Real diploid reads.** On 96 delta smelt `auto` and `fit` were 0.1 points under `off` at
  targets and 0.3 at off-target SNPs.

### Calibration

How far the probability `call` gives its called genotypes lies above the share of them that are
right, in points. Positive is over-confident. The probability is read from `GQ`; simulated
figures are means over the draws.

| Data | `off` | `auto` | `fit` |
|------|-------|--------|-------|
| Simulated tetraploids, random mating, no bias | −0.1 | −0.1 | +0.6 |
| Simulated tetraploids, random mating, biased | +1.9 | −0.3 | +1.1 |
| Simulated hexaploids, random mating, no bias | −0.5 | −0.5 | +4.0 |
| Simulated hexaploids, random mating, biased | +4.6 | +3.7 | +5.0 |
| Simulated tetraploid inbred lines | +0.4 | +0.4 | +14.4 |
| updog's simulator: tetraploids, biased, 150 reads | +5.4 | −0.5 | +0.5 |
| updog's simulator: hexaploids, biased, 150 reads | +10.8 | 0.0 | +2.2 |
| updog's simulator: tetraploid F1, biased, `--prior hwe` | +8.6 | +11.1 | +14.7 |
| Tetraploids pooled from real reads, three draws | +1.1 | +1.2 | +2.3 |
| Hexaploids pooled from real reads, three draws | +0.8 | +0.8 | +8.3 |
| Potato clones against the array | +4.8 | +2.3 | +2.5 |

- Without a bias in the data, `off` and `auto` are within 1.2 points. `fit` is over-confident
  by 4 to 8 points in hexaploids and 14 on inbred lines.
- With a bias in the data, `off` is over-confident by up to 11 points, because a biased marker is
  called wrongly and confidently. `auto` and `fit` are within 2.2 points at 150 reads.
- A prior that does not fit the samples costs calibration whatever the bias setting: F1 families
  under `--prior hwe` are over-confident by 9 to 15 points.
- Against the potato array the gap is not all the caller's: the array and the reads also differ
  for reasons a caller cannot see. A genotype with a `GQ` of 30 or more was wrong 6 to 7 times in
  100 there.

### Overdispersion and Background Reads

- In simulated diploids and tetraploids with 20 or more reads per sample both values were
  recovered (31 of 32 tables). In hexaploids the overdispersion was recovered in 12 of 16 tables
  and read between 14 and 50 for a simulated 40 in the rest. With 6 to 25 reads per sample the
  background rate was still recovered, an overdispersion of 40 read as 30 to 75, and one of 150
  as anything from 100 to 500: scatter that small cannot be measured at that depth.
- The estimate gives every dosage its own free probability at each marker. With Hardy-Weinberg
  proportions it was badly off for samples from two subpopulations, and with those or the `norm`
  shape, neither of which can describe an excess of homozygotes, it failed for inbred samples.
- polyRAD's `TestOverdispersion` chose values well below the true one on simulated tetraploid and
  hexaploid counts. On 85 tetraploid potato clones it chose 14 where `call` estimates 75, and
  polyRAD's HWE model with 14 matched a SNP array for 91.2% of dosages, against 92.3% for `call`
  with its defaults and 87.9% for polyRAD's own default of 9.

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

- **Alfalfa 3K DArTag, design genome:** on two tetraploid F1 families and the panel's haplotype
  database, with chromosome set 1 of XinJiangDaYe as the genome, every marker was placed at its
  CloneID position, and the filter kept or dropped the same rows as bwa mem for 99.7 to 99.9% and
  as blastn for 99.4 to 99.7%.
- **Alfalfa 3K DArTag, inheritance:** mapped genetically in the two families, 96% of the rows the
  filter kept lie at their own marker, and 91 to 94% of the `paralog` rows that the test can
  judge lie somewhere else. The filter removes 3% of the reads of this panel, against 63% in
  lettuce, and the target genotypes are the same with and without it.
- **Alfalfa 3K DArTag, more than one chromosome set:** with a second set of the design assembly
  in the genome the filter dropped 1,969 of 2,753 RefMatch and AltMatch rows where it dropped 296
  with one. A haplotype-resolved or polyploid assembly has to be cut down to one set first.

See [`filter-paralogs`](../commands/filter-paralogs.md#how-it-was-checked).
