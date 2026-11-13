# madcats call

Call genotypes from a MADC and write them as VCFs.

## Synopsis

```bash
madcats call --madc <FILE> --output-prefix <PREFIX> [--vcf {target,all,haplotype} ...] [OPTIONS]
```

## Description

`call` estimates each sample's genotype from the read counts in a MADC. It needs no other tool:
the MADC can come straight from `generate` or `filter-paralogs`, from a vendor delivery (with or
without the starred preamble), or from HapApp (the fixed-allele-ID layout). Ploidy can differ
between samples.

`--vcf` chooses what a record is. Each type is written to its own file, and one run can write
several:

| `--vcf` | File | One record per | Alleles |
|---------|------|----------------|---------|
| `target` (default) | `PREFIX.targets.vcf` | marker | The designed variant. Ref and RefMatch reads count for the reference allele, Alt and AltMatch reads for the alternate allele; reads of Other rows are left out |
| `all` | `PREFIX.all.vcf` | variant | The target records, plus one record for every other tag position where a row with reads shows a base substitution |
| `haplotype` | `PREFIX.haplotypes.vcf` | marker | Whole tag sequences (microhaplotypes): the Ref tag, then every other sequence with reads |

Coordinates are **marker-local**. `CHROM` is the marker (its CloneID) and `POS` is a 1-based
position on its Ref tag, so these files are not placed on a genome. Every marker is declared as a
contig in the header, and a `##madcatsCoordinates` line says so.

### The genotype model

The model is adapted from [polyRAD](https://github.com/lvclark/polyRAD) (Clark, Lipka and Sacks
2019), with the `norm` prior of [updog](https://github.com/dcgerard/updog) (Gerard and Ferrão
2020) and a prior for structured populations as options. Within a record each allele is considered
on its own, against all the other alleles together:

1. **Likelihood.** In a sample of ploidy *P* an allele can be present in 0 to *P* copies (its
   dosage). The probability of the allele's read count, out of the sample's reads in the record,
   is beta-binomial. The expected read fraction of a dosage is dosage / *P*, changed by the
   allelic bias if one is fitted (`--allele-bias`) and then by the share of reads that show the
   wrong allele (`--background-rate`). `--overdispersion` sets how far read fractions may stray
   from the expected fraction.
2. **Prior.** `--prior uniform` makes every dosage equally probable. `--prior hwe` (the default)
   uses Hardy-Weinberg proportions for the allele's frequency, estimated over all the samples
   together. `--prior norm` gives the dosages a normal shape, with a mean and a spread fitted to
   the samples. `--prior structure` uses Hardy-Weinberg proportions at each sample's own allele
   frequency, taken from the samples most like it.
3. **Posterior.** Prior times likelihood, normalised over the dosages.
4. **Iteration** (`hwe`, `norm` and `structure`). The allele frequencies, and with `norm` the
   prior's mean and spread, are re-estimated from the posteriors and steps 1–3 repeat until they
   change by no more than `--tol`, or `--max-iterations` is reached. Each record converges on its
   own.

The Hardy-Weinberg prior assumes random mating and that every chromosome copy is drawn
independently, which holds for diploids and autopolyploids. It is not the expectation under
disomic inheritance in an allopolyploid, in a biparental family, or when samples of different
ploidy come from different populations. Use `--prior norm` or `--prior uniform` when it does not
fit, or `--prior structure` when the samples come from several subpopulations or families and the
read depth is low.

The norm prior makes no assumption about mating. Its dosage probabilities are proportional to a
normal density at the dosages 0 to *P*, so it can be narrower than Hardy-Weinberg proportions (one
family) or wider (subpopulations), and when the dosages have two peaks (inbred lines) it flattens
out. It is updog's `norm` genotype distribution, fitted the way updog fits it; with the same
likelihood, `call` and updog reach the same prior and the same genotypes. It is fitted separately
for each ploidy among the samples, since a shape over the dosages of one ploidy says nothing about
another. A ploidy with fewer than 10 samples uses Hardy-Weinberg proportions instead; `call` says
so and lists those ploidies in the QC JSON. In a haplotype record every allele gets its own fit.

The structure prior is for samples that are not one randomly mating population: subpopulations, a
cline, or a set of families. It keeps Hardy-Weinberg proportions but gives each sample its own
allele frequency in every record: the mean over the samples most like it (its neighbours), with
the record's frequency over all samples counted as one more of them. A sample's own reads never
enter its own prior. The neighbours are found once for the whole table, before any record is
called:

1. The target counts of the markers (at most 1,000, at even spacing) are called with the `hwe`
   prior, which gives every sample an expected share of alternate copies at each marker.
2. Two samples are alike when those shares differ little over the markers. This ranks pairs of
   samples as a genomic relationship matrix would.
3. The number of neighbours is chosen from 0, 5, 10, 20, 40 and 80 by prediction: the neighbours
   found on one half of the markers predict each sample's shares at the other half, and the number
   with the smallest error is used. If none predicts better than the table's own frequencies, the
   number is 0 and the prior is the `hwe` prior; `call` says so.

Records are then called one at a time as with the other priors, each first with the `hwe` prior
and then with the neighbours. The number of neighbours used and the error of each number tried are
in the QC JSON.

- It matters at low read depth. On simulated tetraploids with 8 reads per sample and marker it
  called 0.8 to 5.5 points more genotypes correctly than `hwe`, depending on how strong the
  structure was; with 25 reads 0.5 to 2.5 points, and with 150 reads 0.8. On 85 potato clones from
  ten families, thinned to 10 and 25 reads, the gain against SNP-array genotypes was 0.9 and 0.7
  points; at their full depth of about 300 reads there was none.
- Compared with polyRAD's `IteratePopStruct` on the same thirteen tables it called the same share
  of genotypes correctly or more, by up to 0.7 points. On a table without structure polyRAD's
  prior lost 3.4 points and `call` chose 0 neighbours.
- It needs at least 20 samples and 100 markers with reads to compare samples on; below that it is
  the `hwe` prior, with a note.
- It cannot be combined with `--allele-bias fit` or `--selfing-rate`. Neither combination has been
  tested.
- It takes two to three times as long as `--prior hwe`.
- The method is MADCATS's own, not a port. polyRAD models structure with principal components;
  that approach and a relationship-matrix predictor were tried in its place and were no better
  (see the development record).

All four priors apply to every VCF type. For the haplotype VCF, on simulated microhaplotype tables
(3 to 6 alleles per marker, diploid and tetraploid; one population, groups and families):

- `hwe` gave the same share of correct genotypes as polyRAD's HWE prior, to within 0.4 points.
- `structure` was never below `hwe` and gained 1 to 4 points in tetraploids with groups or
  families.
- `norm` was never below `hwe`, gained at most 1.7 points, and took two to five times as long.
- `uniform` was the lowest throughout, by up to 6 points.
- With the default `--min-marginal-posterior`, 94 to 99% of the genotypes kept were right under
  `hwe`, `norm` and `structure` alike.

These are simulations with one simulator; there is no published accuracy baseline for haplotype
calls to compare with, and no real haplotype data with known genotypes has been tested. The
allelic bias is never fitted for haplotype records.

For inbred lines and partly selfing species, `--selfing-rate` adjusts the Hardy-Weinberg prior: a
share of each generation comes from self-fertilisation, which moves prior probability from the
intermediate dosages to the two homozygous ones (the equilibrium proportions of de Silva et al.
2005, as in polyRAD). At 1 only the two homozygous dosages have prior probability, so a
heterozygote cannot be called; use a value just below 1, such as 0.99, to allow residual
heterozygosity. It applies to samples of even ploidy; samples of odd ploidy keep random-mating
proportions.

### What the model can and cannot say

With **two alleles**, the dosage of the alternate allele fixes the genotype, so its posterior is
the genotype posterior. These records carry `GQ`, `GP` and `PL`. `PL` is computed from the
likelihoods and `GP` from the posteriors, so they differ by the prior.

With **more than two alleles**, there is one posterior per allele, not one over whole genotypes.
The most probable dosages of the alleles need not add up to the ploidy; when they do not, copies
are added or removed one at a time where that costs the least probability (polyRAD's method).
These records carry `MPC`: the lowest posterior among the called dosages. It is an upper bound on
the probability that the whole genotype is right, not an estimate of it, which is why it is not
called `GQ`. Every haplotype record is written this way, and so is an `all` record with three or
more bases.

### Missing genotypes

A genotype is written as missing (`./.`, with one dot per chromosome copy) when:

- the sample has fewer than `--min-dp` reads in the record;
- two dosages are equally probable;
- its `GQ` is below `--min-gq` (two-allele records) or its `MPC` is below
  `--min-marginal-posterior` (other records);
- more than `--max-other-fraction` of its reads were left out of the call (when that option is
  given).

`AD` and `DP` are always kept. A sample with reads but below `--min-dp` keeps its probabilities
too; a sample with no reads has `.` for them, since they would only repeat the prior.

### Finding the target

The target record needs the designed variant's position and alleles. `call` takes them from the
marker's Ref and Alt tags (the Ref and Alt rows of the MADC, or `--reference-tags` when given):

- If the two tags differ by exactly one substitution, insertion or deletion, that is the target.
  A designed indel is written the VCF way, from the base before it.
- If they differ in more places, the target is the substitution at the reference tags'
  `VariantPositionInTagSequence`. Without reference tags the marker gets no target record.
- A marker without an Alt row gets no target record.

The target is never guessed. Skipped markers are counted by reason in the QC JSON, and they still
get haplotype records and off-target records.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--madc <FILE>` | MADC to call from: from `generate` or `filter-paralogs`, a vendor DArT MADC, or a HapApp fixed-allele-ID MADC. Gzipped input is read directly |
| `--output-prefix <PREFIX>` | Path prefix of the outputs |

## Output Selection

| Option | Default | Description |
|--------|---------|-------------|
| `--vcf {target,all,haplotype}` | `target` | VCF types to write, one file each. Give one or more: `--vcf target haplotype` |
| `--qc-output <FILE>` | `PREFIX.qc.json` | QC JSON |

## Samples and Panel

| Option | Default | Description |
|--------|---------|-------------|
| `--ploidy <INT>` | 2 | Ploidy of every sample not listed in `--ploidy-file` |
| `--ploidy-file <FILE>` | none | Two columns, sample name and ploidy, separated by a tab or a comma, with or without a header line. Names are the MADC's sample column names |
| `--reference-tags <CSV>` | none | The panel's reference-tag CSV (from `design` or `infer-design`). Its tags and `VariantPositionInTagSequence` identify each marker's target |
| `--max-indel-length <INT>` | 10 | Widest insertion or deletion considered when a row is compared with the Ref tag. Use the value given to `generate` |

## Model

| Option | Default | Description |
|--------|---------|-------------|
| `--prior {uniform,hwe,norm,structure}` | `hwe` | Genotype prior (see [the genotype model](#the-genotype-model)) |
| `--selfing-rate <FLOAT>` | 0 | With `--prior hwe`, the share of individuals that come from self-fertilisation: 0 is random mating, 1 is fully inbred lines |
| `--allele-bias {off,fit}` | `off` | Allelic bias, where a copy of one allele yields more reads than a copy of the other. `fit` estimates it for each record with two alleles (see below). Not with `--prior structure` |
| `--overdispersion <auto\|FLOAT>` | `auto` | Beta-binomial overdispersion parameter, as in polyRAD. Lower values allow read fractions to stray further from the expected fraction. `auto` estimates it from the read counts (see below); a number sets it |
| `--background-rate <auto\|FLOAT>` | `auto` | Probability that a read shows one allele although it came from another. `auto` estimates it from the read counts (see below); a number sets it |
| `--contamination-rate <FLOAT>` | 0 | polyRAD's contamination: the share of a sample's reads that come from other samples and show an allele at its frequency among them. `--background-rate` covers the same reads and more, so this is off unless set; polyRAD uses 0.001 |
| `--tol <FLOAT>` | 1e-5 | With `--prior hwe`, `norm` or `structure`, stop iterating a record once its allele frequencies change by no more than this (with `norm`, also its prior probabilities; with `structure`, also the samples' own allele frequencies) |
| `--max-iterations <INT>` | 100 | With `--prior hwe`, `norm` or `structure`, the most iterations a record may take; with `structure`, in each of its two stages |

### Overdispersion and background reads

Two numbers describe how read counts stray from what a genotype predicts, and both are estimated
from the table unless you set them.

- **Overdispersion** is how far read fractions scatter around the expected fraction. It matters
  most in polyploids, where neighbouring dosages differ by a small read fraction. polyRAD's
  default of 9 was chosen for RAD-seq; on two tetraploid alfalfa DArTag populations it left 59% of
  target genotypes below the confidence thresholds, against 8–11% with estimated values.
- **The background-read rate** is the share of reads that show an allele they did not come from:
  sequencing and amplification errors, index hopping and contamination together, as they remain
  after the MADC was built. If it is set too low, such reads can only be explained as scatter, and
  the overdispersion estimate drops with it.

`call` chooses the pair that makes the read counts most probable. It tries background rates from
0.0002 to 0.02 and overdispersions from 2 to 500, on each marker's reference and alternate read
counts with the genotypes summed out, using up to 20,000 sample-by-marker cells from markers at
even spacing through the table. The two can be told apart because they show in different samples:
stray reads in samples that do not carry an allele, scatter in samples that carry both. The values
used, and the log-likelihood of every value tried, are in the QC JSON, and the VCF header states
them.

- The estimate does not use `--prior`. It gives every dosage its own free probability at each
  marker, because the estimate is only as good as the genotype proportions it assumes. With
  Hardy-Weinberg proportions it was badly off for samples from two subpopulations, and with those
  or the `norm` shape, neither of which can describe an excess of homozygotes, it failed for
  inbred samples.
- In simulated diploids and tetraploids with 20 or more reads per sample both values were
  recovered (31 of 32 tables). In hexaploids the overdispersion was recovered in 12 of 16 tables
  and read between 14 and 50 for a simulated 40 in the rest. With 6 to 25 reads per sample the
  background rate was still recovered, an overdispersion of 40 read as 30 to 75, and one of 150
  as anything from 100 to 500: scatter that small cannot be measured at that depth.
- A table with no heterozygous samples and no stray reads has no scatter to measure: the
  overdispersion estimate ends at the lowest value tried, and `call` warns and uses 9. Inbred
  lines whose reads do include stray reads of the other allele can be estimated from those.
- An overdispersion of 500, the highest value tried, means the counts scatter no more than random
  sampling of reads would.
- A table with fewer than 100 sample-by-marker cells with reads is too small to estimate from;
  `call` says so and uses an overdispersion of 9 and a background rate of 0.002.
- With `--allele-bias fit` the estimate fits the bias too. Without it, biased markers look like
  scatter and the overdispersion comes out lower (30 and 40 against 150 on the two alfalfa
  populations).
- This is not polyRAD's `TestOverdispersion`, which on simulated tetraploid and hexaploid counts
  chose values well below the true one.

To run polyRAD's model as polyRAD runs it:
`--prior hwe --overdispersion 9 --contamination-rate 0.001 --background-rate 0`.

### Allelic bias

In amplicon data a copy of one allele often yields more reads than a copy of the other, so a
sample with two copies of each allele in four does not show half its reads for each. Called
without regard to this, such a marker gets shifted dosages. `--allele-bias fit` estimates the bias
of every record with two alleles (targets, and off-target sites with two bases) and calls with it.
Haplotype records and sites with three or more bases are called without a bias.

The bias is updog's (Gerard et al. 2018): with bias *h*, an expected read fraction *f* of the
reference allele becomes *f* / (*h*(1 − *f*) + *f*), so below 1 the reference allele is
over-represented and above 1 under-represented. Like updog, `call` pulls the log of the bias
towards 0 with a normal prior of standard deviation 0.7, starts from updog's five starting values
(0.37 to 2.7) and keeps the fit with the highest log-likelihood after that prior. With the other
settings equal, the fitted bias and the genotypes match updog's.

- The bias is found by calling the record at a series of bias values, 17 or more. On 1,000 alfalfa
  markers and 186 samples, calling took eight times as long on one worker (57 s against 7 s), and
  the estimate of the overdispersion and background rate, which fits the bias for every pair of
  values it tries, about 115 s against 17 s. The whole run, 170 s, was still 6.5 times as fast as
  updog. The estimate does not grow with the table, and giving `--overdispersion` and
  `--background-rate` as numbers skips it.
- It is applied to the sample's own allele copies, and background reads are added afterwards.
  updog applies its sequencing error first. In this order, samples that do not carry an allele say
  nothing about its bias, which matters because the background rate is one value for the table.
- The search covers biases from 0.14 to 7.4. A record whose fit ends at either limit is counted in
  the QC JSON and marked in the bias table; such a marker is usually not behaving as one locus.
- The VCF carries the bias as `INFO/BIAS`, and `PREFIX.bias.csv` lists every fitted record.

When to fit it, from read counts drawn with updog's and polyRAD's own simulators, two tetraploid
alfalfa populations and 85 tetraploid potato clones that also have SNP-array genotypes (the tables
are in the development record):

- Use it with `--prior norm`. With `--prior hwe` it lost accuracy wherever the samples were not in
  Hardy-Weinberg proportions (families, inbred lines).
- On biased polyploid data it brought `call` level with updog (within 0.7 points of genotypes
  correct on seven updog-simulated tables), where calling without it was 8 to 20 points behind. On
  the alfalfa populations it raised agreement with published dosages from 90% to 92%.
- Against the potato array genotypes, which do not come from a read-count caller, it raised the
  share of matching dosages from 92.3% to 96.6%, where updog reached 96.7%. At the markers updog
  found biased the share rose from 80% to 95%.
- On outbred tetraploid data without bias it cost up to 1.5 points.
- Do not fit it for inbred lines. Without heterozygous samples a bias cannot be told from a shift
  of every dosage, and 8 points were lost.
- Check the result in hexaploids with widely scattered counts: on one such simulated table 10
  points were lost, as updog lost them.

## Filters

| Option | Default | Description |
|--------|---------|-------------|
| `--min-dp <INT>` | 10 | Write a genotype as missing when the sample has fewer reads than this in the record |
| `--min-gq <INT>` | 10 | Two-allele records: write a genotype as missing when its `GQ` is below this |
| `--min-marginal-posterior <FLOAT>` | 0.9 | Haplotype records and other records with more than two alleles: write a genotype as missing when its `MPC` is below this |
| `--max-other-fraction <FLOAT>` | none | `target` and `all` records: write a genotype as missing when more than this share of the sample's reads was left out of the call |

## Performance

| Option | Default | Description |
|--------|---------|-------------|
| `--jobs <INT>` | up to 4 | Worker processes. Records are called independently of each other, and the output does not depend on this |

## Examples

```bash
# Target genotypes of a diploid project
madcats call --madc results/project_MADC.csv --output-prefix results/project

# Target and haplotype VCFs of a mixed-ploidy project
madcats call --madc results/project_MADC.csv --output-prefix results/project \
  --vcf target haplotype --ploidy 4 --ploidy-file ploidy.tsv

# A tetraploid diversity panel at low read depth: each sample's prior follows the samples most like it
madcats call --madc results/panel_MADC.csv --output-prefix results/panel \
  --ploidy 4 --prior structure
```

## Outputs

| File | Contents |
|------|----------|
| `PREFIX.targets.vcf`, `PREFIX.all.vcf`, `PREFIX.haplotypes.vcf` | The VCFs chosen with `--vcf` (VCF 4.3, uncompressed) |
| `PREFIX.haplotypes.alleles.csv` | With `--vcf haplotype`: for each marker and allele number, the MADC AlleleIDs it was counted from, its sequence and its total reads |
| `PREFIX.bias.csv` | With `--allele-bias fit`: for each fitted record, the bias of the reference allele, the log-likelihood of the fit (after the prior on the bias), the starting value the search continued from, the number of bias values tried, and whether the fit converged or ended at a search limit |
| `PREFIX.qc.json` | Parameters, the model in words, the overdispersion and background-rate estimates, the number of neighbours chosen for the structure prior, a summary of the fitted biases, samples per ploidy, and per VCF: records, genotypes called and missing by reason, reads used and excluded, records that did not converge, targets skipped by reason |

### VCF fields

| Field | In | Meaning |
|-------|----|---------|
| `GT` | all | Genotype. Allele copies are unphased, so they are joined by `/` |
| `AD` | all | Reads of each allele |
| `DP` | all | Reads used for the call: the sum of `AD` |
| `XD` | `target`, `all` | Reads left out of `AD` and `DP`: Other rows at a target; elsewhere, rows that do not show the position (a deletion or an insertion removed it from the tag) |
| `GQ` | two-allele records | Phred-scaled posterior probability that `GT` is wrong, capped at 99 |
| `DS` | all | Posterior mean number of copies of each alternate allele |
| `GP` | two-allele records | Posterior probability of each genotype, from 0 to *P* alternate copies |
| `PL` | two-allele records | Phred-scaled likelihood of each genotype, before the prior |
| `MPC` | records with more than two alleles; every haplotype record | Lowest marginal posterior among the called dosages |
| `INFO/TARGET` | `target`, `all` | The record is a marker's designed variant |
| `INFO/BIAS` | two-allele records, with `--allele-bias fit` | Fitted allelic bias of the reference allele, as updog reports it for reference read counts |
| `INFO/NS` | all | Samples with a called genotype |
| `INFO/AF` | all | Estimated frequency of each alternate allele |

In the haplotype VCF, rows with the same sequence are added together, and an alternate sequence
with no reads in any sample is not an allele. An indel haplotype is written as the tag that was
read, which ends a few bases away from where the Ref tag ends.

In the `all` VCF, the three record kinds come from the same reads, so agreement between them is a
check on consistency, not independent evidence. Off-target insertions and deletions are not
written as records yet. On simulated microhaplotype tables (diploid and tetraploid, 8 to 150
reads) off-target records were called as accurately as target records, and 96 to 99.8% of the
off-target genotypes that passed the default filters were right; real data with known off-target
genotypes has not been tested.

## See Also

- [madcats generate](generate.md)
- [madcats filter-paralogs](filter-paralogs.md)
- [Output Formats](../reference/output-formats.md)
