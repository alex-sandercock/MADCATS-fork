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
2019). Within a record each allele is considered on its own, against all the other alleles
together:

1. **Likelihood.** In a sample of ploidy *P* an allele can be present in 0 to *P* copies (its
   dosage). The probability of the allele's read count, out of the sample's reads in the record,
   is beta-binomial. The expected read fraction of a dosage is dosage / *P*, moved slightly toward
   the allele's frequency to allow for contamination between samples
   (`--contamination-rate`). `--overdispersion` sets how far read fractions may stray from the
   expected fraction.
2. **Prior.** `--prior uniform` makes every dosage equally probable. `--prior hwe` (the default)
   uses Hardy-Weinberg proportions for the allele's frequency, estimated over all the samples
   together.
3. **Posterior.** Prior times likelihood, normalised over the dosages.
4. **Iteration** (`hwe` only). The allele frequencies are re-estimated from the posteriors and
   steps 1–3 repeat until the frequencies change by no more than `--tol`, or `--max-iterations` is
   reached. Each record converges on its own.

The Hardy-Weinberg prior assumes random mating and that every chromosome copy is drawn
independently, which holds for diploids and autopolyploids. It is not the expectation under
disomic inheritance in an allopolyploid, in a biparental family, or when samples of different
ploidy come from different populations. Use `--prior uniform` when it does not fit.

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
| `--prior {uniform,hwe}` | `hwe` | Genotype prior (see [the genotype model](#the-genotype-model)) |
| `--selfing-rate <FLOAT>` | 0 | With `--prior hwe`, the share of individuals that come from self-fertilisation: 0 is random mating, 1 is fully inbred lines |
| `--overdispersion <auto\|FLOAT>` | `auto` | Beta-binomial overdispersion parameter, as in polyRAD. Lower values allow read fractions to stray further from the expected fraction. `auto` estimates it from the read counts (see below); a number sets it |
| `--contamination-rate <FLOAT>` | 0.001 | Expected share of a sample's reads that come from other samples |
| `--tol <FLOAT>` | 1e-5 | With `--prior hwe`, stop iterating a record once its allele frequencies change by no more than this on average |
| `--max-iterations <INT>` | 100 | With `--prior hwe`, the most iterations a record may take |

### Choosing the overdispersion

`--overdispersion` matters most in polyploids, where neighbouring dosages differ by a small read
fraction. polyRAD's default of 9 was chosen for RAD-seq. On two tetraploid alfalfa DArTag
populations it left 59% of target genotypes below the confidence thresholds, where the estimated
values (40 and 50) left 8–11%, with the same agreement with published dosages.

By default (`auto`) the value is estimated by maximum likelihood. For each value on a grid from 2
to 500, `call` adds up the log-likelihood of every marker's reference and alternate read counts,
with each sample's genotype summed out under Hardy-Weinberg proportions, and takes the value with
the highest total. It uses up to 20,000 sample-by-marker cells, from markers at even spacing
through the table, and takes a few seconds. The value used, and the log-likelihood of every value
tried, are in the QC JSON; the VCF header states the value and that it was estimated.

- The estimate assumes Hardy-Weinberg proportions. On simulated F1 families it came out up to 1.5
  times the true value. On inbred lines, which have no heterozygotes to measure scatter from, it
  ends at the lowest value tried; `call` then warns and uses 9.
- When it is the highest value tried (500), the counts scatter no more than random sampling of
  reads would.
- A table with fewer than 100 sample-by-marker cells with reads is too small to estimate from;
  `call` says so and uses 9.
- To reproduce a polyRAD run, give its value: `--overdispersion 9`.
- This is not polyRAD's `TestOverdispersion`, which on simulated tetraploid and hexaploid counts
  chose values well below the true one.

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
| `--jobs <INT>` | up to 4 | Worker processes. Markers are independent, and the output does not depend on this |

## Examples

```bash
# Target genotypes of a diploid project
madcats call --madc results/project_MADC.csv --output-prefix results/project

# Target and haplotype VCFs of a mixed-ploidy project
madcats call --madc results/project_MADC.csv --output-prefix results/project \
  --vcf target haplotype --ploidy 4 --ploidy-file ploidy.tsv
```

## Outputs

| File | Contents |
|------|----------|
| `PREFIX.targets.vcf`, `PREFIX.all.vcf`, `PREFIX.haplotypes.vcf` | The VCFs chosen with `--vcf` (VCF 4.3, uncompressed) |
| `PREFIX.haplotypes.alleles.csv` | With `--vcf haplotype`: for each marker and allele number, the MADC AlleleIDs it was counted from, its sequence and its total reads |
| `PREFIX.qc.json` | Parameters, the model in words, samples per ploidy, and per VCF: records, genotypes called and missing by reason, reads used and excluded, records that did not converge, targets skipped by reason |

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
| `INFO/NS` | all | Samples with a called genotype |
| `INFO/AF` | all | Estimated frequency of each alternate allele |

In the haplotype VCF, rows with the same sequence are added together, and an alternate sequence
with no reads in any sample is not an allele. An indel haplotype is written as the tag that was
read, which ends a few bases away from where the Ref tag ends.

In the `all` VCF, the three record kinds come from the same reads, so agreement between them is a
check on consistency, not independent evidence. Off-target insertions and deletions are not
written as records yet.

## See Also

- [madcats generate](generate.md)
- [madcats filter-paralogs](filter-paralogs.md)
- [Output Formats](../reference/output-formats.md)
