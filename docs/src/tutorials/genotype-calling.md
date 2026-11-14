# Calling Genotypes from a MADC

A MADC holds read counts. Most analyses need genotypes: how many copies of each allele a sample
carries. In a diploid the counts nearly say so by themselves. In a tetraploid or a hexaploid,
neighbouring dosages differ by a small share of the reads, and read counts scatter.
[`madcats call`](../commands/call.md) estimates each sample's genotype from the counts and writes
VCFs. It reads the MADC as it is, from `generate`, `filter-paralogs`, the vendor or HapApp, and
needs no other tool.

The steps below use the [example project](../getting-started/example.md), which is diploid. They
are the same for polyploids.

## 1. Call the Targets

```bash
madcats call --madc results/demo_MADC.csv --output-prefix results/demo
```

```text
note: too few read counts to estimate the read model (96 sample-by-marker cells with reads; 100 needed).
Using background-read rate: 0.002
Using overdispersion: 9
results/demo.targets.vcf: 8 records, 96 of 96 genotypes called (100.0%)
Wrote QC to results/demo.qc.json
```

Eight markers, twelve samples, and every genotype called. All 96 equal the simulated genotypes in
`truth_genotypes.csv`.

The note is about the size of the example. On a real table `call` measures two things from the
counts themselves: how far read fractions scatter (the overdispersion) and how often a read shows
an allele it did not come from (the background-read rate). Eight markers are too few for that, so
it uses fixed values and says so.

## 2. Read the VCF

One record of `results/demo.targets.vcf`, with two of its twelve samples:

```text
#CHROM           POS  REF  ALT  INFO                 FORMAT                    Line_01                            Line_02
Chr01_000412883  34   G    C    TARGET;NS=12;AF=0.5  GT:AD:DP:XD:GQ:DS:GP:PL   0/0:76,0:76:0:48:0:1,0,0:0,51,142  0/1:35,35:70:0:40:1:0.0001,0.9999,0.0001:40,0,40
```

`Line_01` has 76 reads of G and none of C, so it is `0/0`. `Line_02` has 35 of each and is `0/1`.

| Field | Meaning |
|-------|---------|
| `GT` | The genotype, with one allele per chromosome copy: `0/0/0/1` in a tetraploid |
| `AD`, `DP` | Reads of each allele, and their sum |
| `GQ` | The model's confidence in `GT`, Phred-scaled: 20 is one chance in 100 of a wrong genotype, 30 one in 1,000 |
| `DS` | The expected number of alternate copies, for analyses that take a dosage |
| `GP` | The probability of each genotype |

**Coordinates are marker-local.** `CHROM` is the marker and `POS` is a position in its tag, not on
a genome. The header says so in a `##madcatsCoordinates` line.

A genotype is written as missing (`./.`) when the sample has fewer than 10 reads or the `GQ` is
below 10. `--min-dp` and `--min-gq` change that. The read counts stay in the record either way.

## 3. Give the Ploidy

`call` takes every sample to be diploid unless told otherwise:

```bash
madcats call --madc project_MADC.csv --ploidy 4 --output-prefix results/project
```

When ploidy differs between samples, list the exceptions in a two-column table and pass it with
`--ploidy-file`. Samples it does not list take `--ploidy`:

```text
sample	ploidy
Clone_017	2
Clone_044	6
```

## 4. Choose the Settings for Your Samples

Start with the defaults for diploids, and add `--allele-bias auto` for polyploids:

```bash
madcats call --madc project_MADC.csv --ploidy 4 --allele-bias auto --output-prefix results/project
```

`call` estimates how the read counts scatter from the table itself. With `auto` it also fits an
allelic bias at every marker, and uses it only where the bias is clear. Markers left without a
bias are called exactly as without the option.

Change this when one of these describes your samples:

| Your samples | Change | Why |
|--------------|--------|-----|
| Inbred lines, doubled haploids, selfing species | `--selfing-rate`: `0.95` to `0.99` for inbred lines, the species' own rate for partly selfing ones. Leave the bias off for inbred lines | The default prior expects the heterozygotes of random mating, and a bias has nothing to correct where there are none |
| Families (F1, half-sibs) | Add `--prior norm` | The default prior assumes random mating. It takes two to five times as long |
| Hexaploids, or fewer than about 50 reads, on a panel known to be biased | `--allele-bias fit` | `auto` keeps a bias only where it is clear, which takes depth; `fit` uses every bias |
| Subpopulations at 25 reads or fewer, on a panel without bias | `--prior structure` | Each sample's prior comes from the samples most like it. It cannot be combined with a bias |

[Which settings to use](../commands/call.md#which-settings-to-use) has the measurements behind
each row.

If you come from updog or polyRAD:

| You used | The same model in `call` |
|----------|--------------------------|
| updog, `norm` | `--prior norm --allele-bias fit`. updog always fits a bias |
| polyRAD, `IterateHWE` | The defaults. polyRAD's own settings are `--overdispersion 9 --contamination-rate 0.001 --background-rate 0` |
| polyRAD, `IteratePopStruct` | `--prior structure`. The method differs, and the purpose is the same |

## 5. Write Other Record Types

`--vcf` takes one or more record types, and each goes to its own file:

```bash
madcats call --madc results/demo_MADC.csv --output-prefix results/demo --vcf target all haplotype
```

| `--vcf` | File | A record is |
|---------|------|-------------|
| `target` | `PREFIX.targets.vcf` | The variant a marker was designed for |
| `all` | `PREFIX.all.vcf` | A target, or any other base substitution seen in the tags |
| `haplotype` | `PREFIX.haplotypes.vcf` | A marker, with whole tag sequences as its alleles |

Records with more than two alleles, which is every haplotype record, carry `MPC` in place of `GQ`.
It is the lowest probability among the allele dosages called, and `--min-marginal-posterior`
(default 0.9) is its threshold.

## 6. Check the Run

`call` prints what to look at, and `PREFIX.qc.json` keeps it:

- **How many genotypes were called.** A low share at the default thresholds usually means low
  depth.
- **Markers without a target record.** `call` does not guess a target when a marker's Ref and Alt
  tags differ in more than one place. It says how many markers that left out; pass the panel's
  reference tags with `--reference-tags` to settle them.
- **The overdispersion estimate.** A warning that the estimate is the lowest value tried means the
  table has no heterozygous samples to measure scatter from, as with inbred lines.
- **The allelic bias**, when one is fitted. The `allele_bias` block counts the records called with
  their bias and those called without it, by reason. A warning that most fitted biases ended at a
  limit of the search means the samples have few heterozygotes: call them without a bias.

## 7. Filter Markers Afterwards

No genotype threshold removes a marker that does not behave as one locus. On two potato panels
with SNP-array genotypes for every marker:

- **Filter markers on call rate first.** The share of samples that pass the default thresholds
  separated the markers that disagreed with the array from the rest better than depth, `GQ`, the
  fitted bias or allele frequency. Removing markers with a call rate under 0.9 raised agreement
  with the array from 91.3% to 94.1% on one panel and from 93.3% to 95.0% on the other.
- **Raise `--min-gq` to 20 at most.** Beyond that, most of the genotypes lost were right ones.
- **`INFO/BIAS` is a weak filter.** A strongly biased marker is called right once its bias is
  fitted. Removing markers with a bias outside 0.5 to 2 took out 62 good markers among 162 on one
  panel and 121 among 179 on the other.

These figures are from runs with the bias fitted.

## Know the Limits

- **Allopolyploids.** `call` takes every chromosome copy to be able to carry any allele. A marker
  that amplifies all subgenomes but varies on one is outside that model. A marker that amplifies
  one subgenome is a diploid marker: call it with `--ploidy 2`.
- **Pedigrees.** There is no prior for a cross with known parents, and linkage between markers is
  not used. polyRAD has both.
- **Haplotype records.** Their accuracy has been measured on simulated tables only.

## Next Steps

- [madcats call](../commands/call.md): every option, the model and its measurements
- [Output Formats](../reference/output-formats.md): the VCF fields and the QC JSON
- [Validation & Accuracy](../reference/validation.md): what `call` was checked against
