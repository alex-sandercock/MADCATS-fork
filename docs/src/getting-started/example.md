# Try the Example

MADCATS comes with a small simulated DArTag-style project that runs in seconds. The true
haplotypes and genotypes behind its reads are included, so you can check every row of the output
against what was simulated.

## Get the Data

Download it from the latest release (MADCATS 0.13.0 and later):

```bash
curl -LO https://github.com/Breeding-Insight/MADCATS/releases/latest/download/madcats-demo.tar.gz
tar -xzf madcats-demo.tar.gz
cd madcats-demo
```

The same files are in the repository under `examples/madcats-demo`.

| File | Contents |
|------|----------|
| `primers.csv` | DArT-style primer CSV for the 8 markers |
| `probe_design.tsv` | Probe design with one bracketed SNP per marker |
| `fastq/` | One FASTQ per sample, `S01.fastq.gz` … `S12.fastq.gz` |
| `sample_names.csv` | Renames the samples from `S01` to `Line_01`, and so on |
| `truth_haplotypes.csv` | Every simulated haplotype's 81-bp tag sequence |
| `truth_genotypes.csv` | Each sample's simulated diploid genotype at each marker |
| `make_demo_data.py` | The script that made all of the above; rerunning it reproduces the files exactly |

## What Was Simulated

Twelve diploid samples, eight markers with 81-bp tags, and about 550 reads per sample. Beyond each
marker's designed Ref and Alt alleles:

- **Chr02** and **Chr04** carry an extra microhaplotype: a second SNP on the Ref or Alt background.
- **Chr06** carries a haplotype with an undesigned 1-bp insertion.
- **Chr08** has no Alt allele in any sample.
- **Chr03** and **Chr07** have their primer on the reverse strand of the probe-design target.
- Every read has a 0.2% chance per base of a sequencing error, about 8% of reads are shorter than
  the tag, and each sample has a few off-target reads.

## Run It

Reconstruct the panel from the design files, then count:

```bash
madcats design \
  --primers primers.csv \
  --probe-design probe_design.tsv \
  --output results/reference_tags.csv

madcats generate \
  --fastq-dir fastq \
  --reference-tags results/reference_tags.csv \
  --sample-lookup sample_names.csv \
  --output results/demo_MADC.csv
```

`generate` writes the MADC table, `results/demo_MADC.csv.qc.json` and
`results/demo_MADC.csv.variants.csv`.

## Check the Results

The table has 19 rows: a Ref and an Alt row for each of the 8 markers, plus the 3 microhaplotypes
MADCATS discovered. With MADCATS 0.13.0 the discovered rows are:

| AlleleID | Variants | MarkerType | Reads (all samples) |
|----------|----------|------------|--------------------:|
| `Chr02_001257730\|RefMatch` | `67:A>T` | `snp` | 150 |
| `Chr04_002318006\|AltMatch` | `25:G>T;57:A>T` | `snp` | 122 |
| `Chr06_001644492\|RefMatch` | `29:G>GC;clip=1` | `indel` | 128 |

Things to look for:

- **Every simulated haplotype is a row, and nothing else is.** Compare the `AlleleSequence` column
  with `truth_haplotypes.csv`. The sequencing errors were merged back into their true haplotypes by
  denoising, so none of them became rows.
- **The insertion was recovered.** Its reads failed the substitution-only identity check and were
  placed by gapped alignment. `Variants` reports the insertion left-aligned (`29:G>GC`), and
  `clip=1` notes that it pushed the last Ref base out of the 81-bp window.
- **Chr08 still has an Alt row.** Canonical Ref and Alt rows are always written. Its single read
  comes from a sequencing error that happened to hit the SNP base.
- **Counts follow the genotypes.** In each `Line_` column, a heterozygous sample splits its reads
  between two haplotypes and a homozygous sample puts nearly all of them on one. Compare with
  `truth_genotypes.csv`.
- **Every read is accounted for.** The `read_accounting` section of the QC JSON shows 6,989 reads in
  total, 129 assigned by gapped alignment, 466 partial reads added to counts, 851 merged away by
  per-sample denoising, and 64 in error sequences that failed the evidence rules. Each sample's
  entry in `sample_qc` also counts its off-target reads (`reads_no_primer`).

## Try Changing Something

Each of these reruns `generate` with one setting changed:

| Change | What happens |
|--------|--------------|
| `--reproducible-min-samples 5 --depth-min-reads 0` | The 3 discovered rows disappear: each extra haplotype is carried by only 4 samples, and the depth route is off |
| `--gapped-rescue off` | The Chr06 insertion row disappears; its reads are no longer assigned |
| `--partial-read-mode off` | Same 19 rows, but reads in rows fall from 6,661 to 6,195 |

To compare two runs, use `validate`:

```bash
madcats validate --generated results/other_MADC.csv --expected results/demo_MADC.csv
```

## Next Steps

- [Choosing a Workflow](workflows.md): which commands your own data needs
- [Calling Genotypes from a MADC](../tutorials/genotype-calling.md): genotypes and VCFs from this example's MADC
- [Evidence Model](../reference/evidence-model.md): how reads are assigned, discovered and denoised
- [Output Formats](../reference/output-formats.md): every column of the MADC table
