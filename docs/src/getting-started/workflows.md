# Choosing a Workflow

Every `generate` run consumes one **panel definition**: a reference-tag CSV holding each marker's
fixed-length Ref and Alt tag sequences plus an assignment anchor. Which commands you run before
`generate` depends on two questions: *how will you get the panel?* and *are your reads
primer-started?*

## Step 1: Define the Panel

| You have | Use | Notes |
|----------|-----|-------|
| Vendor primer CSV + probe-design file | [`madcats design`](../commands/design.md) | Uses no reads and gives the same result every time. The other routes were validated against it. |
| A vendor-delivered MADC from the same panel | `scripts/madc_to_reference_tags.py` | Reproduces the vendor's tags exactly, including its IUPAC curation. |
| Probe-design file only | [`madcats infer-design`](../commands/infer-design.md) | Selects each tag window from read evidence. Freeze the result and reuse it. |
| A curated Ref/Alt tag table from any source | Use it directly | `generate` accepts any reference-tag CSV in the documented format, whatever made it. See [Input Formats](input-formats.md). |
| Reads that are not primer-started | [`madcats normalize-reads`](../commands/normalize-reads.md) | Produces the panel *and* pipeline-ready FASTQs together. |

## Step 2: Match the Read Geometry

| Read geometry | Platforms | Route |
|---------------|-----------|-------|
| Single-end, every read starts with the primer | DArTag, GT-seq R1 | Straight to `generate` |
| Raw DArT delivery (barcode + common sequence before the amplicon) | Older DArTag orders | `generate --leading-anchor auto` |
| Bidirectional, variable start positions | AgriSeq (Ion Torrent) | `normalize-reads`, then `generate` |
| Paired-end, mates overlap | GT-seq PE, rhAmpSeq | `merge-pairs`, then `normalize-reads` or `generate` |
| Paired-end, mates cannot overlap | Flex-Seq (dual-anchored capture) | Pool both mates per sample, then `normalize-reads` |

MADCATS does not detect paired-end data automatically: check how your reads are laid out and choose
the route yourself. Adapter trimming is not required on any route.

## Common Pipelines

### DArTag with design files

```bash
madcats design --primers primers.csv --probe-design probes.txt --output tags.csv
madcats generate --fastq-dir fastqs/ --reference-tags tags.csv --output MADC.csv
```

See [DArTag from Design Files](../tutorials/dartag.md).

### DArTag without the primer file

```bash
madcats infer-design --probe-design probes.txt --fastq-dir fastqs/ --output tags.csv
madcats generate --fastq-dir fastqs/ --reference-tags tags.csv --output MADC.csv
```

See [Panels Without a Primer File](../tutorials/primer-free.md).

### AgriSeq / Ion Torrent

```bash
madcats normalize-reads --probe-design markers.tsv --fastq-dir raw/ --out-dir norm/
madcats generate --fastq-dir norm/ --fastq-glob "*.fastq" \
  --reference-tags norm/window_reference_tags.csv --tag-length 41 \
  --min-supporting-freq 0.10 --output MADC.csv
```

See [AgriSeq and Other Non-Primer-Started Reads](../tutorials/agriseq.md).

### Overlapping paired-end reads

```bash
madcats merge-pairs --fastq-dir raw/ --out-dir merged/
madcats normalize-reads --probe-design markers.tsv --fastq-dir merged/ --out-dir norm/ \
  --group-routing covered
madcats generate --fastq-dir norm/ --fastq-glob "*.fastq" \
  --reference-tags norm/window_reference_tags.csv --output MADC.csv
```

See [Paired-End Reads and Tiling](../tutorials/paired-end.md).

## Optional: Clean Paralogs With a Genome

If co-amplified paralogous copies inflate heterozygosity in your panel, clean the finished MADC
(from `generate`, from the vendor, or from HapApp) against a reference genome:

```bash
madcats filter-paralogs --madc MADC.csv --genome genome.fa --genome-index genome.sr.mmi --output MADC_filtered.csv
```

See [Filtering Paralogs with a Reference Genome](../tutorials/paralog-filter.md).

## Optional: Call Genotypes

`call` estimates each sample's genotype from the MADC (from `generate`, `filter-paralogs`, the
vendor, or HapApp) and writes VCFs. Ploidy can differ between samples:

```bash
madcats call --madc MADC.csv --ploidy 4 --output-prefix genotypes
```

See [Calling Genotypes from a MADC](../tutorials/genotype-calling.md).

## Freeze Windows Across Runs

A primer-defined tag is fixed by the chemistry, but `normalize-reads` places each window based on
the reads, so the same panel re-run on new data can place windows differently and make identical
biology look like new haplotypes. After the first run of a panel, pass its window CSV to every later
run with `--reuse-tags`. On two honeybee chips this raised cross-chip haplotype sharing from 42% to
88%.
