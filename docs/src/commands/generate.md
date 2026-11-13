# madcats generate

Generate a MADC microhaplotype count table from per-sample FASTQs and a panel.

## Synopsis

```bash
madcats generate --fastq-dir <DIR> --reference-tags <FILE> --output <FILE> [OPTIONS]
madcats generate --fastq-dir <DIR> --primers <FILE> --probe-design <FILE> --output <FILE> [OPTIONS]
```

## Description

Every workflow ends with `generate`. It runs in seven stages, and the option tables below are
grouped by stage:

1. Load the panel from a reference-tag CSV (or build it from DArT design files) and find one FASTQ
   per sample.
2. Assign each read to a marker by its left-primer/anchor prefix, keep the fixed-length tag, and
   accept it if it reaches `--min-identity` to the marker's Ref or Alt tag. Reads that fail for
   every marker get one more try with a gapped alignment under the same edit limit, so undesigned
   indels are counted. A read that fits several tiled markers of one amplicon equally is credited
   to each of them; other ties are dropped. Reads shorter than the tag are added to a full-length
   allele whose start they match. Results are cached per sample.
3. Count exact tag sequences per sample.
4. Apply the evidence rules. Canonical Ref and Alt rows are always written; any other sequence is
   reported only if it passes an evidence route.
5. Denoise. Within each sample, rare sequences one edit away from a much more common one are merged
   into it. A pooled pass then merges sequences that failed the evidence rules into a parent that
   carries the same target allele.
6. Classify each row as Ref, Alt, RefMatch, AltMatch or Other, and annotate it.
7. Write the MADC, the QC/provenance JSON and the variants sidecar.

The [Evidence Model](../reference/evidence-model.md) explains stages 2–6 in detail.

> **Process a cohort in one invocation.** The reproducibility route counts supporting samples
> within the run, so splitting a cohort into batches weakens discovery. `--cache-dir` makes
> re-running one large invocation cheap.

## Required Arguments

| Argument | Description |
|----------|-------------|
| `--fastq-dir <DIR>` | Folder containing one FASTQ(.gz) per sample |
| `--output <FILE>` | MADC CSV to write |
| `--reference-tags <FILE>` | Reference-tag CSV, **or** |
| `--primers <FILE>` + `--probe-design <FILE>` | DArT design files |

## Panel and Samples (stage 1)

| Option | Default | Description |
|--------|---------|-------------|
| `--fastq-glob <GLOB>` | `*` | FASTQ filename glob within `--fastq-dir` |
| `--target-metadata <FILE>...` | `targets_*.csv` in the FASTQ folder | Target-manifest CSV(s) for sample metadata |
| `--sample-name-field <FIELD>` | targetid / FASTQ stem | Unique metadata field used for sample columns |
| `--sample-lookup <FILE>` | none | Two-column CSV/TSV mapping FASTQ-derived ids to sample names. Every run sample must be mapped |
| `--tag-length <INT>` | 81 | Tag length |

## Read Assignment (stage 2)

| Option | Default | Description |
|--------|---------|-------------|
| `--min-identity <FLOAT>` | 0.80 | Minimum substitution-only identity to a Ref or Alt tag; a gapped alignment uses the same edit limit |
| `--leading-anchor <SEQ\|auto>` | none | Constant sequence before the amplicon in raw deliveries (DArT: `GATGGATCCAGTGC`). Each read is cut where the anchor ends (it must end within 40 bases beyond its own length). `auto` detects it from the first FASTQ |
| `--primer-mismatches {0,1}` | 1 | Mismatches allowed while assigning a primer/anchor. With `0`, reads from Match haplotypes whose distinguishing SNP falls inside the anchor are lost |
| `--gapped-rescue {on,off}` | `on` | Give reads that fail substitution-only identity one gapped alignment, so undesigned indels are counted |
| `--max-indel-length <INT>` | 10 | Widest gap (alignment band) a gapped alignment may use |
| `--tile-groups {auto,declared,off}` | `auto` | Link tiled markers (markers that read the same tag window of one amplicon), so a read is credited to every tile it fits instead of being dropped as a tie. `auto` links tiles declared with `LocusGroup` and detects undeclared ones; `declared` links only the declared ones; `off` links none. See [Tiled Markers](../getting-started/input-formats.md#tiled-markers) |
| `--partial-read-mode {on,off}` | `on` | Add a read shorter than the tag to the one full-length allele (already passing the evidence rules) whose start it matches, only when it is longer than every prefix shared between that locus's alleles. `off` drops such reads |
| `--min-partial-length <INT>` | tag length − 9 | Minimum length of a credited partial read |
| `--cache-dir <DIR>` | none | Reusable per-sample count cache |
| `--jobs <INT>` | up to 4 | Samples processed in parallel (also parallelizes denoising) |

## Evidence (stage 4)

| Option | Default | Description |
|--------|---------|-------------|
| `--min-global-count <INT>` | 1 | Minimum project-wide reads for a de novo allele |
| `--reproducible-min-reads <INT>` | 10 | Reproducibility route: reads required per supporting sample |
| `--reproducible-min-samples <INT>` | 2 | Reproducibility route: supporting samples required; `0` disables the route |
| `--depth-min-reads <INT>` | 30 | Depth route: reads within one sample that establish an allele on their own; `0` disables |
| `--quality-min-clean-reads <INT>` | 0 | Quality route: project-wide clean reads that establish an allele; `0` (default) disables |
| `--discovery-min-quality {10,20,30,40}` | 20 | Phred score every tag base must reach for a read to count as clean |
| `--min-supporting-freq <FLOAT>` | 0 | Minimum share of a sample's locus reads before that sample can support the reproducibility or depth route |

## Denoising (stage 5)

| Option | Default | Description |
|--------|---------|-------------|
| `--denoise-distance {0,1,2}` | 1 | Per-sample denoising: the most differing bases allowed between a rare sequence and its parent; `0` disables |
| `--global-denoise-distance {0,1,2}` | 1 | Most differing bases allowed in the pooled same-allele merge; `0` disables |
| `--denoise-alpha <FLOAT>` | 2.0 | UNOISE-style abundance-skew exponent |

## Output (stages 6–7)

| Option | Default | Description |
|--------|---------|-------------|
| `--qc-output <FILE>` | `OUTPUT.qc.json` | QC/provenance JSON |
| `--annotations {on,off}` | `on` | Repurpose `PICSnp` as `Variants` and `AvgPIC` as `MarkerType`, and write `OUTPUT.variants.csv`. `off` keeps the exact vendor header |
| `--number-match-ids` | off | Append stable per-file serials to Match IDs |
| `--dart-preamble` | off | Prepend the vendor-style seven-row starred preamble that HapApp expects |
| `--preamble-project <STR>` | output name before its first dot | Label for the first preamble row |

## Examples

### DArTag with a reference-tag panel

```bash
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/project_MADC.csv
```

### Raw DArT delivery, HapApp-ready output

```bash
madcats generate \
  --fastq-dir raw_fastqs/ \
  --reference-tags reference_tags.csv \
  --leading-anchor auto \
  --dart-preamble \
  --sample-lookup sample_names.csv \
  --output results/project_MADC.csv
```

### High-depth AgriSeq after normalization

```bash
madcats generate \
  --fastq-dir normalized/ \
  --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 41 \
  --min-supporting-freq 0.10 \
  --output results/project_MADC.csv
```

### Stricter discovery

```bash
madcats generate ... --reproducible-min-reads 20
```

### Exhaustive discovery

Setting every route's threshold to 0 emits every sequence with at least `--min-global-count` reads:

```bash
madcats generate ... --reproducible-min-reads 0 --depth-min-reads 0 --quality-min-clean-reads 0
```

## Output

| File | Contents |
|------|----------|
| `OUTPUT` | MADC table; see [Output Formats](../reference/output-formats.md) |
| `OUTPUT.qc.json` | Parameters, provenance, per-sample read accounting |
| `OUTPUT.variants.csv` | One variant event per line (`--annotations on`) |

`OUTPUT.qc.json` and `OUTPUT.variants.csv` append to the full output name, so `--output MADC.csv`
writes `MADC.csv.qc.json` and `MADC.csv.variants.csv`.

## Tuning by Depth

Evidence thresholds are absolute read counts, so they do not scale with sequencing depth.

| Depth regime | Recommendation |
|--------------|----------------|
| Low (~60x) and standard DArTag depth | Defaults |
| Very high (above ~500x per sample per marker) | `--min-supporting-freq 0.10` |
| Unreplicated designs relying on the depth route | Consider `--min-supporting-freq 0.01` |
| Low-depth diploid, more sensitivity needed | `--reproducible-min-reads 5 --depth-min-reads 15 --min-supporting-freq 0.10`, confirming extra rows by segregation or replicates |

## See Also

- [Evidence Model](../reference/evidence-model.md)
- [Output Formats](../reference/output-formats.md)
- [madcats validate](validate.md)
