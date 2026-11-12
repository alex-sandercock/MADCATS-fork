# Command Overview

MADCATS is a single executable with nine subcommands. Each handles one stage: defining the panel,
adapting reads, counting, cleaning, checking the result, or calling genotypes.

## Command Structure

All commands follow the pattern:

```bash
madcats <command> [options]
```

Get help for any command with `--help`:

```bash
madcats --help              # List all commands
madcats generate --help     # Help for a specific command
madcats --version           # Print the version
```

When running from a source checkout without installing, use `python3 madcats.py <command>`.

## Available Commands

### Panel Commands

Each produces the reference-tag CSV that `generate` consumes.

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`design`](design.md) | Reconstruct Ref/Alt tags from primer + probe-design files | DArTag panels with full design files |
| [`infer-design`](infer-design.md) | Infer tags from the probe design plus reads | Primer file withheld |

### Read Preparation Commands

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`merge-pairs`](merge-pairs.md) | Merge overlapping paired-end reads | GT-seq PE, rhAmpSeq |
| [`normalize-reads`](normalize-reads.md) | Orient and window non-primer-started reads; emits a window panel | AgriSeq, capture, merged amplicons |
| [`window-scan`](window-scan.md) | Recommend a window length from the project's reads | Before the first `normalize-reads` of a new panel |

### Counting and Checking Commands

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`generate`](generate.md) | Build the MADC count table | Every workflow ends here |
| [`validate`](validate.md) | Compare a generated MADC with a reference MADC | Checking against a vendor delivery |

### Cleaning Command (Optional)

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`filter-paralogs`](filter-paralogs.md) | Drop or flag haplotypes a reference genome places at another locus; needs `pip install 'madcats[genome]'` | Panels whose co-amplified paralogs inflate heterozygosity; MADCATS, vendor or HapApp fixed-ID MADCs |

### Genotype Calling Command

| Command | Description | Primary Use |
|---------|-------------|-------------|
| [`call`](call.md) | Call genotypes from a MADC and write target, all-variant or haplotype VCFs | Diploid and polyploid genotypes without leaving MADCATS; MADCATS, vendor or HapApp fixed-ID MADCs |

## Typical Workflows

### DArTag

```bash
madcats design --primers primers.csv --probe-design probes.txt --output tags.csv
madcats generate --fastq-dir fastqs/ --reference-tags tags.csv --jobs 8 --output MADC.csv
madcats validate --generated MADC.csv --expected vendor_MADC.csv --report validation.json
madcats call --madc MADC.csv --output-prefix genotypes --vcf target haplotype
```

### Paralog-heavy panel (optional genome filter)

```bash
madcats generate --fastq-dir fastqs/ --reference-tags tags.csv --min-identity 0.95 --output MADC.csv
madcats filter-paralogs --madc MADC.csv --genome genome.fa --genome-index genome.sr.mmi --output MADC_filtered.csv
```

### AgriSeq

```bash
madcats window-scan --probe-design markers.tsv --fastq-dir raw/ --lengths 31,41,61,81 --out-dir scan/
madcats normalize-reads --probe-design markers.tsv --fastq-dir raw/ --out-dir norm/ --tag-length 41
madcats generate --fastq-dir norm/ --fastq-glob "*.fastq" \
  --reference-tags norm/window_reference_tags.csv --tag-length 41 \
  --min-supporting-freq 0.10 --output MADC.csv
```

See [Choosing a Workflow](../getting-started/workflows.md) for every route.

## Common Options

These options appear on several commands:

| Option | Description |
|--------|-------------|
| `--help`, `-h` | Display help information |
| `--version` | Display the version (top level only) |
| `--fastq-dir` | Folder holding one FASTQ(.gz) per sample |
| `--fastq-glob` | Filename glob within `--fastq-dir` (default `*`) |
| `--probe-design` | Probe/marker design TSV or CSV with bracketed variants |
| `--tag-length` | Tag or window length in bases |
| `--jobs` | Samples processed in parallel (default: up to 4) |
| `--output` / `--out-dir` | Output file or folder |
| `--report` / `--qc-output` | QC/provenance JSON path (each command has a default beside its output) |

## QC and Provenance

Every command that processes reads writes a JSON report recording its parameters, version, inputs
and read accounting, beside its main output. Keep it with the result, since it records exactly how
the table was made.

## Exit Codes

| Code | Meaning |
|------|---------|
| `0` | Success |
| `2` | Error: invalid arguments, missing or malformed input, or a design inconsistency (the message names the cause) |
