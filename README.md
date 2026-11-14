# MADCATS <img src="https://raw.githubusercontent.com/Breeding-Insight/MADCATS/main/docs/src/images/madcats-logo.png" align="right" width="250" alt="MADCATS logo" />

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![PyPI](https://img.shields.io/pypi/v/madcats.svg)](https://pypi.org/project/madcats/)
[![Tests](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml/badge.svg)](https://github.com/Breeding-Insight/MADCATS/actions/workflows/tests.yml)
[![Documentation](https://img.shields.io/badge/docs-readthedocs-blue.svg)](https://madcats.readthedocs.io)

**Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing.**

DArTag and similar targeted-sequencing panels are usually delivered as a MADC table of
microhaplotype read counts, made by the vendor's own pipeline. `madcats` builds that table from the
FASTQs and the panel design, with every threshold on the command line and a JSON record of each
run. It can then call genotypes from the table.

- **Counts alleles and finds new ones.** The alleles a panel was designed for are counted, and a
  new microhaplotype is reported when several samples carry it or one sample has many reads of it.
  Indels are handled, including ones the panel was not designed for.
- **One code path for every platform.** Tested on DArTag, AgriSeq (Ion Torrent), GT-seq, rhAmpSeq
  and Flex-Seq capture.
- **Close to the vendor's table.** Read counts correlate at r = 0.9997 with the vendor MADC on a
  sweetpotato DArTag panel, and dosage calls agree 95–99% across five DArTag panels.
- **Genotypes from the same tool.** `madcats call` writes diploid and polyploid genotypes to VCFs
  of the designed targets, of all tag variants, or of whole microhaplotypes. The MADC can also go
  through HapApp to BIGr and other tools, as before.
- **No dependencies.** Python 3 and its standard library. One optional command, the genome-based
  paralog filter (`madcats filter-paralogs`), needs `mappy`.

> **A prior MADC is not an input.** Generation uses only the FASTQs plus the assay design: the
> DArTag primer/probe files, a probe design alone, or a platform-neutral Ref/Alt tag manifest.

## Installation

### From PyPI (Recommended)

Requires Python 3.9 or newer.

```bash
pip install madcats
madcats --help
```

To keep it isolated from other Python packages, use [pipx](https://pipx.pypa.io/):
`pipx install madcats`.

### Pre-built Binaries

Standalone executables that need no Python installation are attached to each
[release](https://github.com/Breeding-Insight/MADCATS/releases):
- `madcats-linux-x86_64.tar.gz`: Linux
- `madcats-macos-x86_64.tar.gz`: macOS (Intel)
- `madcats-macos-aarch64.tar.gz`: macOS (Apple Silicon)

```bash
# Download and extract (example: Linux)
curl -LO https://github.com/Breeding-Insight/MADCATS/releases/latest/download/madcats-linux-x86_64.tar.gz
tar -xzf madcats-linux-x86_64.tar.gz
./madcats --help

# Optional: install to PATH
mkdir -p ~/bin && mv madcats ~/bin/
# Add to PATH if not already (add to ~/.bashrc or ~/.zshrc):
# export PATH="$HOME/bin:$PATH"
```

### From Source

```bash
git clone https://github.com/Breeding-Insight/MADCATS.git
cd MADCATS
pip install .            # or: python3 madcats.py --help (no install needed)
```

## Quick Start

To try MADCATS on a small simulated project first, see [Try the Example](https://madcats.readthedocs.io/en/latest/getting-started/example.html).

```bash
# 1) Reconstruct each marker's Ref/Alt tags from the DArTag design files
madcats design \
  --primers panel_primer_sequence.csv \
  --probe-design panel_probe-design.txt \
  --output reference_tags.csv

# 2) Generate the MADC count table from one FASTQ per sample
madcats generate \
  --fastq-dir fastqs/ \
  --reference-tags reference_tags.csv \
  --cache-dir work/cache \
  --jobs 8 \
  --output results/project_MADC.csv

# 3) Optional: compare against a vendor-delivered MADC
madcats validate \
  --generated results/project_MADC.csv \
  --expected vendor_MADC.csv \
  --report results/validation.json

# 4) Optional: call genotypes and write a VCF
madcats call \
  --madc results/project_MADC.csv \
  --ploidy 4 \
  --output-prefix results/project
```

No primer file? Infer the panel from the probe design and the reads:

```bash
madcats infer-design \
  --probe-design panel_probe-design.txt \
  --fastq-dir fastqs/ \
  --output reference_tags.csv
```

Reads that are not primer-started (for example AgriSeq on Ion Torrent) are normalized first:

```bash
madcats normalize-reads --probe-design markers.tsv --fastq-dir raw_fastqs/ --out-dir normalized/ --jobs 8
madcats generate --fastq-dir normalized/ --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 41 --min-supporting-freq 0.10 --output project_MADC.csv
```

## Documentation

Full documentation is available at [madcats.readthedocs.io](https://madcats.readthedocs.io), including:

- [Getting Started](https://madcats.readthedocs.io/en/latest/getting-started/installation.html): installation and first steps
- [Tutorials](https://madcats.readthedocs.io/en/latest/tutorials/dartag.html): walkthroughs for DArTag, primer-free, AgriSeq and paired-end data, the paralog filter and genotype calling
- [Command Reference](https://madcats.readthedocs.io/en/latest/commands/overview.html): every subcommand and option
- [Evidence Model](https://madcats.readthedocs.io/en/latest/reference/evidence-model.html): how alleles are counted, discovered and denoised

The design rationale and version-by-version measurements behind every default are recorded in the
[development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md).

## Commands

| Command | Description |
|---------|-------------|
| `madcats design` | Reconstruct Ref/Alt tags from vendor primer + probe-design files |
| `madcats infer-design` | Infer a reference-tag panel from the probe design plus FASTQ evidence |
| `madcats merge-pairs` | Merge overlapping paired-end reads into single fragments |
| `madcats normalize-reads` | Orient and window non-primer-started amplicon reads (e.g. AgriSeq) |
| `madcats window-scan` | Sweep candidate tag lengths and recommend one from the project's own reads |
| `madcats generate` | Generate the MADC count table |
| `madcats validate` | Compare a generated MADC against a reference MADC |
| `madcats filter-paralogs` | Drop or flag MADC rows that a reference genome places at another locus (needs `mappy`) |
| `madcats call` | Call genotypes from a MADC and write target, all-variant or haplotype VCFs |

Run `madcats <command> --help` for full options.

## Input Formats

### Reference tags (CSV)
The platform-neutral panel definition that `generate` consumes. Produced by `design`,
`infer-design` or `normalize-reads`, or curated by hand.
```
MarkerName,LeftPrimerSequence,RefSequence,AltSequence
Chr01_000123,ACGTTGCAGGTCAAGT,ACGTTGCAGGTCAAGTTCA...,ACGTTGCAGGTCAAGTTCG...
```
`VariantPositionInTagSequence`, `VariantRef`, `VariantAlt`, `RightPrimerSequence` and
`LocusGroup` are optional.

### Probe design (TSV/CSV)
`MarkerName` and `TargetSequence` with one bracketed variant: a SNP such as `[C/T]` or an
anchored indel such as `[A/ACTTG]`.

### FASTQs
One FASTQ (optionally gzipped) per sample in a folder. Sample columns are named by the FASTQ stem,
or by a `--sample-lookup` table.

## Output

| File | Contents |
|------|----------|
| `project_MADC.csv` | MADC table: `AlleleID`, `CloneID`, `AlleleSequence`, summary columns, one read-count column per sample |
| `project_MADC.csv.qc.json` | Parameters, provenance, and per-sample read accounting |
| `project_MADC.csv.variants.csv` | One line per variant event per row (with `--annotations on`, the default) |

## Validation

Each platform was checked against the best independent comparison available for its dataset.
"Genotype agreement" means updog dosage calls made from the MADCATS table (through HapApp and BIGr)
compared with the reference genotypes.

| Platform | Accuracy | Validated against |
|----------|----------|-------------------|
| DArTag (sweetpotato, red clover, blackberry, alfalfa, pecan) | count r = 0.9997; dosage agreement 95.4–98.8% | vendor MADC, vendor/updog dosage, families, skim-seq VCF |
| AgriSeq, Ion Torrent (honeybee, oyster) | 94.9% vs Torrent VariantCaller; 98.0% replicate concordance | vendor caller + replicates |
| GT-seq (delta smelt) | 99.49% | independent BWA/bcftools truth VCF |
| Flex-Seq capture (cranberry) | 95.13% genotype agreement | published population VCF + 20 trios |
| rhAmpSeq (8 grape families) | 99.92–99.99% family transmission | Mendelian constraints + same-reads head-to-head |

Details: [Validation & Accuracy](https://madcats.readthedocs.io/en/latest/reference/validation.html).

## Development

```bash
python3 -m unittest discover -s tests -v
```

See [Contributing](https://madcats.readthedocs.io/en/latest/development/contributing.html) for the
docs build and release process.

## Citation

If MADCATS is useful in your work, please cite:

> Sandercock, A.M. (2026). MADCATS: Marker-Anchored Discovery and Counting of Alleles from Targeted
> Sequencing. GitHub repository: https://github.com/Breeding-Insight/MADCATS

*(A formal article/DOI is planned; this section will be updated when available.)*

For the MADC format, see Breeding Insight's
[MADC file format overview](https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/).

**AI assistance:** MADCATS was developed in collaboration with Claude Opus 5.5 and Fable 5.1
(Anthropic), which helped write its code and documentation. The authors directed the work, reviewed
and tested the results, and are ultimately responsible for the application.

## License

GPL-3.0-or-later
