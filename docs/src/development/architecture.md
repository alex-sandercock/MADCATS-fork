# Architecture

MADCATS is three Python modules. `madcats.py` holds every command but two and uses only the Python
standard library. `madcats_paralogs.py` holds `filter-paralogs`, the one command with a dependency
outside it: the `mappy` aligner, imported when the command runs (`pip install 'madcats[genome]'`),
so no other command depends on it. `madcats_call.py` holds `call`, the genotype caller, which also
uses only the standard library and is kept apart because it is a separate subject: a statistical
model on a finished table, not read processing. There is nothing to compile, and the same files run as a script
(`python3 madcats.py`), as the installed `madcats` command, or frozen into a standalone binary.

## Repository Layout

| Path | Contents |
|------|----------|
| `madcats.py` | Every command but `filter-paralogs` and `call`: library functions and the command line |
| `madcats_paralogs.py` | `filter-paralogs`: genome alignment, marker placement, labels and the table filter |
| `madcats_call.py` | `call`: the genotype model, what each VCF type counts, and VCF writing |
| `scripts/make_polyrad_fixture.R` | Regenerates `tests/fixtures/polyrad_parity.json`, the polyRAD values the `call` tests compare against (needs R and polyRAD; the tests do not) |
| `tests/` | Unit and end-to-end tests (`unittest`) |
| `scripts/` | Diagnostic and conversion scripts (e.g. `madc_to_reference_tags.py`, `validate_inferred_design.py`) |
| `docs/` | This documentation (mdBook), plus `DEVELOPMENT.md`, the long-form development record |
| `examples/madcats-demo/` | Simulated example project and the script that generates it |
| `assets/` | Full-size logo (2028 × 2028 PNG); the docs use a 600 px copy in `docs/src/images/` |
| `packaging/` | PyInstaller entry point for the standalone binaries |
| `pyproject.toml` | Package metadata for PyPI |
| `.github/workflows/` | Tests on every push; PyPI and binary releases on version tags |

## Reading Map of `madcats.py`

The module is laid out in the order a read travels. Each section opens with a banner comment:

| Section | Contents |
|---------|----------|
| 0. Constants, errors and data types | `Locus`, `Sample`, `SampleResult`, MADC column definitions |
| 1. Sequence and file utilities | IUPAC matching, mismatch counts, CSV/FASTQ helpers |
| 2. Alignment engine | Myers bit-vector prescreen, banded DP, `accepted_gapped_reading` / `row_reading` |
| 3. Panel definition | Design-file reconstruction (`build_panel`), reference-tag tables, primer-free inference |
| 4. Platform intake | Pair merging, window location and freezing, tiling and routing, window-scan |
| 5. Per-sample read assignment | `PrimerIndex`, `select_marker`, `process_sample` and its cache |
| 6. Evidence, denoising and partial-read salvage | `EvidenceRule`, `collect_evidence`, `finalize_counts`, `apply_partial_reads` |
| 7. Classification and annotation | `classify_sequence`, `variant_events`, `marker_type_label` |
| 8. MADC input/output and QC | `read_madc`, `write_madc`, `validate_madc` |
| 9. Commands and command line | One `run_*` per subcommand, `build_parser`, `main` |

`madcats_paralogs.py` follows the same order on a smaller scale: scoring and alignment
(`GenomeAligner`, `WindowIndex`), marker placement and labels (`place_marker`, `score_at_sites`,
`genome_label`), the two-pass table filter (`filter_paralogs`), then its command line. It imports
the MADC helpers it shares with `validate` from `madcats`, so `build_parser` imports it when the
parser is built rather than at the top of `madcats.py`.

`madcats_call.py` is imported the same way. It opens with the genotype model (`call_marker` and the
functions it is built from, adapted from polyRAD), then the step from per-allele posteriors to one
genotype (`correct_dosages_to_ploidy`), input (`read_count_table`, `read_ploidy_table`), what each
VCF type counts (`resolve_target`, `RowReading`, `process_marker`), VCF text (`format_record`,
`vcf_header`) and its command line. It reads rows against the Ref tag with `variant_events`, the
same reading that fills the MADC's `Variants` column.

## The `generate` Pipeline

```
FASTQs ──► 5. assign reads ──► (per-sample cache)
              │  anchor prefix → identity check → gapped alignment
              ▼
           6. partial reads
              ▼
           6. evidence (reproducibility / depth / quality routes)
              ▼
           6. denoise (per-sample skew, then pooled same-allele merge)
              ▼
           7. classify + annotate each row
              ▼
           8. write MADC + QC JSON + variants sidecar
```

Rows hold the exact sequences observed in reads, with no consensus correction, so the output can be
checked against the raw FASTQs.

Decisions that several stages need are made in one place. The same `EvidenceRule` decides
discovery, which sequences denoising protects, and partial-read crediting, and the same
`row_reading` is used for classification and annotation, so those stages always agree.

## Parallelism and Caching

- Per-sample read assignment, pair merging, normalization and inference run across a
  `ProcessPoolExecutor` sized by `--jobs`.
- Since 0.12.1, `generate`'s finalisation also runs across the workers: per-sample denoising per
  sample, and the pooled pass per locus. The result is identical to the sequential passes.
- The per-sample cache (`--cache-dir`) is keyed by the panel and every read-assignment parameter,
  so downstream thresholds never invalidate it. `CACHE_FORMAT_VERSION` is bumped whenever the
  serialized cache or raw-read processing changes.

## Vocabulary

| Term | Meaning |
|------|---------|
| *tag* | The fixed-length sequence a read is measured as |
| *window* | The position of that sequence on a longer read (intake code) |
| *marker* / *locus* | The same panel entry (`Locus` is the record) |
| *canonical* | The design's Ref and Alt tags |
| *category* | Ref, Alt, RefMatch, AltMatch or Other, read at the target base |
| *anchor* | The prefix a read is assigned by: the left primer or a pseudo-primer |
