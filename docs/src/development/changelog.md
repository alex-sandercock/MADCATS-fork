# Changelog

Notable changes by version. Measurements behind each change are in the
[development record](https://github.com/Breeding-Insight/MADCATS/blob/main/docs/DEVELOPMENT.md) on GitHub.

## Unreleased

- New command `call`, in its own module, `madcats_call.py`. It calls genotypes from a MADC (from
  `generate` or `filter-paralogs`, a vendor MADC, or HapApp's fixed-allele-ID output) and writes
  VCFs: `--vcf target` (the designed variant of each marker; the default), `all` (targets plus
  other base substitutions in the tags) and `haplotype` (whole tag sequences as alleles), one file
  per type. Ploidy is set per sample (`--ploidy`, `--ploidy-file`). The model is adapted from
  polyRAD: a beta-binomial likelihood with a uniform or Hardy-Weinberg prior, iterated per record.
  Two-allele records carry `GQ`, `GP` and `PL`; records with more alleles carry `MPC`, the lowest
  marginal posterior among the called dosages. Coordinates are marker-local. It uses only the
  standard library. The tests compare its likelihoods, posteriors, converged allele frequencies
  and genotypes with polyRAD 2.0.1 on stored reference values.
- New optional command `filter-paralogs`, in its own module, `madcats_paralogs.py`. It cleans a
  MADC (from `generate`, a vendor DArT MADC with or without its preamble, or HapApp's
  fixed-allele-ID output) against a reference genome. Each marker is placed by its Ref tag
  and every haplotype is labelled `on_target`, `paralog` or `ambiguous`. Paralog and ambiguous
  rows are dropped (default) or flagged in a `GenomeLabel` column; canonical rows are always kept.
  It writes a labels sidecar, a markers sidecar and a QC JSON. It needs `mappy`
  (`pip install 'madcats[genome]'`), imported only by this command, so MADCATS still has no
  required dependencies. On the lettuce 3K DArTag panel the filtered vendor MADC agreed with a
  bwa-mem-based prototype for 99.4% of reads.
- `filter-paralogs --clone-positions ignore` keeps CloneID coordinates out of marker placement, for
  a genome other than the one the CloneIDs refer to.
- The MADC readers find a header row whose `AlleleID` is not the first cell, as in HapApp's
  fixed-allele-ID output, which starts with `Code_version`.
- `validate` and `filter-paralogs` warn when a MADC has negative read counts. Some vendor MADCs
  store counts in a signed byte: counts of 128–255 appear as count − 256, and counts of 256 or more
  wrap silently. A negative count with a maximum of exactly 127 is reported as this overflow, and
  both QC reports record it. Counts are never corrected, since the wrapped ones cannot be recovered.
  The lettuce 3K vendor MADC has 90,913 such cells across 1,108 markers.

## 0.13.0

- Removed the `madc_pipeline` compatibility module. Import `madcats` instead; the bundled
  scripts already do.
- `--partial-read-mode` takes `on`/`off` like the other switches; the old value `assign` still
  works as `on`, with a warning.
- Clearer command-line help: every command has a description, every option has help text, and
  the program is named `madcats` in usage lines.
- The development record (`docs/DEVELOPMENT.md`) is no longer published on the docs site; it
  stays in the repository.
- A small simulated example project, with its true haplotypes and genotypes, is included in
  `examples/madcats-demo` and attached to each release as `madcats-demo.tar.gz`.
- The test suite checks that the command reference pages match each command's `--help`.

## 0.12.1

- Packaged for PyPI: `pip install madcats` installs a `madcats` command.
- Standalone binaries for Linux and macOS are built for each release.
- New documentation site; the former README is now the development record, `docs/DEVELOPMENT.md`.
- Licensed under GPL-3.0-or-later.
- `generate`'s finalisation (per-sample denoising and the pooled pass) runs across the `--jobs`
  workers. The selection is identical to the sequential passes. On a 44-run rhAmpSeq family the
  pooled pass alone had taken 1,123 s in one process.

## 0.12.0

- Raw DArT deliveries: `generate` and `infer-design` take `--leading-anchor SEQ|auto`. Reads are
  cut where the common sequence ends; runs on an undeclared constant prefix stop and name it.
- Design route: a reconstructed tag must contain the bracketed variant, so primer occurrences
  in duplicated flanks are no longer candidates.
- Denoiser scaling: at distance 1, parents are found by enumerating single-substitution
  neighbours once centroids outnumber them, fixing quadratic time at highly paralogous loci.
  Output is identical.

## 0.11.0

- Partial reads (`--partial-read-mode assign`, now `on`; default): reads shorter than the tag are
  added to an allele that already passed the evidence rules and whose start they match, only when
  this cannot change an allele ratio.
  `--min-partial-length` defaults to tag length − 9. The experimental `extend` mode was removed.
- One reading per row: classification and annotation share `row_reading`, so they cannot
  disagree.
- Freeze health: `normalize-reads --reuse-tags` scores frozen windows against the new run's reads.
- `window-scan --scan-draws` reports recommendation stability across seeded sample subsets.

## 0.10.0

- Undesigned indel recovery: reads that fail substitution-only identity get one gapped
  alignment under the same edit limit (`--gapped-rescue`, `--max-indel-length`), with substitution
  and gap limits that keep paralogs out and a bit-vector prescreen for speed.
- Target base read through the gapped alignment when that reading is cheaper and verified.
- Annotations (`--annotations on`, default): `Variants` and `MarkerType` columns and the
  `OUTPUT.variants.csv` sidecar.
- The pooled denoising pass merges only into a same-allele parent; N-bearing tags are credited or
  dropped, never emitted.
- `window-scan` gains `--scan-seed` and `--generate-args`.
- Cache format 0.5.0 (caches rebuild once).

## 0.9.0

- Paired-end intake: the `merge-pairs` command, and pooled mates for dual-anchored designs.
- Multi-window tiling: the `LocusGroup` column, `--group-routing nearest|covered`, and the
  check that tiles excise distinct windows.
- The `window-scan` command.
- `--primer-mismatches` default changed from 0 to 1.

## 0.8.3

- `normalize-reads` assigns anchors per locus, growing past `--anchor-length` only where prefixes
  collide.

## 0.8.2

- Uncomputed statistics columns emit `.` instead of blank cells.
- `--dart-preamble`, `--preamble-project` and `--sample-lookup` for HapApp-compatible output.

## 0.8.1

- `normalize-reads --reuse-tags` freezes window placement across runs.
- Auto placement penalizes tandem repeats (`--str-penalty`).

## 0.8.0

- The `normalize-reads` command for bidirectional, variable-start platforms such as AgriSeq.

## 0.7.0

- Indel markers and degenerate (IUPAC) primers: one general window rule for SNPs, insertions,
  deletions and homopolymer shifts; tags built as primer overlay plus genomic suffix.

## 0.5.0

- Evidence model redesign: reproducibility (primary), depth and optional quality routes,
  combined by OR. Flags renamed so old commands stop with an error:
  - `--min-quality-support` → `--quality-min-clean-reads` (route now off by default)
  - `--novel-min-reads-per-sample` / `--novel-min-samples` → `--reproducible-min-reads` /
    `--reproducible-min-samples`
  - `--depth-min-reads` added
  - `--allele-catalogue` / `--known-only` removed (now downstream HapApp responsibilities)
- Per-sample processing 7.1x faster through memoized assignment decisions.

## Earlier

- `infer-design` for panels without a primer file.
- 0.3.x–0.4.x: quality-gated discovery and the exploratory audits that motivated 0.5.0.
- The tool was renamed from `madc_pipeline` to MADCATS, with `madc_pipeline.py` kept as a
  compatibility shim until 0.13.0.
