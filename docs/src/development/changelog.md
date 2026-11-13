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
- `call --allele-bias fit` estimates updog's allelic bias for every record with two alleles and
  calls with it. With the other settings equal, the fitted bias is within 0.1% of updog 2.1.7's and
  the genotypes are the same, on stored reference values. On two tetraploid alfalfa populations it
  raised agreement with published dosages from 90.9% and 90.0% to 92.2% and 91.8% (updog: 93.0% and
  92.7%). Against SNP-array genotypes of 85 tetraploid potato clones it raised the share of matching
  dosages from 92.3% to 96.6% (updog: 96.7%). The VCF gains `INFO/BIAS`, and `PREFIX.bias.csv`
  lists each fit.
- `call` estimates a background-read rate (`--background-rate auto`, the default) together with
  the overdispersion: the share of reads that show an allele they did not come from. It replaces
  polyRAD's contamination rate, which is now 0 unless set, because contamination cannot account for
  stray reads of an allele no other sample carries. The estimate gives every dosage its own
  probability at each marker, so it does not depend on `--prior` and also works for inbred and
  partly inbred samples, and it uses the worker processes.
- Two-allele records are computed for one allele and mirrored for the other, and terms shared by
  the samples of a ploidy are computed once. Results are unchanged.
- `call --prior norm` adds updog's `norm` genotype distribution: dosage probabilities with a normal
  shape, fitted per record and per ploidy, which does not assume random mating. With the same
  likelihood, its fitted priors and genotypes match updog 2.1.7 on stored reference values. A
  ploidy with fewer than 10 samples uses Hardy-Weinberg proportions.
- `call --selfing-rate` adjusts the Hardy-Weinberg prior for self-fertilisation, for inbred lines
  and partly selfing species. Its priors and iterated results match polyRAD's `selfing.rate`.
- `call` estimates the beta-binomial overdispersion from the read counts by maximum likelihood
  (`--overdispersion auto`, the default; a number sets it). On simulated diploid, tetraploid,
  hexaploid and mixed-ploidy populations it recovered the simulated value. polyRAD's default of 9
  left 59% of target genotypes uncalled on two real tetraploid alfalfa populations, against 8–11%
  with the estimated values.
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
