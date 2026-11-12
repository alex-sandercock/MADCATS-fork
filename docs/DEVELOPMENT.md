# MADCATS development reference

**Marker-Anchored Discovery and Counting of Alleles from Targeted Sequencing**

> This is the long-form development record: design rationale, version-by-version
> changes, and the measurements behind every default. It was the project README
> through 0.12.1. For installation and everyday usage, start with the
> [README](../README.md) or the [documentation site](https://madcats.readthedocs.io).
> Command examples below use `python3 madcats.py`; after `pip install madcats`
> the same commands run as `madcats <command>`.

A platform-independent tool for recovering design-specified alleles and
discovering additional alleles from targeted sequencing reads, with MADC
and VCF-ready output. `madcats.py` uses only Python 3's standard library,
exposes all analysis thresholds on the command line, and writes a QC and
provenance JSON beside every result. (The name also reflects the tool's
origins as a FASTQ-to-MADC workflow, where MADC denotes the Missing
Allele Discovery Counts format.)

**A prior MADC is not an input to the workflow.** Generation uses the FASTQs
plus the assay design: either the DArTag primer/probe files or a
platform-neutral Ref/Alt tag manifest.

Every generation run consumes one panel definition: a platform-neutral
reference-tag CSV holding each marker's primer-started, fixed-length Ref and
Alt tag sequences. Three routes produce that definition, differing only in
which inputs they require:

1. **Design route** (`design` command): the vendor's primer CSV plus the
   bracket-notation probe-design file. Deterministic reconstruction from
   design files alone — no sequencing data involved. Use when both files
   are available; it is the reference implementation the other routes are
   validated against. A tag is the primer-started window that carries the
   bracketed variant, so a primer occurrence in a duplicated flank that
   never reaches the variant is not a candidate (the pecan 3K panel has one
   such marker; before this rule it was reported as ambiguous). Any DArT
   marker table with `MarkerName` and a bracket-notation `TargetSequence`
   serves as the probe design, extra markers and columns included.
   Raw DArT deliveries (older orders) keep an inline barcode and the
   library's common sequence `GATGGATCCAGTGC` in front of every amplicon;
   `--leading-anchor GATGGATCCAGTGC` cuts each read where that sequence
   ends (it must end within 40 bases beyond its own length) and counts
   reads without it as `reads_no_leading_anchor`. Preprocessed deliveries
   are primer-started and need no anchor. The anchor is part of the cache
   key and the provenance line; `infer-design` takes the same option, so
   the primer-free route reads raw deliveries too. Both commands check the
   first FASTQ before any read pass: when almost no read starts with a
   panel tag but most carry one after a constant prefix, the run stops and
   names that prefix (the suffix shared by the bases in front of the tag),
   and `--leading-anchor auto` applies the detected sequence instead,
   recording it in the QC as detected.
2. **Vendor-MADC bootstrap** (`scripts/madc_to_reference_tags.py`): any
   delivered MADC from the same panel. The vendor's canonical Ref/Alt rows
   carry the exact tags, their first difference gives the variant position,
   and a tag prefix substitutes for the primer as the assignment anchor.
   Use when a prior vendor report exists but design files are private; the
   result reproduces the vendor's tags exactly, including their IUPAC
   curation. Markers that cannot become loci — above all Ref and Alt rows
   that carry the same sequence, which means the target lies one base past
   the window (or the design is presence/absence) — are announced with a
   `WARNING:` per reason and listed with their reason in
   `OUTPUT.skipped.tsv`; on one DArTag project that was 46 of 3,086
   markers, invisible until the sidecar existed.
3. **Primer-free inference** (`infer-design` command): the probe-design
   file plus the project's own FASTQs, with the primer file never consulted.
   Candidate tag windows are enumerated from the design targets and
   selected by exact read-prefix votes (an identity-envelope fuzzy pass
   covers loci with no exact support). Use when the primer file is withheld
   and no vendor MADC exists; accuracy is measured per project in the
   validation sections below.

A hand-curated CSV from any other source (for example a public panel
database) is equally valid input — the reference-tags CSV is the contract,
not its producer.

Platforms whose reads are not primer-started take one extra stage first:
`normalize-reads` locates each raw read against the probe-design targets by
k-mer, orients it forward, and excises a fixed tag-length window containing
the variant, emitting pipeline-ready per-sample FASTQs together with the
matching window reference-tags CSV and a QC ledger. The window is the
platform-neutral substitute for a primer-defined tag; after normalization,
generation runs unchanged. Because the window is an analytical choice
rather than a physical one, every later run of the same panel should pass
`--reuse-tags` with the first run's window CSV — freezing placement is
what keeps microhaplotype sequences comparable across sequencing runs.
See the AgriSeq pathway section below.

Paired-end platforms add one geometry decision before normalization —
merge overlapping pairs (`merge-pairs`) or pool non-overlapping mates —
and panels may tile several windows onto one locus (`LocusGroup`); both
are 0.9.0 features documented in their own section below. Adapter
trimming is not required on any route: reads are located middle-out, and
adapter content beyond the fragment cannot enter a fully covered window.

## Performance at a glance

Five platforms validated through one engine, each against the strongest
comparator its data permits. "Usable windows" is the fraction of raw
reads that survive to clean, counted evidence at a panel window (see
the read-accounting guardrails below — some loss is correct rejection
of off-target and dimer reads, and external pipelines' percentages use
different denominators).

| platform (dataset) | reads -> usable windows | accuracy | validated against |
|---|---|---|---|
| DArTag (hexaploid sweetpotato; diploid red clover; mixed-ploidy blackberry) | 56.5% of delivered records — vs the vendor's own 50.9% on the same FASTQs (19.6% of records are shorter than the 81-bp tag and unusable by any pipeline) | red clover 3-bin dosage agreement **96.98%** with the vendor table (DP>=20, 682,226 cells, 0.11.0; updog on both tables 99.73% at DP>=10, 0.9.x); sweetpotato count correlation r=0.9997 on 1.20M shared-row cells (0.11.0); blackberry dosage agreement 95.39% (967,299 cells); blackberry design tags 3,015/3,015 match the vendor canonicals (2,881 identical, 134 differ only by the vendor's IUPAC codes); alfalfa panel-paper validation set (two 2021 orders, 54-bp design, raw deliveries cut at the common sequence, 0.12.0): category agreement 99.60% and 99.83%, dosage agreement 98.40% and 97.40%, F1 and BC1 progeny haplotypes within the parental sets at 98.6–98.7% (vendor table 98.9%); pecan 3K (2023 order, 81-bp design, 0.12.0): category agreement 98.89%, dosage agreement **98.78%** on 1.02M cells, F1 progeny haplotypes within the parental set at 99.82% (vendor 99.86%), 93.3% agreement with skim-sequencing genotypes on 11 shared cultivars (vendor 93.3%) | vendor MADC + vendor/updog dosage calls + families + skim-sequencing VCF |
| AgriSeq, Ion Torrent (2 honeybee chips; 372-sample oyster) | 88.4% | **94.9%** vs Torrent VariantCaller (152,698 cells; core 97.8%; 94.5% before 0.10.0's indel recovery); honeybee replicate concordance **98.0%**, cross-chip 98.2% | vendor caller + replicates |
| GT-seq (96-sample delta smelt) | 84.5% merged | **99.49%** R1-only (98.80% before 0.10.0 — the gain is one real insertion allele the fixed-offset reading misclassified); 98.67% with tiles on +54% cells; end-of-chain (HapApp->BIGr->updog) **98.42%** targets + **98.79%** on 38K off-target SNPs; `infer-design` 368/368 primer-free | independent BWA/bcftools truth VCF |
| Flex-Seq capture (192-sample cranberry) | 81.9% (110-bp windows) | end-of-chain **95.13%** targets (295,068 cells); trio incompatibility **3.6%** vs 6.2% for the comparison VCF's own source; locus recovery 80.7%/99.84% vs vendor 83.1%/99.8% | published population VCF + 20 pedigree trios |
| rhAmpSeq (608 runs, 8 biparental grape families) | 86.6% (cohort-discovery panel; 79.3% on the original 16-sample-discovery panel) | family transmission **99.92-99.99%** (~14.9K locus-family tests); head-to-head vs the platform's own `analyze_amplicon` pipeline on identical reads: **93.0%** exact haplotype-genotype agreement, +19%/+79% more callable cells at 10x/30x, and MADCATS is the more Mendelian-consistent of the pair (99.947% vs 99.850% within the inferred parental four). Cohort-scale discovery (per-family union across all 8 families) designs 1,993/2,000 panel amplicons and delivers +43% genotype calls over the subset-discovery panel at held three-nines transmission, with every added locus 100.000% family-consistent (84,757 calls, zero violations) | Mendelian family constraints + same-reads tool head-to-head (no public truth matrix) |
| synthetic truth suite | — | dosage **100%** on clean/linked/paralog/repeat classes; 50.8% on engineered null alleles (the correct negative); paired-end merge 100% with 98.6% molecular recovery | exact simulated truth |

Version 0.10.0 adds recovery of *undesigned* indel haplotypes through a
general gapped reading (see its section below). On a private in-house
109-bp DArTag project with 447 assembly-predicted indel microhaplotypes
that the vendor MADC lacked, 78% became MADC rows at 97% read recall
(0.4% under 0.9.x), and the SNP calls on every truth-scored set above
held or improved.

Accuracy tracks comparator quality more than platform: where truth is
clean (GT-seq, DArTag) MADCATS scores 98-99%+; against noisier
comparators (TVC, raw mpileup) the *mutual* agreement is 94-95% while
internal instruments (replicates, trios) show MADCATS' own calls at the
cleaner end of each pair — on cranberry, measurably more
Mendelian-consistent than the raw calls behind the comparison VCF.

## Important scope

DArT describes its analytical pipeline as proprietary. Exact vendor
equivalence therefore cannot be guaranteed because its adapter/quality
filtering, partial-read handling, consensus construction, and error-correction
rules are not documented. This implementation is an open and deterministic
replacement, not a claim to be DArT's undisclosed implementation.

The generated file preserves the MADC fields needed for microhaplotype work:
`AlleleID`, `CloneID`, `AlleleSequence`, `ClusterConsensusSequence`, and one
integer read-count column per sample, following the vendor conventions
observed in the supplied report (canonical rows carry their own sequence as
the consensus; Match/Other rows carry the Ref consensus). `AvgCountRef` and
`AvgCountSnp` are nonzero means, and `CallRate` is filled with the fraction of
samples having canonical Ref or Alt reads at the locus. The remaining vendor
dosage/genotype summary fields are downstream calls rather than FASTQ
evidence and are emitted as `.` (missing value; see the output-format
section).

Matching is substitution-only first: reads are compared to same-length
tags at each position. **Designed indel markers are fully supported** —
vendors report indel alleles as fixed-length windows whose suffix shifts at
the indel, and both tags are ordinary same-length sequences to the counting
layer. Since 0.10.0, **undesigned indels are recovered too**: a read that
frameshifts against every canonical tag is offered one gapped reading under
the same edit budget (gaps of up to `--max-indel-length` bases, 10 by
default, counted as edits), and is kept as an ordinary fixed-length
shifted-suffix row when that reading is composed like an indel haplotype
(few substitutions, few gap events). The target base is read through the
alignment. Out of scope remain indels longer than the band and gapped
readings that need many substitutions — those are co-amplified paralogs,
not alleles, and stay off-target.

## Quick start for the supplied sweetpotato project

The defaults are the recommended configuration; no evidence flags are needed:

```bash
python3 madcats.py generate \
  --fastq-dir DSp25-10607_leafRoot_1plate \
  --primers 'SweetPotato_DArTag_BI_Cornell_University_(1.0)_primer_sequence.csv' \
  --probe-design 20221216-BI-Sweetpotato_3K_SNPs_DArTag-probe-design.txt \
  --cache-dir work/count-cache-v05 \
  --jobs 8 \
  --output results/DSp25-10607_v05_default_MADC.csv
```

A stricter, higher-precision view for high-depth projects raises the
reproducibility depth: add `--reproducible-min-reads 20`.

To inspect or reuse the reconstructed 3,120 Ref/Alt tags, run the separate
design command:

```bash
python3 madcats.py design \
  --primers 'SweetPotato_DArTag_BI_Cornell_University_(1.0)_primer_sequence.csv' \
  --probe-design 20221216-BI-Sweetpotato_3K_SNPs_DArTag-probe-design.txt \
  --output sweetpotato_reference_tags.csv
```

Other targeted platforms can pass this same reference-tag format with
`--reference-tags`. The provided DArT report is used only after generation,
as an external validation target:

```bash
python3 madcats.py validate \
  --generated results/DSp25-10607_v05_default_MADC.csv \
  --expected Report-DSp25-10607/DSp25-10607_MADC.csv \
  --report results/validation.v05_default.json
```

## The evidence model

For each marker, the pipeline locates the left primer in the supplied target
sequence in either orientation, extracts the primer-started tag, verifies the
target-SNP position, and constructs the Alt tag by substitution. Each FASTQ
read is assigned by its left-primer prefix, retained as a raw full tag of the
configured length (81 bp by default), and accepted when it has at least 80%
substitution-only identity to Ref or Alt — or, when every candidate rejects
it, under one gapped reading that spends the same 20% edit budget
(`--gapped-rescue on`, the default; see the 0.10.0 section).

Canonical Ref and Alt rows are always emitted, including at zero counts. A
**non-canonical sequence** becomes an output row when it reaches
`--min-global-count` project-wide reads and passes at least one enabled
evidence route:

- **Reproducibility route** (primary; default on): at least
  `--reproducible-min-reads` reads (default 10) in each of at least
  `--reproducible-min-samples` samples (default 2). Identical 81-mers
  recurring independently across samples are strong evidence even when a base
  position is systematically low-quality.
- **Depth route** (default on): at least `--depth-min-reads` reads (default
  30) within a single sample. This keeps private alleles discoverable in
  unreplicated designs such as diversity panels of singleton accessions; on
  the replicated project below it adds a single row, because deep sequences
  there always recur across the paired extracts.
- **Quality route** (default off): at least `--quality-min-clean-reads` reads
  whose every tag base is at or above `--discovery-min-quality` (default
  Q20). On the supplied project this route's marginal precision was 3%
  (2,159 extra rows for 73 additional vendor-matching rows), and whole-tag
  cleanliness vetoes real alleles with one systematically low-quality cycle,
  so it is no longer a default gate.

An optional HaploCharmer-style dosage-aware floor, `--min-supporting-freq`
(default 0 = off), additionally requires a sample's reads for the sequence to
reach that share of the sample's locus total before the sample can support
the reproducibility or depth route. Measured on the supplied project, a 0.01
floor raises precision 2.2 points and costs 0.35 recall points (17,319 rows,
86.2%/77.2%), but it removes 33 of the vendor-confirmed rows the evidence
redesign recovered: amplicon pools inflated by co-amplified paralogs leave
real minor haplotypes at sub-percent shares, so the floor stays off by
default. It is most useful as an artifact guard on the depth route in
unreplicated designs, where around 0.01 is a reasonable starting point.

Evidence thresholds are absolute read counts and therefore do not scale
with sequencing depth. At the honeybee AgriSeq run's ~750x per sample per
marker, 10 reads sits below the Ion error noise floor, and the defaults
admitted 1,752 novel rows of which roughly 92% were error haplotypes
(median carrier share near zero). Sweep-backed guidance by depth regime:
at low depth (~60x) and at standard DArTag depth the defaults stand; at
very high depth (above roughly 500x per sample per marker) enable
`--min-supporting-freq 0.10` — on the honeybee run it collapsed the novel
rows to 147 while raising caller concordance from 89.8% to 90.6%, and
locus-share floors dominate rescaled absolute thresholds because
per-marker depth varies tenfold within a single run.

Setting every route's threshold to 0 emits every sequence at or above
`--min-global-count` (exhaustive discovery). The same rule protects
qualifying sequences from abundance-skew denoising and gates partial-read
extension, so discovery, protection, and extension stay consistent.
Abundance-skew denoising (UNOISE-style, `2^(alpha*d+1)`, closest parent
first) then absorbs low-abundance neighbors within each sample; a pooled
pass (`--global-denoise-distance`, default 1) afterwards folds an
evidence-failing sequence into the protected parent the project-wide skew
names — but only when both carry the same target allele, so no move can
change a sample's Ref/Alt dosage (0.9.x moved reads across alleles, which
could turn a private 25/25 heterozygote into Ref=50). Evidence-passing and
canonical sequences are never absorbed. Reads whose sequences pass no route
remain visible in the QC read accounting rather than silently vanishing. At the default distance of 1 the eligible parents of a sequence are its single-substitution neighbours, which the denoiser enumerates and looks up rather than finding by scanning every centroid; the selection is identical, and a paralogous locus with thousands of unique sequences per sample (the pecan 3K panel has one) no longer makes the pass quadratic.

Classification against the design is IUPAC-aware on the design side: with an
ambiguity-coded design allele, a compatible concrete observed base classifies
as RefMatch/AltMatch. An ambiguous observed base such as N is missing data
and classifies as Other, matching the vendor's convention. Ref/Alt category
assignment and per-row counts always use exact sequence equality, so an
ambiguity-coded canonical tag reports zero exact counts; the pipeline does
not guess which concrete variants a curated IUPAC consensus was meant to
absorb (DArT uses such codes to encode real secondary SNPs, so collapsing
them would merge genuine microhaplotypes — measured on the blackberry
panel). A read whose tag contains an N is never counted as a haplotype of
its own: it is credited to the one N-free sequence at its locus it is
compatible with, or dropped (`reads_n_mapped` / `reads_n_dropped`).

Reads shorter than the tag length are salvaged by default
(`--partial-read-mode assign`): a trimmed prefix is credited to the one
learned full tag it matches, and only when it is longer than every prefix
shared between two learned alleles of its locus (see the 0.11.0 section
"Truncated reads: partial-read salvage with a locus guard"). `off` excludes such
reads. A partial read never establishes an allele: a prefix of nothing
learned is left unmapped. Note that a partial read's quality is judged over
its observed bases only.

## Recovering a panel without the primer file (`infer-design`)

Some providers treat the primer file as proprietary. Its only irreplaceable
content is where each tag starts inside the probe-design target: the tag must
be a tag-length window of the bracket-expanded target containing the SNP, so
inference selects among a closed candidate set (~4x tag-length windows per
locus) rather than reconstructing sequences freely. Because reads are
primer-started, exact read prefixes identify the true window almost
digitally, and Ref/Alt identity comes from bracket order (the first allele is
the reference-genome base — verified at all 3,120 loci of this panel with
zero swaps).

```bash
python3 madcats.py infer-design \
  --probe-design 20221216-BI-Sweetpotato_3K_SNPs_DArTag-probe-design.txt \
  --fastq-dir DSp25-10607_leafRoot_1plate \
  --jobs 8 \
  --output inferred_reference_tags.csv
```

Two evidence tiers, strictly ordered: exact read-prefix votes decide first
and are authoritative; an identity-envelope fuzzy pass (seeded matching
within the same 80% identity rule as generation, with already-resolved tags
acting as competition sinks) runs only for loci with zero exact support and
can never overrule exact evidence. Each locus needs a dominant window
(default 10x the runner-up); the written anchor (`LeftPrimerSequence`) is the
shortest panel-unique tag prefix, since the true primer boundary is
unobservable in seamless reads and unnecessary for assignment. `--tag-length`
sets the panel's reported tag length (54, 81, 109, ...); reads must reach it
and the design flanks must contain a full window, both validated per locus.

Measured blinded on this project (primer file untouched, ~11 minutes):
3,068/3,120 loci inferred — 2,742 by exact evidence (100.0% correct against
the primer-derived truth) and 326 by the fuzzy tier (319 correct, 97.9%);
zero Ref/Alt swaps or variant-position errors; 52 loci excluded with stated
reasons (no reads within identity, support below threshold, or ambiguous
windows). Running `generate` on the inferred panel reproduces 99.1% of the
full-input MADC's rows (novel Match/Other rows: 99.65%) with 100%
classification agreement; count totals run 0.6% high because the short
anchors reject fewer error-carrying reads than the truer, longer primers.
The `*.qc.json` report records each locus's route, support, margin, anchor
length, and flags. The fuzzy tier is the review tier: its errors win by
large margins (paralog haplotype pools), so no internal statistic can catch
them — cross-check `fuzzy_anchor` loci against any delivered MADC's
canonical rows or a public panel asset, and prefer those sources outright
when available (a single delivered MADC reconstructs the panel exactly).
Inference is a one-time operation per panel; freeze its CSV and reuse it.
`scripts/validate_inferred_design.py` compares an inferred panel against a
ground-truth reference-tag CSV, with per-route accuracy breakdowns.

## Migration from 0.4.x

Version 0.5.0 changed the discovery semantics, so the evidence flags were
renamed to make old commands fail loudly rather than silently change meaning:

- `--min-quality-support` became `--quality-min-clean-reads` (and the route
  is now optional and off by default);
- `--novel-min-reads-per-sample`/`--novel-min-samples` became
  `--reproducible-min-reads`/`--reproducible-min-samples`, and the
  reproducibility test moved from a post-hoc AND filter into the core
  evidence rule (routes combine by OR);
- `--depth-min-reads` is new;
- `--allele-catalogue`/`--known-only` were removed: fixed allele-ID
  assignment and the cross-project allele database are downstream HapApp
  responsibilities, and reprocessing against a frozen catalogue is no longer
  this pipeline's job;
- the per-sample cache format version rose to 0.4.0 (fuzzy primer matches now
  resolve across all matched depths by tag identity instead of preferring the
  deepest primer), so caches are rebuilt on first use. Rebuilding is cheap:
  assignment decisions are memoized per unique sequence, which made
  per-sample processing 7.1x faster; the full 80-sample, 154.6M-read project
  processes from scratch in about 4.5 minutes with 8 workers and regenerates
  from a warm cache in about a minute.

## Validation on the supplied project

Against the vendor MADC (19,319 rows), with every retained shared row having
the identical Ref/Alt/RefMatch/AltMatch/Other classification:

| configuration | rows | shared | precision | recall |
|---|---|---|---|---|
| 0.4.0 benchmark (`5/10x4`, quality AND sample filter) | 16,424 | 14,611 | 89.0% | 75.6% |
| 0.5.0 default (reproducibility 10x2 OR depth 30x1) | 17,844 | 14,988 | 84.0% | 77.6% |
| 0.5.0 `--reproducible-min-reads 20` | 16,517 | 14,839 | 89.8% | 76.8% |
| 0.5.0 `--min-supporting-freq 0.01` | 17,319 | 14,921 | 86.2% | 77.2% |

The strict view now dominates the old benchmark on both axes. Across shared
row/sample cells the default's Pearson correlation is 0.99928, 68.2% of cells
are exactly equal, and the shared count total is 2.0% above DArT's.

### Why the evidence model changed

The 0.4.x discovery gate required project-wide reads whose **every** tag base
reaches Q20. Auditing the raw caches showed this vetoed 286 vendor rows that
are present as exact 81-bp FASTQ prefixes with strong support — up to 6,626
reads spread across as many as 70 samples — because one base position is
systematically low-quality (for many, below Q10, so no clean-read threshold
can rescue them). Symmetrically, the old `10x4` post-filter removed 598
vendor rows that had clean-read support but narrower sample spread. The
reproducibility-primary rule recovers both classes; every vetoed-but-
reproducible vendor row is present in the 0.5.0 default output.

The vendor's own file shows DArT applies no within-project support threshold
of this kind: 27.3% of its positive-count rows have fewer than 10 reads in
fewer than 2 samples of its own counts, 4,081 rows have under 10 total reads,
and 1,417 rows are all-zero. Its row list reflects proprietary clustering
plus a cross-project haplotype database — reproducibility evidence at the
database level, which a single-project de novo run cannot see. The
within-project evidence rule is therefore this pipeline's explicit,
configurable substitute, not an imitation of a vendor threshold.

### Independent validation: red clover (diploid, 109-bp tags, low depth)

The DCl25-10856 red clover project (376 samples, 3,000 loci, 109-bp tags,
mean read length ~88 bp so only ~65% of reads reach the tag length, ~60x mean
locus depth) was processed with no primer file: the probe-design targets plus
`infer-design` recovered 2,968 loci — 2,902 exactly matching the vendor's
canonical tags, 56 differing only by the vendor's IUPAC curation, and the
remaining handful confined to loci where the vendor's own tags have zero
reads. `scripts/madc_to_reference_tags.py` converted the delivered MADC into
the equivalent panel (3,000/3,000) as the cross-check route.

Generation at the defaults gave 85.4% precision and 62.9% raw recall with
100% classification agreement and a 0.9975 shared-count total ratio. The raw
recall is dominated by a structural ceiling: 36.7% of the vendor's
positive-count rows fail every within-project evidence route on the vendor's
own counts (their cross-project catalogue at low depth); against
evidence-bearing vendor rows, recall is 94.1%. Softening to
`--reproducible-min-reads 5 --depth-min-reads 15` captured nearly the whole
recoverable band (67.2% recall against a 67.0% ceiling) at a 10-point
precision cost — but carrier-frequency profiling shows why the bare soft
view is not recommended: its vendor-shared additions look like real alleles
(median 44% of their carrier's locus reads, matching the 48% of established
rows), while its generated-only additions sit at a median 7.5% — not a
diploid genotype. The defaults remain the recommendation for diploid and
polyploid projects alike. When a low-depth diploid project needs more
sensitivity, add the diploid-justified frequency floor:
`--reproducible-min-reads 5 --depth-min-reads 15 --min-supporting-freq 0.10`
recovers about 40% of the available recall band while excluding the
artifact-profile rows, and its extra rows should still be confirmed by
segregation or replicates downstream. On this project the pre-0.10 partial
salvage (minimum length 30, no locus guard) added the ~24M primer-bearing
short reads to counts without changing the row set, and lowered count
concordance with the vendor (Pearson 0.976 to 0.941), showing the vendor
does not count short fragments; the 0.10.0 default (minimum length
`--tag-length` minus 9) excludes those fragments, and `--partial-read-mode
off` restores the exact pre-salvage counts when vendor comparability matters.

Segregation-based validation closed the loop on this project. The run
contains a 152-progeny bi-parental mapping family with both parents, and
both MADCs were pushed through the same updog (diploid) + BIGr workflow.
Updog imputes dosages at low or zero depth (about a third of calls sit at
DP=0), so all comparisons are depth-filtered. Dosage concordance between
the two pipelines' calls is 99.73% at DP>=10 in both (99.77% at DP>=20).
Mendelian-violation rates in the family are statistically identical —
3.38% for this pipeline versus 3.35% for the vendor MADC (parents DP>=20,
progeny DP>=10, MPP>=0.9) — so the shared violation floor reflects the
assay and family, not either counting implementation. Family transmission
of microhaplotype rows (present at >=5 reads in >=2 progeny must appear in
a parent) adjudicates the discovery tiers directly: vendor-shared novel
rows are parent-confirmed at 97.4% (the assay background), the default
configuration's generated-only rows at 72.0% — majority genuinely
transmitted alleles absent from the vendor's catalogue — and the soft
band's generated-only rows at only 28.4%, independently confirming the
carrier-frequency verdict that the bare soft view admits mostly artifacts.
The symmetric test on the vendor's own counts completes the picture: only
313 of their 10,175 exclusive novel rows segregate in the family at all
(the rest are read-starved catalogue content this project can neither
validate nor falsify), and those confirm at 67.7% against a 97.8%
consensus background on the same count matrix. Each pipeline's exclusive
rows therefore carry the same real-but-imperfect profile (72.0% versus
67.7%, indistinguishable at n=313), while this pipeline contributes about
three times as many transmission-validated alleles (~633) as it forgoes
(~212) — and the forgone ones sit in known, tunable territory just below
the evidence thresholds.

The canonical Ref/Alt rows — the input to downstream dosage calling —
compare to the vendor's as follows. Of our 5,936 canonical rows, 5,827
(98.2%) are sequence-identical to the vendor's; the 109 differences are 89
vendor-IUPAC-curated sequences plus the Ref/Alt pairs of the 10 loci whose
vendor tags carry zero reads in this project. Across the 2.19 million
matched count cells, Pearson correlation is 0.975 with a mean absolute
error of 1.29 reads and totals within 0.8%. On the dosage-relevant
Ref/(Ref+Alt) ratio, evaluated over the 499,414 locus-by-sample cells where
both files have at least 20 reads of canonical depth, the median ratio
difference is 0.0000 (95th percentile 0.015) and 99.68% of cells imply the
identical naive diploid genotype class. The zero-count canonical structure
also mirrors the vendor's (644 all-zero rows versus their 604) — the
germplasm-divergent loci both pipelines observe.

### Independent validation: blackberry (mixed ploidy, 109-bp tags, indels, degenerate primers)

The DBla25-10927 blackberry project (376 samples of mixed *Rubus* species
and ploidies, 3,015 loci, 109-bp tags, ~59-70% of reads reaching tag
length) is the panel that motivated version 0.7.0's indel support and is
the first project validated through **all three panel-definition routes**.
Its design contains 10 indel markers — encoded with anchored VCF-style
brackets in the probe-design file but unanchored dash notation in the
primer CSV — and, as the data revealed, primers synthesized with
degenerate IUPAC bases over flanking SNPs, whose oligo sequence the
vendor's tags carry verbatim. DArT's delivered sample names were
misassigned for this project (a correction crosswalk exists); this
pipeline's targetid-named columns are immune to that class of error.

Version 0.7.0 replaced the SNP-only window construction with one general
rule — the divergence point of a tag pair is the shared prefix of the two
full allele expansions — which handles SNPs, insertions, deletions,
anchored brackets, and homopolymer shifts identically, and reconstructs
tags as primer-sequence overlay plus genomic suffix so degenerate oligo
bases persist. SNP-panel outputs are byte-identical through the change.

Route-by-route results against the vendor's canonical tags:

- **Design route**: 2,881/3,015 exact; all 134 differences are vendor-IUPAC
  loci; zero unexplained. All 10 indel loci reconstruct exactly, including
  divergence positions shifted past homopolymer runs.
- **Bootstrap route**: exact by construction (3,015/3,015 converted, none
  skipped, indels included). Generation from this panel: 92.9% precision,
  73.7% raw recall against a 75.7% within-project evidence ceiling (23.7%
  of the vendor's positive rows carry no such evidence on their own
  counts), **94.7% fair recall** — the third project inside the 94-95%
  band — with 100% canonical row recovery, Pearson 0.981, and 99.2%
  classification agreement (557 of the 639 disagreements sit at
  IUPAC-canonical loci and 68 at indel loci, where the positional variant
  heuristic and vendor curation legitimately differ). Family transmission
  in the reciprocal Columbia Star x Kotata families (both parents
  sequenced, 92 progeny) reproduces the red clover profile:
  vendor-shared novel rows 98.5% parent-confirmed (n=12,064),
  generated-only rows 74.0% (n=2,027); the test asserts only that progeny
  alleles exist in a parent, so unknown hexaploid inheritance and
  segmental allopolyploidy do not affect it.
- **Primer-free inference**: 2,943/3,015 loci inferred (72 excluded with
  stated reasons); 2,796 exact plus 126 concrete realizations of
  vendor-IUPAC tags — identical under IUPAC compatibility, the
  information-theoretic limit when the oligo's degeneracy is unknowable
  without the primer file — and 18 wrong windows (0.6%). Of the 10 indel
  loci, 9 resolve exactly by read votes (one via the fuzzy tier at a locus
  whose canonicals carry no reads) and 1 is excluded as ambiguous rather
  than guessed.

The 18 wrong windows decompose into three classes, none an inference
defect in read-bearing data: 8 are read-starved loci anchored on scraps;
2 are short-amplicon loci (72-bp and 105-bp amplicons, real reads that
never reach tag length); and ~8 are loci where the vendor counted
thousands of reads for tags occurring zero times in the delivered FASTQs
in any form. The leading explanation for the last class is the vendor's
pre-delivery overrepresented-read filtering (quantified for sweetpotato
above): dominant amplicons are precisely what such a filter removes, the
vendor counts predate it, and the affected loci are then unknowable to any
FASTQ-based pipeline.

A germplasm note from the same run: the project's two "Marion" samples
are genetically different individuals (canonical-count correlation 0.47
where technical replicates approach 1); one is near-clonal with the
"Newberry" sample (r=0.958) and the other matches nothing in the panel —
a plate-provenance question for curation, caught by targetid-level
fingerprinting.

### AgriSeq (Ion Torrent) pathway: honeybee validation

AgriSeq reads violate the primer-started model three ways: amplicons are
sequenced from both ends (~50/50 orientation split), primer digestion and
primers placed outside the submitted targets disperse read starts, and Ion
Torrent's dominant error is the homopolymer indel (a median 24% of reads
per marker carry one across ~100 bp on the honeybee run, up to 77% at
run-rich loci — the rate climbs monotonically with the longest run in the
window). The pathway is two commands:

```bash
# first run of a panel: learn window placement from observed coverage
madcats.py normalize-reads --probe-design markers.tsv \
  --fastq-dir raw_fastqs/ --out-dir normalized/ --jobs 8

# every later run of the same panel: freeze the first run's windows
madcats.py normalize-reads --probe-design markers.tsv \
  --fastq-dir raw_fastqs_run2/ --out-dir normalized_run2/ \
  --reuse-tags normalized/window_reference_tags.csv --jobs 8

madcats.py generate --fastq-dir normalized/ --fastq-glob "*.fastq" \
  --reference-tags normalized/window_reference_tags.csv \
  --tag-length 41 --primer-mismatches 1 --min-supporting-freq 0.10 \
  --output project_MADC.csv
```

The 41-bp default window is measured, not assumed: sweeping 41/61/81 bp on
the 81-marker honeybee validation (24.2M reads, 81 samples of mixed drone
and worker material) showed longer windows lose coverage and indel-free
reads faster than they add context — end-to-end read yield 64%/59%/44% and
genotype call rate 96.6%/93.9%/85.6% — with the 81-bp concordance uptick a
survivor artifact of scoring only the easiest cells. Data-driven window
placement (per marker, from observed read coverage, keeping the variant
out of the anchor) raised adapter yield from 66% to 85.5% of raw reads.
Since 0.8.1, auto placement scores coverage minus a soft tandem-repeat
penalty (`--str-penalty`, default 0.002 per repeat bp beyond 4, unit
sizes 1-3): the original homopolymer tiebreak provably never fired,
because real coverage fractions are never exactly tied, and it was blind
to dinucleotide repeats entirely.

Genotypes called from target-base dosage agree with the vendor's Torrent
caller at 90.6%, and both pipelines agree with an independent pre-existing
reference genotype set at 94.4% versus 94.5% — a dead heat, with disputed
cells splitting 48:48 and the residual disagreement concentrated at a
handful of three-way-inconsistent markers — most of which the second-chip
study below traced to window placement rather than to the loci
themselves. Row accuracy against the Amel_HAv3.1 genome: all 81
real markers' Ref tags occur verbatim at their named coordinates and every
Alt tag is structurally exact; the only four failures are the vendor's
four design-failed markers, whose submitted target sequences do not match
the genome at their named positions — a pre-submission target-versus-
genome check would have caught them before synthesis. Caveats: a single
project, two imperfect arbiters, and mixed-ploidy material, which the
row-level checks tolerate by construction.

### Second chip, replicates, and window freezing (0.8.1)

A second honeybee chip on the same panel (96 wells: 65 unique samples, 30
within-chip technical replicate pairs, one no-template control; 40.8M
reads at 2.5x the first chip's depth) ran the recipe unchanged: 84.8% of
raw reads normalized, the control clean at 119 reads against a
160,000-read median sample, vendor-caller concordance 92.4%, and — the
first truth-free accuracy instrument — 98.2% genotype concordance across
the 30 replicate pairs. The depth-scaling law reproduced exactly:
default evidence thresholds admitted 4,410 novel rows at this depth
(1,752 on chip 1), while `--min-supporting-freq 0.10` pinned the count
near 145 on both chips and improved every concordance metric, settling
the frequency floor as the AgriSeq default.

The chips also exposed a stability failure with consequences for
microhaplotype databases: re-learning window placement per run moved 42
of 85 windows, which made 58% of the second chip's haplotype sequences
spuriously "novel" against the first chip — technical artifacts
indistinguishable from biology in any downstream database. Freezing the
panel with `--reuse-tags` (the emitted window CSV is byte-identical to
the source panel, every row validated against the design, partial
mismatches refused) raised cross-chip haplotype sharing from 42% to 88%,
with the remainder consistent with genuinely new germplasm on the second
chip.

Window placement also turned out to be the dominant cause of apparently
bad markers, through allele-biased read dropout: where a variant allele
extends an adjacent tandem repeat, its reads carry longer runs and die
preferentially inside any window spanning the repeat, faking
monomorphism. Surgical single-marker re-windowing rescued four of the
five worst markers on both chips (9%->100%/96%, 27%->99%/96%,
48%->100%/99% vendor concordance), leaving only two markers plus one
sub-amplicon-spaced marker cluster genuinely bad on every window tried —
the true redesign list. The settled production workflow is: first chip
with penalized auto placement, patch the few coverage-forced windows the
cross-run comparison exposes, freeze the curated panel, and run every
subsequent chip with `--reuse-tags`. The curated frozen honeybee panel
scores 93.4%/93.9% vendor concordance on the two chips with 98.0%
replicate concordance.

Re-sweeping window sizes under the new scoring (31/41/61 bp, frequency
floor on) confirmed 41 bp at the knee of the concordance-versus-
information curve: 31 bp genotypes marginally better (92.7% vs 91.8%)
but carries 20% fewer distinct discovered microhaplotypes (110 vs 137),
while 61 bp adds essentially none (138) and loses two concordance points
— on this chemistry, indel exposure destroys the extra flanking
variation as fast as a longer window adds it. Curation matters more than
size: the curated 41-bp panel outscores every uncurated size.

### External benchmarks on public datasets (GT-seq, rhAmpSeq)

Public datasets extend the validation record beyond in-house data, with
intake converters archived beside each dataset (the Flex-Seq capture
benchmark is documented in the 0.9.0 section below).

**Pacific oyster AgriSeq** (PRJNA1118744; 372 Ion S5 samples, 591
amplicons of median 170 bp, design published as Ion AmpliSeq Region +
Hotspot BEDs converted against the oyster_v9 genome with zero REF
mismatches): with the honeybee-derived defaults unmodified, 88.4% of
144.6M raw reads normalized and genotypes agreed with the vendor's
Torrent VariantCaller at **94.52%** (143,383/151,693 cells across 371
name-matched samples, depth >= 20), 97.85% on the 460/542 markers above
90% — within half a point of the curated in-house honeybee panel, on a
foreign lab, species, chip generation, and target length. One panel-size
effect appeared on schedule: 12-bp assignment anchors collide at 591
markers (the birthday bound predicts it), resolved at the time with
`--anchor-length 13`. Since 0.8.3, normalize-reads assigns anchors per
locus instead (the same allocator infer-design uses): every marker keeps
the `--anchor-length` minimum, only colliding prefix groups grow, and
windows identical all the way to the variant share an anchor with a
surfaced warning (their reads resolve by full-tag identity downstream).
A 17,500-marker genome-wide stress panel measured the necessity — 111
colliding 12-mer groups, far above the random-sequence bound — and runs
in ~820 MB with 17,291 markers at the minimum and 209 grown to 13-32 bp.

**Delta smelt GT-seq** (Zenodo 10.5281/zenodo.19069550; 96 PE100
samples, 410 loci, and — uniquely — the panel's real primer file):
GT-seq R1 reads are primer-started, so this platform needs no
normalize-reads stage and exercised the primer-informed `design` route
on non-DArT data for the first time, unchanged. Genotypes from the
R1-only MADC agreed with the archive's independent BWA/bcftools
reference VCF at **98.80%** (33,928/34,340 cells; 99.72% on the
well-behaved core; per-sample median 98.9% across all 96 samples) — the
pipeline's best platform result, on indel-free Illumina chemistry. The
same dataset settled the primer-free question: `infer-design`, run with
the primer file withheld, recovered **368/368 loci's Ref and Alt tags
exactly** as the physical primers define them — a 100% inference score
against ground truth rather than a vendor report.

**Vitis rhAmpSeq** (PRJNA281110, Cornell; the published 2,000-marker
collinear-core panel, NextSeq PE150, 200-300-bp amplicons): the panel
was designed on the collinear core of multiple *Vitis* genomes, so many
markers have no exact home on any single reference — intake is
**reference-free**: each marker's target frame is the consensus of its
own primer-anchored read pileup. An R1-only first pass recovered 1,405
markers at 63% read assignment; the paired-end route then measured the
geometry (94.8% of pairs overlap and merge — the platform's own
published pipeline merges with FLASH, independently confirming
merge-first) and rebuilt consensus from full merged amplicons: **1,759
markers (+25%), frames of median 228 bp, up to three 110-bp tiles per
amplicon** (2,951 markers, 1,030 multi-SNP amplicons), 79.3% read
extraction with `covered` routing. The tiled windows capture ~90% of
all segregating sites (production's full-insert haplotypes span ~174 bp
and capture all sites plus cross-tile phase; the trade buys fixed-length,
frozen, cross-run-comparable MADC haplotypes). The merger itself was
validated against FLASH 1.2.11, vsearch, and fastp on a four-regime
exact-truth pair set — equal to the best (vsearch) in every regime,
zero false merges, and on real NextSeq data the only inter-tool
differences are equal-quality overlap ties, where a per-locus pileup
jury rules our R1-preferring pick correct 62.6% of the time.

No public genotype matrix exists for these accessions, but family
structure is its own truth: progeny of a biparental cross can carry at
most four haplotypes per locus (two per parent), each individual at
most two. All **eight biparental families in the submission (763 runs,
~17.5 GB) were processed family-by-family under a 56 GB disk budget**
against one frozen 2,951-window panel (`--reuse-tags`), so haplotypes
are identical objects across families:

| family | progeny | informative loci | <=4 family haps | <=2 per progeny | calls within inferred parental 4 |
|---|---|---|---|---|---|
| MM531 | 165 | 1,324 | 99.17% | 99.17% | 99.94% |
| MM532 | 166 | 2,055 | 99.22% | 99.29% | 99.92% |
| Y533 | 166 | 2,632 | 99.51% | 99.66% | 99.98% |
| Y519 | 49 | 2,181 | 99.54% | 99.62% | 99.99% |
| Y312 | 44 | 1,852 | 99.57% | 99.62% | 99.98% |

(Three additional tiny families, n = 2-7, are consistent at 99.5-100%.)
Across ~14,900 locus-family tests spanning eight independent crosses
and three sequencing batches, **99.92-99.99% of all haplotype calls in
every family are alleles of that locus's four inferred parental
haplotypes** — at three nines, the microhaplotypes the MADC reports in
763 offspring are the parental alleles Mendelian transmission requires.

### Read accounting

The QC JSON now reports where every read goes. For the default run:
154,584,628 FASTQ records; 88,126,464 assigned to a locus at full length;
87,373,009 (99.1% of assigned) retained in emitted rows; 753,455 dropped with
their 558,194 low-evidence sequences (mean 1.3 reads each); 3,294,541 moved
between sequences by within-sample denoising and 103,628 by cross-sample
denoising (moved reads stay in the table under their centroid). Per-sample
`reads_in_emitted_rows` joins the existing per-sample QC. Twelve samples have
fewer than 100,000 assignable 81-bp tags (lowest: `4404075` with 23); these
FASTQs are dominated by short or primer-less reads and downstream
completeness should be reviewed for them. (Those figures are the 0.5.0 run;
since 0.10.0 the pooled pass folds same-allele sequences only —
`reads_moved_by_global_denoise` keeps its name — and the ledger gains
`reads_assigned_gapped`, `reads_n_mapped` and `reads_n_dropped`.)

### Historical exploratory analyses (0.3.x)

The `results/v03_*` files document the audits that motivated the redesign:
characterization of generated-only rows, right-primer support, cross-extract
reproducibility of paired leaf/root samples, and a diagnostic BLAST against
the public NSP306 *I. trifida* reference following the open `dartag_alleles`
pattern (best hit at the designed locus for most high-support extra rows;
off-target best hits flagged for paralog review, without automatic exclusion
because cultivated sweetpotato is hexaploid and the reference is diploid).
Those analyses used the 0.3.x quality-gated outputs; their row lists are
superseded by the evidence model above, but the biological conclusions —
filter the low-support tail, keep reproducible on-target microhaplotypes,
flag multimappers — are what the current defaults implement.

### Evidence about prior FASTQ filtering and tag length

Every vendor `AlleleSequence` and `ClusterConsensusSequence` is exactly 81
bp; the panel design reconstructs a 97-bp primer-to-primer amplicon at all
3,120 loci, so an 81-bp tag is an observed prefix of a longer amplicon. The
two target manifests report 175,642,223 `tagcounttotal` observations versus
154,584,628 delivered FASTQ records; applying each manifest's `overreppct`
predicts 18,989,753 removed overrepresented reads, and per-sample predicted
and observed deficits are strongly correlated. The delivered FASTQs were
therefore already filtered before release, which supports the preprocessing
hypothesis without identifying the exact vendor commands. A raw-data audit
finds 17,334 of the vendor's 18,677 pure-ACGT rows as exact 81-bp FASTQ
prefixes; the vendor also reports 569 Ref/Alt rows with curated IUPAC
ambiguity (their concrete groupings are not derivable from the A/C/G/T probe
design) and 73 Match/Other rows containing `N`. All 6,238 81-bp prefixes of
the public Breeding Insight PanelHub 109-bp design sequences agree exactly
with this pipeline's independently reconstructed tags.

## Output format and HapApp compatibility (0.8.2)

Statistics columns the pipeline does not compute (`OneRatioRef` through
`AvgPIC`, and `RatioAvgCountRefAvgCountSnp`) emit `.` rather than a blank
cell in every data row (with `--annotations on`, the default, `PICSnp` and
`AvgPIC` instead carry the `Variants` and `MarkerType` annotations — see
the 0.10.0 section). This is a 0.8.2 format change: byte comparisons
against MADCs generated by earlier versions will differ in exactly those
nine columns and nothing else.

Three `generate` flags reproduce the raw-MADC framing that
preamble-expecting consumers such as HapApp and HapApp_utils rely on:

- `--dart-preamble` prepends the vendor-style seven starred rows above the
  header: service/project label, plate barcode, plate number, plate row,
  plate column, plate name, and the targetid row, which carries the
  FASTQ-derived sample ids directly above the header row's sample names —
  the same fastq-id-over-sample-name adjacency as vendor files. Every
  preamble cell without a value is `*`, exactly as in vendor output, with
  one exception: the file's very first cell holds a full provenance string
  (pipeline, version, and every generation parameter as `key=value`,
  prefixed with `* ` so consumers that skip preamble lines by their
  leading star still skip it), so the file records how it was made even
  when separated from its QC sidecar. No timestamp is included — output
  stays byte-reproducible.
- `--preamble-project` sets the label repeated across the first preamble
  row's sample cells (default: the output filename before its first dot).
- `--sample-lookup` takes a two-column CSV/TSV (FASTQ-derived sample id,
  desired sample name; a header line is detected and skipped) and replaces
  the output column names while the targetid row keeps the FASTQ ids.
  Validation is deliberately strict, because a rename table with holes
  silently mislabels count columns: every run sample must be mapped, ids
  matching no sample are errors, and duplicate names are errors.

`read_madc` round-trips preambled output: the header is found by its
`AlleleID` cell, and a numeric, unique targetid row is adopted for sample
identity exactly as when reading vendor files.

## Paired-end reads and multi-window tiling (0.9.0)

The single-window model reads one excised window per locus, which caps
usable data whenever the sequenced fragment is larger than one read: on
a dual-probe capture platform (Flex-Seq, ~300-bp probe-to-probe
fragments, NovaSeq 2x150) only ~52% of reads could ever overlap the one
window, by geometry rather than by quality. 0.9.0 adds the two missing
pieces — paired-end intake and multiple windows per locus — as explicit,
separately validated mechanisms.

**Paired-end intake is one explicit decision, never auto-detected.**
Measured read geometry decides it:

- **Overlapping pairs** (insert shorter than R1+R2, e.g. GT-seq PE100
  over 130-200-bp amplicons): merge first with the `merge-pairs`
  subcommand, then feed merged fragments to `normalize-reads` or
  `generate`. The merger reverse-complements R2, scans overlaps longest
  first, and handles both geometries — the merged read is always R1's
  head plus the overlap plus R2's tail, so in the extended geometry that
  is the whole fragment, and in engulfed read-through R1's adapter tail
  past the fragment end falls outside the overlap and is discarded by
  construction (no adapter trimming step) while R2 bases past R1's end
  are kept when the reads differ in length; disagreeing overlap bases
  take the higher-quality call, `--unmerged r1|pool|drop` states what
  happens to the remainder, mates are paired by read ID (out-of-step
  names or unequal counts are a hard error), and samples are independent
  so `--jobs` merges them in parallel processes. On a synthetic paired-end truth set with modeled
  adapter read-through, 100% of 513,715 pairs merged and molecular count
  recovery was 98.6% — one fragment, one count.
- **Non-overlapping dual-anchored designs** (Flex-Seq: R1s anchor the
  left probe, R2s the right, mates cover opposite halves of a fragment
  neither read spans): pool both mates as independent reads. Merging is
  geometrically impossible there and pooling cannot double-count a
  window, because no window is covered by both mates.

**Multi-window tiling** puts several markers on one physical locus: an
optional `LocusGroup` column (probe-design TSV and reference-tags CSV)
declares markers as tiles of one target region. Tiles must share an
identical Ref expansion frame (validated at load), share one k-mer
territory, and reads locate to the group before routing to member
windows. Panels without the column behave exactly as before, to the
byte.

**Routing between a group's windows** is `--group-routing`:

- `nearest` (default): each read feeds exactly one tile — the fully
  covered window whose center is nearest the read's center.
- `covered`: the read's excision is emitted at every member window it
  fully covers. Each window still sees each fragment at most once in
  both modes; the QC separates the extra per-window emissions
  (`tile_extra_emissions`) from the read-fate accounting.

The distinction matters because merged full-span fragments all share
the amplicon's center, so `nearest` degenerates to winner-take-all: on
the GT-seq benchmark, one tile took effectively all depth at 163 of 164
tiled loci and the other went dark. Use `covered` whenever reads span
more than one tile (merged fragments, or long reads over closely spaced
tiles); a fragment that covers two windows genuinely observes both
SNPs.

**Distinct-window guard.** Tiles share a coverage profile, so fresh
window placement can collapse two members onto the same start — and two
markers with identical Ref tags make every ref-carrying read an exact
tie, which the assignment layer drops as ambiguous. Measured on GT-seq,
that silently zeroed a locus's ref depth and flipped its calls from
100% to 3% concordance. Fresh placement therefore forces member windows
to excise distinct sequences (a colliding member takes its best-ranked
distinct start; the nudge is listed in QC `locus_groups.window_dedup_nudged`;
no distinct start at all is an error), and `--reuse-tags` refuses a
frozen panel containing twin windows instead of replaying the defect.

**Window freezing is load-bearing for paired-end modes.** A four-way
GT-seq experiment on identical reads and truth (96 samples, 409
original + 206 added tile/rescued markers, DP >= 20 against the
archive's independent BWA/bcftools VCF):

| variant | original markers | added markers | scored cells |
|---|---|---|---|
| R1-only | 98.80% | — | 34,340 |
| pooled mates + fresh windows | 95.18% | 87.92% | — |
| merged + fresh windows | 92.99% | 90.06% | — |
| merged + frozen tags + tiles, `covered` | **98.80%** | **98.43%** | **52,971** |

Every fresh-placement regression was window drift (all 368 windows
moved; the worst landed on twin-window ambiguity), and the recipe that
holds is: `merge-pairs`, then `normalize-reads --reuse-tags` with the
proven panel plus tile rows, then `generate` — original markers
reproduce the R1-only result on the same cells while the tiles add 54%
more truth-validated genotype calls at parity.

On the dual-probe capture benchmark (cranberry Flex-Seq, 24 samples,
17,242 reconstructed clusters scored against the population truth VCF),
the same mechanisms lift the geometry ceiling: R1-only single-window
92.39% on 32,675 cells with median 7,733 loci recovered per sample;
pooling both mates raises it to 93.24% on 57,577 cells and 13,421 loci
(the published vendor pipeline recovered 14,550 of 17,502); a
three-tile 41,033-marker panel with `covered` routing reaches 82.9%
read extraction and **94.59% on 120,261 cells** — 3.7x the callable
genotypes of the single-window R1 run, at higher concordance. A
subsequent full-cohort campaign (192 samples, K=4 tiles at 45-bp
spacing with 110-bp windows, one cohort-level `generate`) raised
capture of the population's segregating-site catalogue from 64.5% to
78.4% while *improving* trio incompatibility from 4.5% to 3.6% —
window geometry solved the coverage gap (outside-window misses fell
from 18.6% to 1.9% of sites) and cohort-scale evidence produced a
cleaner haplotype catalogue, not just a larger one.

**Adapter trimming is not required, and adapter cannot reach the
emitted haplotypes.** Three structural guarantees, one per read mode:
windows are defined in frame coordinates and bounded by the frame
(`window start + length <= frame length`), so on amplicon platforms —
where the frame is the amplicon — a window cannot reach past the
biology into adapter by construction; merge mode removes read-through
adapter geometrically (the engulfed merge drops R1's tail past the overlap;
proven on the synthetic paired-end set with modeled read-through: 100%
merged, zero length errors); and pooled mode exists precisely for
long-fragment designs whose reads end before the fragment does, so
there is no read-through to leak. The residual theoretical path
(capture platforms whose genome-slice frames outrun a short fragment,
on the unmerged tail) funnels into identity-gate rejection, because
variable shear junctions never reproduce as one sequence. The census
agrees: across 127,667 emitted rows in six MADCs spanning five
platforms (GT-seq at two window lengths, rhAmpSeq pooled and
family-scale, the 41K-tile capture panel, Ion AgriSeq), zero rows
contain TruSeq or Ion P1 adapter at a strict 13-mer probe. (A loose
10-mer probe flags one capture locus whose genome happens to contain
`AGATCGGAAG` diverging from TruSeq at base 11 — a chance genomic
k-mer inside the reference tag itself, not contamination.) Measured on
the capture benchmark, trimming adapters before the run changed dosage
concordance by 0.05 points (92.24% untrimmed vs 92.19% trimmed) while
costing 2 points of extraction from the shortened reads.

## Choosing the window length (`window-scan`)

Tag/window length was historically chosen per platform against external
truth (the honeybee 41-bp knee came from a concordance-versus-length
sweep against vendor calls). `window-scan` makes the same decision from
the project's own reads: one coverage pass predicts **retention**
(fraction of located reads fully covering the best window) for a dense
grid of lengths, and each candidate length then runs the real
normalize-and-generate machinery on a sample subset to measure
**cleanliness** (assigned reads not absorbed by denoising or dropped as
failed evidence — per-base error and indel exposure accumulate with
length) and **information** (evidence-passing novel haplotype rows per
100 markers — rises with length until the local linked variation is
exhausted). The recommendation is **survival-first**: survival =
retention x assignment x cleanliness; candidates within `--info-band`
(default 10%) of the best survival are eligible, and the most
informative eligible length wins, ties breaking short. Information
deliberately cannot lead the decision, because it can be inflated by
exactly the rows one fears: systematic errors (Ion homopolymer indel
families recur across samples) pass the reproducibility route and
masquerade as haplotype diversity. An earlier score that weighted
information directly followed that inflation long — on Ion oyster data
it recommended 81 bp against the externally validated 41-bp knee — and
the survival-first rule was the correction. Output: a TSV of the
curves, a QC JSON, and per-length run directories; a candidate the
distinct-window guard refuses (tiles inseparable at that length) is
reported failed rather than aborting the sweep.

Under the survival-first rule the scan reproduces every externally
arbitrated choice: Ion oyster 41 bp (assignment falls 15 points from
31 to 81 bp as indel-bearing reads fail the substitution-only identity
gate — the Ion damage shows in assignment, not cleanliness), GT-seq
81 bp (the frozen-validated length; the 96-sample truth test measured
110 bp at only -0.08 points, so the band's conservatism is cheap), and
capture 61 bp. On merged rhAmpSeq amplicons — clean chemistry,
fragments dwarfing every candidate, information unsaturated through
150 bp — the rule recommends 81 while the family instrument cleared
110 and the platform's own pipeline uses full inserts: where richer
haplotypes are the goal and an external or genetic arbiter is
available, widening `--info-band` expresses that choice explicitly.

**Seed sensitivity (0.11.0).** The scan runs on a sample subset (six files
by default), and a recommendation can be a property of that subset rather
than of the panel. `--scan-draws N` repeats the whole sweep on N
independent seeded draws (seeds `--scan-seed`, +1, +2, ...) into
`draw_1/`, `draw_2/`, ... and reports each draw's recommendation, the
votes per length, and their agreement (`recommendation_stability`);
the consensus is the most-voted length, ties breaking short. A
recommendation that flips between draws is the signal to scan more
samples or to arbitrate externally — the one check a user without truth
can make on the length decision itself.

Falling cleanliness deserves a precise reading: a churned read is one
whose exact sequence the denoiser recognized as an error neighbor and
folded back into its parent — the count survives, no error row is
emitted, so churn measures error-handling load, not leakage. The
narrow leak path is a systematic error (same position, same
substitution, recurring across samples) clearing the reproducibility
route and emitting as a false Match row — and the scan's information
metric cannot distinguish that from real linked variation, so the scan
is a screening instrument whose recommendation still wants an external
or genetic arbiter. On GT-seq the arbitration is done: the scan
recommended 110 bp over the historical 81, and against the archive's
independent truth VCF on all 96 samples (identical merged reads,
identical fresh placement), 110 bp scores 98.62% versus 81 bp's 98.70%
— at most 0.08 points of dosage cost while emitting 27% more Match
haplotypes (1,568 vs 1,231 rows). The information the scan saw is not
corrupting target genotypes; whether each added Match row is linked
biology versus a systematic error family is exactly what the family
instrument settles (on rhAmpSeq at 110 bp: under 1% family-bound
violations across 763 progeny). The same experiment refined the
window-drift picture: fresh 81-bp placement WITH the distinct-window
guard scores 98.77/98.57 (old/new markers) versus the pre-guard 92.99/
90.06 — the drift "catastrophe" was entirely the twin-window defect,
and honest residual drift costs only 0.03-0.14 points against frozen
windows. Freezing remains mandatory for cross-run haplotype identity;
it is no longer load-bearing for accuracy.

## General indel recovery and review fixes (0.10.0)

Two independent code reviews and a colleague's finding on a DArTag project
— genome assemblies predicted indel-bearing microhaplotypes (1–13 bp,
often several events per tag) that the vendor MADC never reported, and
neither did 0.9.x — motivated this release. The 0.9.x identity gate was
substitution-only, so an undesigned indel frameshifts the tag against
every canonical sequence and the read was dropped as off-target; on that
project only 6 of 447 read-supported indel microhaplotypes fell inside the
envelope, while 423 lay within the same 20% budget once gaps count as
edits.

**One general mechanism, not a catalogue of cases.** When every primer
candidate rejects a read, `select_marker` offers one gapped reading: a
start-anchored banded alignment (band = `--max-indel-length`, default 10)
whose end may overhang by up to a band on either side (a read tail past
the reference is a deletion's signature, a reference tail past the read an
insertion's), scored under the same edit budget as substitutions. The
reading is believed only when it is composed like an indel haplotype —
at most 6 substitutions and 6 gap events alongside the gaps — because a
co-amplified paralog reaches the budget through many substitutions (the
captured reads on the motivating project carried a median of 13
substitutions per read against 2 for true indel haplotypes, and one locus
would otherwise have swollen from 351 to 38,653 reads). Accepted reads are
counted as ordinary fixed-length rows in the vendor's shifted-suffix
convention; the counting layer, evidence rules and denoiser are unchanged.
`classify_sequence` reads the target base at the position the *cheaper*
explanation puts it: the fixed offset unless a gapped reading is strictly
cheaper, well composed, and verified by at least 12 bases after its last
gap (clustered errors at a read end must not masquerade as a shift); a
target inside a deletion or pushed past the read is `Other`. This also
corrects in-envelope indel reads that 0.9.x classified at a shifted offset
— the delta smelt locus `NC_061083`, whose real 1-bp insertion scored 26%
against the truth VCF, is now concordant. Cost is bounded by an exact
lower bound: a Myers/Hyyrö bit-vector edit distance (C-speed integer
operations) proves most hopeless tags cannot meet the budget before any
alignment runs, so a 15-sample DArTag run takes 75 s against 41 s
(a prototype without the prescreen took 150–290 s). `--gapped-rescue off`
restores the 0.9 assignment path; QC reports `reads_assigned_gapped`.

**Review fixes adopted.** The pooled denoising pass now folds an
evidence-failing sequence only into a protected parent carrying the same
target allele (0.9.x moved reads across alleles); parent ranking is
closest-first, then deepest; an observed N counts as a mismatch and
N-bearing tags are credited to their unique compatible sequence or dropped,
never emitted; a k-mer repeated at two positions inside one target frame no
longer anchors reads (it kept only its last position); `merge-pairs` keeps
R2's extension when the reads differ in length and errors on mates out of
step; `window-scan` takes `--scan-seed` for a random sample subset and
`--generate-args` to forward evidence flags to the candidate runs. The
per-sample cache format is 0.5.0, so existing caches rebuild once.

**Considered and not adopted, with the measurement.** Collapsing reads onto
IUPAC-coded canonical tags: DArT uses those codes to encode real secondary
SNPs, and on the blackberry panel the collapse lost 272 vendor-matching
rows and moved count agreement from 0.982 to 0.976 — so rows stay exact
and ambiguity-coded canonicals keep reporting zero exact counts. A
synthesized-primer mask (worse still on blackberry), an unconditional
merge of indel families in long homopolymer runs (it absorbed 66 real
indel alleles on the DArTag project at its Ion-safe default), a
hard-coded ubiquity filter, and a single-gap-only indel tier with its
downstream-context columns were all left out in favour of the one
mechanism above.

**Validation, identical inputs, 0.9.1 → 0.10.0 — every platform in the
corpus.** Delta smelt GT-seq vs the BWA/bcftools VCF 98.80% → **99.49%**
R1-only and 98.80% → **99.49%** on the merged tiled panel (core 99.80%;
358/361 markers ≥90%); oyster AgriSeq vs Torrent VariantCaller 94.51% →
**94.91%** (markers ≥90%: 460 → 472); honeybee AgriSeq chip 2 vs the
vendor caller 91.60% → **93.84%**, replicate concordance across 30 pairs
97.64% → **98.49%**, depth per sample-marker +28%; cranberry Flex-Seq
(24-sample subset, pooled mates, frozen K4 panel) vs the truth VCF 94.50%
→ **94.78%** on 3.8% more callable cells; vitis rhAmpSeq MM531, 165
progeny, 193,737 → **202,811** calls (+4.7%) at within-parental-four
transmission 99.890% → **99.894%**, and the 16-sample head-to-head against
`analyze_amplicon` 93.0% → **94.9%** exact agreement at 99.953% Mendelian
consistency; the three in-house DArTag projects `validate`d against their
vendor MADCs kept every substitution-side number (sweetpotato count
Pearson 0.9992 → 0.9993, red clover 0.9733 → 0.9733, blackberry 376
samples 0.9814 → 0.9814) while adding 1,544 / 2,393 / 8,313 indel-type
rows respectively — rows the vendor pipeline structurally cannot report —
and the private DArTag project (15 checked samples): 2 → **348 of 447**
assembly-predicted indel microhaplotypes recovered at 97% read recall,
including 72–82% of those with two or more indel events, while the new
heterozygotes those rows create have allele fractions indistinguishable
from vendor-agreeing heterozygotes and a third are confirmed by the same
sample's genome assembly. Indel-rich loci gain reads (about 60 markers
more than doubled on that project, the same set two independent indel
implementations found), and Ion row counts rise about 9% (5,890 → 6,443
on oyster) as slippage haplotypes are now assigned rather than dropped;
their target base is still read correctly, which is why SNP accuracy
rose. Two vendor-agreement effects are expected and explained: DArT emits
N-bearing rows that 0.10.0 never does (29 of the 30 sweetpotato rows that
"disappeared"), and at indel-shifted rows inside the substitution envelope
the vendor classifies at the fixed offset while 0.10.0 reads through the
gap (153 sweetpotato rows, mostly vendor `Other` becoming Ref/AltMatch or
Ref/AltMatch swapping) — the smelt truth VCF settled that question in
favour of the gapped reading. Where an indel allele carries the alternate
base, the recovered rows create heterozygotes the vendor calls homozygous
(red clover: updog dosage agreement with the vendor MADC 97.86% → 95.95% at
DP≥10, 13,797 new indel-driven heterozygotes). The 152-progeny red clover
mapping family with both parents adjudicates them: with every MADC pushed
through the same updog (diploid, norm) fit, the Mendelian-violation rate is
1.59% under the vendor MADC, 1.61% under 0.9.1 and **1.30% under 0.10.0**
on 2% more tests (parents DP≥20, progeny DP≥10, MPP≥0.9) — the new
heterozygotes are inherited, and the vendor's disagreement is the vendor
missing real alleles.

**Band width.** `--max-indel-length` 10 is the validated default. On the
private DArTag project's full 558-sample order, band 10 recovered 374 of
the 447 read-supported assembly-predicted indel microhaplotypes (84%; the
extra samples let 26 more clear the reproducibility route than in the
15-sample run) and band 20 recovered 393 (88%) by adding exactly 76 rows —
all 11–20 bp indel haplotypes, median ~5,600 reads each — with the same
row-plausibility profile, vendor SNP concordance within 0.04 pp
(97.87% vs 97.83% over 1.5M cells) and equal runtime. That order also
carries a 264-progeny F2 family and three genotypes in technical
quadruplicate, which arbitrate the rows directly: indel-haplotype calls
fall inside the founders' four haplotypes at 99.86–99.96%, the same rate
as substitution haplotypes (99.91–99.94%; the vendor's own rows 99.95%),
and an indel haplotype called in one replicate is re-called in its twins
98.4% of the time against 99.3% for substitution haplotypes — inherited
alleles, with the replicate gap consistent with indel alleles amplifying a
little less efficiently. Widen the band on
panels where long indels are expected; the composition guard, not the
band, is what keeps paralogs out.

**Annotations: `MarkerType`, `Variants`, and the variants sidecar
(`--annotations on`, the default).** Two vendor statistics this pipeline
never computes are repurposed, header included, and a sidecar is written
next to the table; `--annotations off` restores the exact vendor header
and writes no sidecar. The column positions are unchanged, so positional
consumers see the same 16-column layout.

- `AvgPIC` becomes **`MarkerType`**: what the row carries *beyond its
  canonical*, derived from its Variants. The Ref row is `.`; the Alt row
  is labelled by the design itself (`snp` for a designed SNP, `indel` for
  a designed indel); every other row by its events once the marker's
  designed events are removed: `snp` (substitutions only), `indel` (gaps
  only), `both`, or `.` (no concrete event beyond the canonical — a row
  that merely realises an ambiguity code). A downstream user who wants
  SNP-only microhaplotypes drops `indel`/`both` rows by this label alone
  and keeps designed-indel markers intact.
- `PICSnp` becomes **`Variants`**: the row's events against the Ref tag,
  one token per event, `;`-separated, in 1-based coordinates on the tag as
  printed (design strand; a consumer that knows the strand flips them, as
  it already does for substitutions). Substitutions are `91:A>T`;
  insertions and deletions are VCF-anchored on the preceding reference
  base and left-normalised (`55:C>CAAAA`, `70:TAG>T`); then comes the
  window geometry a fixed-length row cannot show: `clip=k` when an
  insertion pushed the last *k* Ref bases out of the window, `overhang=k`
  when a deletion pulled *k* downstream bases in. The Ref row is `.`, the
  Alt row of a designed SNP is exactly the design (`46:C>G`), and an
  AltMatch row lists the designed substitution among its events. A
  position where the Ref tag carries an ambiguity code is never an event
  (a concrete base realises the code), matching the IUPAC handling
  elsewhere.
- **`OUTPUT.variants.csv`** lists the same events one per line —
  `CloneID, AlleleID, AlleleSequence, Position, Ref, Alt, Type` with Type
  `snp`, `ins`, `del`, `clip` (Ref holds the unobserved bases) or
  `overhang` (Alt holds the pulled-in bases) — keyed by CloneID and
  sequence so it survives allele-ID renumbering downstream.

Every row is described by its cheapest reading against Ref — a gapped
reading only when strictly cheaper than substitutions — through the
alignment with the fewest gap openings among equal-edit ones (a
three-state banded table), so an inserted block is one event even where
its bases happen to match the reference, as an affine-gap aligner would
describe it. Rows of a marker whose own Alt differs from Ref by a gap
are always read with gaps, so an AltMatch of a designed insertion that
adds one SNP lists the inherited insertion plus the SNP (and is labelled
`snp`) rather than a run of shifted substitutions; because MarkerType is
derived from the same events, the two annotations cannot disagree.
Canonical tags that share an ambiguity code are scored as matching at
that column by the alignment, its prescreen and the substitution count
alike — the relation `tag_mismatches` uses. Rerunning an output with
`--annotations off` removes a sidecar left by an earlier run. A
microhaplotype with several SNPs and indels is simply a longer
list: phase lives in the row, and a per-site consumer emits one record
per event, which is what `BIGr::madc2vcf_all` does today for
substitutions by realigning each row to Ref. Checks: the Alt-row token
reproduces the design's own `position:Ref>Alt` for 3,120/3,120
sweetpotato, 3,005/3,015 blackberry-design, 2,982/3,000 red clover and
3,030/3,040 private-DArTag-panel loci (every remainder is an Alt tag that
is not a single SNP — ambiguity codes, or a vendor Alt tag that carries
an 8-bp insertion, reported as such); replaying each annotated row's
events onto its Ref tag returns the row's own sequence for all 9,430
annotated rows of the public sweetpotato plate (16,954 rows; 43,766
substitutions, 1,700 insertions, 1,830 deletions), the delta smelt GT-seq
run (1,775 rows) and the private project's 15-sample run (9,430 rows; a
median of two events per annotated row, 41 on the most heavily
substituted one); the annotations add no measurable runtime (the 15-sample
run: 78.7 s off, 76.8 s on). Row counts by
label (`marker_types`) and by event kind, plus the sidecar path
(`variants`), are in the QC summary.

Consumers: HapApp's fixed-allele-ID step keeps columns 1–3 plus the
samples, so these columns do not travel through it to
`BIGr::madc2vcf_all` (the sidecar does not need to); BIGr's raw-MADC
readers (`get_countsMADC`, `filterMADC`) drop statistics columns by
name, so until `MarkerType` and `Variants` are on their list, files
headed there should be generated with `--annotations off` or the two
names added; polyRAD's `readDArTag` and HapApp locate samples
positionally and are unaffected.

## Raw deliveries, tag reconstruction and denoiser scaling (0.12.0)

Three changes made while setting up the alfalfa (DAl21-5779, DAl21-6024)
and pecan (DPec23-8107) panel-paper datasets, none of which alters a table
produced from a preprocessed delivery:

**Raw DArT deliveries.** Older DArT FASTQs keep the inline barcode and the
library's common sequence `GATGGATCCAGTGC` in front of every amplicon
(alfalfa 2021: 12 + 14 bases, then the 54-bp tag, on 83-bp reads; pecan
2023: the same 26 bases, the 81-bp tag, then adapter and poly-G on 138-bp
reads), so nothing is primer-started and a 0.11.0 run assigned no read.
`generate --leading-anchor SEQ` and `infer-design --leading-anchor SEQ`
cut each read where the anchor ends (it must end within 40 bases beyond
its own length) and count reads without it as `reads_no_leading_anchor`;
the anchor is part of the cache key, the provenance line and the QC
parameters (so assignment caches written by 0.11.0 rebuild once, as they
did when the partial-read length joined the key). Both commands check the first files before any read pass:
when almost no read starts with a panel tag but most carry one after a
constant prefix, the run stops and names that prefix (the suffix shared by
the bases in front of the tag at the modal offset, sampled across files so
per-sample barcodes do not survive), and `--leading-anchor auto` applies
the detected sequence instead (`leading_anchor_detected` in the QC). A
preprocessed delivery passes the check silently (red clover: 87% of reads
start with a panel tag, 0.3% carry one further in), and `auto` refuses on
it rather than cutting anything. Measured: alfalfa assigns 5–6% more reads
per sample than the vendor counted, as the other DArTag panels do; the
primer-free route on the raw alfalfa reads, with no primer file, resolves
2,999 of 3,000 loci in the exact tier, all identical to the design tags.

**Design route: a tag carries its variant.** Tag reconstruction from the
primer table and probe design now keeps only the primer-started windows
that contain the bracketed variant. The pecan design has one marker
(Chr04_45131234) whose target carries a duplicated primer-bearing flank;
both placements showed the Ref base at the variant offset and the marker
was reported as ambiguous. With the rule, the pecan design tags match the
vendor's Ref/Alt rows at 3,100 of 3,100 markers (3,084 identical, 16 under
vendor IUPAC codes). Any DArT marker table with `MarkerName` and a bracket
`TargetSequence` serves as the probe design, extra markers included.

**Denoiser scaling.** At the default distance of 1 the eligible parents of
a sequence are its single-substitution neighbours, which `denoise_counts`
now enumerates and looks up once the centroids outnumber those neighbours
(the scan stays cheaper at an ordinary locus with a few dozen centroids;
the choice is made per sequence); the selection is identical (a randomized test checks it against the scan, and
the sweetpotato table regenerates byte-identically from cache, and the alfalfa tables re-stamped with 0.12.0 are byte-identical in body to the ones made with the 0.11.0 denoiser).
The pecan panel has one paralogous locus (Chr04_23494647: reads 10–16
mismatches from the design tags around a GGGGGG homopolymer, about 6,000
unique sequences per sample, 409 rows in the vendor table) that made the
old pass quadratic: its finalisation ran for more than 90 minutes before
being stopped, and completes in minutes now. The locus is not flagged in
the generate QC yet; a per-locus diversity and paralog warning is the
natural next addition.

80 unit tests.

**0.12.1: the finalisation runs across the workers.** `generate`'s
finalisation (protection, per-sample denoising, the pooled pass, the
bookkeeping) ran in one process; on a family of 169 rhAmpSeq runs with
150-bp windows and 637K reads per sample (311K distinct sequences per
sample, 67 per locus, so the lookup rule above never triggers) that stage
took 107 minutes against 32 for the eight-worker assignment. Stage timers
(`finalize: ... in N s` on stderr) showed where it went on a 44-run family:
protection 30 s, per-sample denoise 154 s, pooled pass 1,123 s. Both
denoising passes are independent per unit of work, the per-sample pass per
sample and the pooled pass per locus (its aggregate, the protected set and
the gapped classifier that checks a move keeps the target allele), so
`finalize_counts(..., jobs=N)` now maps each over a process pool with the
panel and protected set sent once per worker. The selection is identical
to the sequential passes (a test asserts it on two loci and five samples
with reads moved by both passes, and the 0.12.0 tables regenerate
byte-identically from cache). Protection, the final bookkeeping and the
table write stay in the main process. The table is unchanged with
`--jobs 1`.

## Truncated reads: partial-read salvage with a locus guard (0.11.0)

A DArTag order whose delivered files were trimmed to 109–113 bp per file,
with several amplicons 1–9 bases shorter than the 109-bp window, showed
what the 0.10 rule "reads shorter than the tag are dropped" costs: two loci
with 101- and 102-bp amplicons had no reads at all in the table (the vendor
counted 104,507 and 69,381), a Y-specific locus lost 46 of 211 carrier
samples because one plate's reads were 106 bp, and in cells where one
allele's amplicon is shorter than the other's every read of that allele was
discarded, so the cell was called homozygous for the other allele (546
Ref reads of 109 bp against 478 Alt reads of 108 bp in one sample).
**Truncated reads support alleles but never invent them (partial-read
salvage, on by default).** Reads shorter than the window — adapter- or
quality-trimmed reads at amplicons a few bases short of the tag, or files
trimmed to fewer bases than the tag — are primer-assigned and then credited
to the one full-length allele that has already earned a row and of which the
read is an exact prefix. Two guards
keep this from moving a call: the read must be longer than the longest
prefix shared between any two learned alleles of its locus (a locus credits
truncated reads for both target alleles or for neither, so a ratio can
never shift, and a read that does not reach the target base is always
ambiguous), and `--min-partial-length` defaults to `--tag-length` minus 9,
so at most nine tail bases are unobserved and primer-dimer fragments never
qualify. Defaults: `--partial-read-mode assign`; `off` restores the drop
behaviour. An `extend` mode that borrowed the design Ref suffix to turn an
unmatched prefix into a provisional allele existed through 0.11.0 and was
removed: on the project above it manufactured a RefMatch row from a 30-bp
primer-dimer. A tag table that carries identical Ref
and Alt tags (presence/absence designs) shares its whole tag between the two
alleles, so no partial is ever credited at such a locus. The QC reports `partial_reads_assigned`,
`partial_reads_ambiguous` and `partial_reads_unmapped` per sample. The
minimum partial length is part of the cache fingerprint, so runs that relied
on the old default of 30 rebuild their cache once.

**Measured on that project (564 samples, 3,040 markers, band 20).** The
default salvage credits 3.04M reads (+0.98%) to 37,359 cells at 979 loci,
makes 15,987 more cells callable for dosage, and raises count agreement with
the vendor on the 24,755 shared rows (Pearson r over nonzero cells) from
0.961 to 0.971, because the vendor counts these near-full reads too; dosage
agreement with the vendor moves from 97.83% to 97.85% on 1.53M cells, and
BIGr `madc_summary` missingness (cells under 10 reads) falls from 7.55% to
6.62% overall against the vendor's 6.74% (per-sample median 7.20% → 6.22%
against 6.38%). The locus guard, run against the same cache without it, removes 504,677 of the
salvaged reads from 7,914 cells at 269 loci — loci with a learned Match row
that differs from Ref or Alt only beyond the truncated length — and reverts
469 of the 472 dosage-class changes that unguarded crediting introduces.
The 2,316 cells (0.15%) whose dosage class differs from the no-salvage table
are the allele-specific amplicon-length cases above, which the no-salvage
table calls homozygous; two independent instruments say the new calls are
right: their haplotypes fall within the locus's four founder haplotypes in
the F2 families at 99.8% (GEN, n=1,067) and 100% (DV, n=558), the F2 rate,
and re-call in technical replicates at 98.6% against 98.5% for the calls
they replace, while the vendor's own call agrees with the new call in 1,276
of 2,164 cells and with the old one in 737. Reads that stop before the
target (a 35-bp junction product at one Y-specific locus, 90-bp amplicons
whose target sits beyond them, primer-dimers) stay out: they are below the
length minimum, and would be ambiguous between Ref and Alt even above it.
The assign stage adds under a minute to a 50-minute run. The preamble now
records `min_partial_length=` beside `partial_read_mode=`.

**Measured on three other DArTag panels (salvage on against off; same code,
scripts and band 10).** Sweetpotato (80 samples, 81-bp tags): +40,438 reads
(0.04%) at 304 loci; 5 of 196,494 called dosage cells change class and the
vendor's own call agrees with the salvaged class in all 5; leaf/root
replicate re-call 90.11% either way; count correlation with the vendor table
0.9993 → 0.9997. Red clover (376 samples, 109-bp tags, primer-free panel):
+941,538 reads (1.3%) at 823 loci, 8,733 more callable cells, 786 of 717,768
called cells change class (the vendor sides with the salvaged class in 714,
with the unsalvaged in 14); dosage agreement 96.87% → 96.98%, count
correlation 0.9765 → 0.9845. Blackberry (376 samples, 109-bp tags):
+1,552,191 reads (0.76%) at 634 loci, 6,160 more callable cells, 2,112 of
1,006,967 change class (1,994 against 50); dosage agreement 95.18% → 95.39%,
count correlation 0.9823 → 0.9867. Rows are identical in every pair: a
partial read never establishes one. The credited reads sit at loci whose
amplicons run short of the tag, where the counts land on the vendor's totals
(one red clover locus: 2,606 → 83,324 reads against the vendor's 83,303);
loci whose amplicons are shorter than the partial minimum (78–84 bp at the
largest red clover examples) or whose reads end before the distinguishing
base stay uncredited, by design.

**One reading per row (0.11.0).** 0.10.0 read a row twice: classification
confirmed a gapped reading on a unit-cost traceback with the composition
guard and the 12-base verified suffix, while `Variants`/`MarkerType` took
the fewest-openings alignment unconfirmed, so a row could list an indel its
category ignored, or be classified through a gap its annotation placed
elsewhere. Both now use `row_reading`: the cheaper of the two
substitution-only readings is the bound, the fewest-openings alignment is
the reading, and it is accepted only when verified by GAP_VERIFIED_SUFFIX
bases after its last gap and well composed — composition being judged on
the equal-edit explanation that needs the fewest substitutions
(`fewest_substitution_alignment`), because the fewest-openings explanation
trades a gap opening for a substitution wherever it can, and against the
canonical the row is read from, so a designed substitution never counts
against an AltMatch row's indel reading. A reading that fails is not a
reading of the row at all: the row is classified at the fixed offset and
described as substitutions, by both consumers. On the project above this
re-labelled 76 of 28,793 rows (53 `Other` rows became RefMatch/AltMatch,
16 swapped between RefMatch and AltMatch, 7 became `Other`) and rewrote the
annotation of 554 rows whose listed indel sat inside the last 12 bases or
whose composition failed the guard (`both`/`indel` → `snp`); counts moved by
97 reads in 66 rows through the category-aware pooled denoise.

**Freeze health (0.11.0).** Freezing a window panel (`--reuse-tags`) is
what keeps haplotype identity comparable across runs, but read geometry
drifts run to run and a frozen window can slide off the reads a new run
delivers. `normalize-reads --reuse-tags` therefore still runs the coverage
pass fresh placement would use and scores every frozen window against it:
the QC's `markers` table carries each window's predicted full coverage on
this run beside the window fresh placement would choose today and its
coverage, and `freeze_health` summarises the read-weighted coverage of
the frozen panel against the fresh one, the markers losing more than ten
points, and the ten worst. One line on stderr says what freezing costs
this run; `--window-placement center` skips the pass. The excised FASTQs
and the tag table are unchanged by the check.

**Vendor-only markers, re-derived.** Of the 53 markers in the tag table that
had no reads in the 0.10 table, 46 carry vendor counts (40,101 called
cells): 33 markers (29,861 cells) have only primer-bearing fragments that
end before the target (the vendor counts them as Ref), 10 (8,478 cells)
have 44–77-bp amplicons that do contain the target — recoverable only by a
per-locus window, not by salvage — or, for two of them, reads that are
largely absent from the delivered files, 2 (1,762 cells) have 101/102-bp
amplicons and are counted since 0.11.0, and 7 have no reads anywhere and no
vendor counts. The 46 markers whose vendor Ref and Alt rows are identical
are a separate class: their target sits at position 110, one base past the
109-bp window (the vendor's SNP report lists them without a position, and
the 110th base segregates in the plates whose reads are longer than 109 bp).

## Genome paralog filter (`filter-paralogs`, unreleased)

On the public lettuce 3K DArTag data about two thirds of the reads come from
paralogous copies. A paralog read that carries the target base becomes a
RefMatch/AltMatch row, and `BIGr::madc2vcf_all` sums it into the target SNP,
so inbred lines look heterozygous. Raising `generate --min-identity` removes
the diverged copies but not the close ones. `filter-paralogs` is an opt-in
cleaner for such panels. It reads any MADC, from MADCATS or the vendor, with
or without the preamble, so it also cleans vendor deliveries, and HapApp's
fixed-allele-ID output (`Code_version`, `AlleleID`, `CloneID`,
`AlleleSequence`, samples; IDs such as `|Ref_0001`; merged counts written as
`12.0`). The command lives in `madcats_paralogs.py`, which imports the MADC
helpers it shares with `validate` from `madcats`; `build_parser` imports it
when the parser is built. Run as `python3 madcats.py`, the file is
`__main__`, so it registers itself as `madcats` before the filter module
imports it; otherwise the module would load a second copy whose `MadcError`
`main` does not catch.

Design:

- **Aligner.** minimap2 through `mappy`, an optional extra imported only by
  this command; MADCATS still has no required dependency. The index is a
  standard `-x sr` index (k=21, w=11). At map time single-seed chains
  (`min_cnt=1`, `min_chain_score=20`) and 50 hits (`best_n`) are allowed.
  mappy reports neither AS nor XS, so each hit's DP score is recomputed from
  CIGAR, NM, matching bases and block length with the pinned short-read
  scoring (+2/−8, gaps 12+2 or 24+1, ambiguous −1). Ties are exact.
- **Marker placement.** Each marker's Ref tag is aligned, and every
  equal-best site the search returns is kept. The CloneID `CHROM_POS` is
  read as a position on the genome given: it picks among equal-best Ref-tag
  sites and stands in (±tag length) when the Ref tag does not align. When no
  equal-best site covers it, the Ref tag is realigned to a k=11 index of the
  window around the position (±tag length + 20 bp), and an alignment
  covering the position becomes the marker's site (`clone_id_site`), even
  when the tag scores better elsewhere; without one the marker is
  `unplaced`. `--clone-positions ignore` turns these uses off, for a genome
  other than the one the CloneIDs refer to. On the lettuce vendor MADC
  against v8 it changes 11 markers: the 9 fallback windows become
  `unplaced`, keeping 1,032 rows labelled `paralog`; 1 tie stays
  multi-copy; and `Chr2_001081815` goes back to Chr4, dropping 11 rows that
  `use` keeps. The id95 labels do not change.
- **Labels.** Each haplotype gets its best score at its marker's site and
  its best score elsewhere: higher at the marker is `on_target`, higher
  elsewhere is `paralog`, a tie is `ambiguous`. Every haplotype is also
  realigned to a k=11 index of its own marker's sites (±20 bp), and the site
  score is the higher of that and the genome-wide hit there. Without that
  realignment, 18% of id95 haplotypes had no hit at their own marker, and
  326 that bwa places on target were unmapped. Each marker's sites get an
  index of their own, so markers sharing a site never compete for the
  preset's hit limit. Hits are one locus when they share at least half of
  the shorter interval, both when a sequence's hits are merged and when a
  hit is matched to its marker's site.
- **Actions.** Paralog and ambiguous rows are dropped by default or flagged
  in the `PICRef` column, renamed `GenomeLabel`. Canonical rows are never
  dropped. Every input row, dropped or kept, gets a line in the labels
  sidecar.

Measurements (lettuce v8 genome, 2.4 Gb; 12-core laptop):

| | id95 MADCATS MADC | DArT vendor MADC |
|---|---|---|
| Rows / distinct sequences | 17,847 / 17,846 | 732,534 / 732,154 |
| Run time, prebuilt index | 12 s (8 threads) | 2 min 15 s (8 threads), 9.4 GB |
| Rows kept | 12,695 | 113,667 |
| Reads in paralog rows | 16.0% | 63.3% |
| Keep/drop agreement with bwa-mem labels | 98.2% | 94.9% |
| Keep/drop agreement with the bwa prototype's output, rows / reads | 98.2% / 99.7% | 94.9% / 99.4% |
| Het rate (proxy): unfiltered → prototype → `filter-paralogs` | 12.6% → 6.44% → 6.40% | 19.9% → 6.79% → 7.02% |

Indexing took 50 s and 10.6 GB at peak; the `.mmi` is 5.7 GB. The het
proxy sums Ref+RefMatch and Alt+AltMatch reads per marker and sample, calls
cells with ≥10 reads, and counts het where the minor allele has ≥10%. The
bwa-mem labels came from aligning the same haplotypes with `bwa mem -T 20`
(AS > XS unique). Variants tried and rejected:

- Default sr chaining: keep/drop agreement on id95 was 97.3%, against
  98.2% with single-seed chains and `best_n` 50.
- A k=15, w=10 index: DArT agreement 94.4%, against 94.9% with k=21, and
  slower mapping.
- Also dropping site-realigned DArT haplotypes with ≥6 mismatches:
  het proxy 6.93%, against 7.02%. The rest of the gap to bwa is minimap2
  finding fewer distant copies of highly diverged vendor haplotypes.

Count overflow. The lettuce vendor MADC stores counts in a signed byte:
90,913 cells across 1,108 markers are negative (minimum −128) and the
largest count is exactly 127. Counts of 128–255 appear as count − 256
(adding 256 matches MADCATS's counts for the same cells), and counts of 256
or more wrap silently. `read_madc` (and so `validate`) and `filter-paralogs`
now warn and record `count_overflow` in their reports. They never correct
the counts, since the wrapped ones cannot be told from real small counts.
The figures in the table above use the counts as stored.

HapApp layout. The id95 MADC rewritten in HapApp's layout (fixed IDs,
leading `Code_version`, every tenth row's counts as floats) gets the same
label, action and scores for all 17,847 rows as the raw MADC, and the same
QC totals. Moving the command into its own module left the id95 output, both
sidecars included, byte-identical.

Review fixes. An external review of the first version (cd2d7cb) gave a
counterexample for each of the defects below; each is now a test
(`ReviewCaseTests`, and `FilterTests` with the fake aligner). Effects on
the lettuce MADCs, default `use` mode:

| Fix | Effect |
|---|---|
| Ref tag realigned at its CloneID position on disagreement | `Chr2_001081815`: the genome-wide search returned 51 hits, best on Chr4 (96), none at its Chr2 design site (156). Placed on Chr2, 11 of its 14 dropped DArT rows become `on_target` |
| Every haplotype scored at its own site; the higher of the two scores kept | 530 DArT rows that the shared site index had crowded out get a site score; 1 row (`Chr5_323370531`, genome-wide 86 at the site, 88 elsewhere, 94 realigned) becomes `on_target` |
| One locus rule for merging hits and for target membership | No label changes; 850 DArT sequences had adjacent hits merged before |
| Every equal-best Ref-tag site kept for placement | No lettuce Ref tag has more than 8 |
| Multi-part `.mmi` refused | The lettuce index is one part; reading its layout takes 0.03 s |

Together, 12 DArT rows (994 reads) go from dropped to kept and none the
other way. Keep/drop agreement with bwa-mem on the 726,693 DArT Match/Other
rows with bwa labels rises from 94.8128% to 94.8136% (+6 rows); id95
labels do not change (97.3484% on 11,842 rows). The DArT run took 76 s,
against 78 s for cd2d7cb in the same session.
For `Chr2_001081815`, both the hit limit and the ambiguity codes matter: the
vendor tag, or the tag with R and M read as N, misses the design site at
`best_n` 50, 200 and 1,000 and finds it at 5,000; with the codes resolved to
bases it is found at `best_n` 50. bwa agreement is concordance, not truth:
source-known fixtures for both errors (target evidence removed, paralog
evidence kept) are still to be built.

## Genotype calling (`call`, unreleased)

`madcats call` (in `madcats_call.py`) calls genotypes from a MADC and writes
VCFs. The model is adapted from polyRAD 2.0.1: a beta-binomial likelihood of
each allele's read count against the rest of the record, a uniform or
Hardy-Weinberg prior, and (for `hwe`) iteration of the allele frequencies.
Only the parts needed for this were adapted. polyRAD's population-structure,
linkage and mapping-population priors, its ploidy tests and its import and
export code were not.

**Differences from polyRAD.**

- Likelihoods are kept as logarithms. polyRAD exponentiates them and sets a
  sample's likelihoods to 1 when they all underflow to 0.
- Each record iterates until its own allele frequencies settle. polyRAD stops
  all loci together on their mean change.
- There is an iteration limit (`--max-iterations`).
- Ties between dosages, and corrections with no possible step, give a missing
  genotype.

**Agreement with polyRAD.** `scripts/make_polyrad_fixture.R` runs polyRAD on
three simulated loci (2, 3 and 4 alleles; 16 samples of ploidy 2, 3, 4 and 6)
and stores the values in `tests/fixtures/polyrad_parity.json`:

| Compared | Agreement |
|---|---|
| Likelihoods, priors and posteriors from the same starting frequencies | 9 or more decimal places |
| Allele frequencies and posteriors after `IterateHWE(tol = 1e-10)` | 7 decimal places |
| Genotypes after the correction to the ploidy | identical (144 sample-by-allele values) |

**Correction to the ploidy.** With more than two alleles the most probable
dosage of each allele is found separately, and the dosages need not add up to
the ploidy. polyRAD adds or removes one copy at a time. An exact search for
the dosages with the highest product of posteriors
(`best_dosages_summing_to_ploidy`) gave the same genotype in all 4,104
corrected genotypes of a simulation (hexaploid, five alleles, 40 markers, 200
samples, 20–200 reads), where 51% of genotypes needed correcting. `call` uses
polyRAD's method and counts in the QC JSON how often the exact search differs.

**Overdispersion.** polyRAD's default of 9 allows more scatter in read
fractions than the gap between neighbouring polyploid dosages. On a simulated
mixed-ploidy table without extra scatter (500 markers, 300 samples, target
records, defaults otherwise):

| `--overdispersion` | 2x called / correct | 4x called / correct | 6x called / correct |
|---|---|---|---|
| 9 | 93.1% / 99.83% | 11.4% / 83.71% | 2.5% / 65.46% |
| 50 | 97.6% / 99.98% | 76.7% / 99.63% | 11.1% / 99.25% |
| 200 | 97.7% / 99.98% | 86.2% / 99.39% | 51.5% / 98.70% |
| `auto` (chose 300) | 97.7% / 99.98% | 86.9% / 99.33% | 55.2% / 98.34% |

`--overdispersion auto` is a maximum-likelihood estimate over a grid
(`estimate_overdispersion`). polyRAD's `TestOverdispersion` was ported first
and reproduced polyRAD's scores to four decimal places on three simulated
tables, but both chose 10 for a tetraploid table simulated with 50 and 4 for a
hexaploid table simulated with 20 (and 25 for a diploid table simulated with
20). It tests the read counts of samples that a rough estimate assigns one
allele copy, and in polyploids that estimate is wrong often enough to make the
counts look more scattered than they are. The port was dropped. The
maximum-likelihood estimate, on 200 simulated markers and 150 samples per
case, with the simulated value in the first column:

| Simulated | 2x | 4x | 6x | mixed 2x/4x/6x | 4x, 8–40 reads | 4x F1 | 6x F1 | 2x inbred |
|---|---|---|---|---|---|---|---|---|
| 5 | 5 | 5 | 5 | 5 | 5 | 6 | 7 | 2 |
| 10 | 10 | 10 | 10 | 10 | 10 | 12 | 14 | 2 |
| 20 | 20 | 20 | 20 | 18 | 20 | 25 | 30 | 2 |
| 50 | 50 | 50 | 50 | 50 | 50 | 50 | 75 | 2 |
| none (binomial) | 500 | 500 | 500 | 500 | 500 | 500 | 500 | 2 |

Inbred lines have no heterozygotes to measure scatter from and do not fit the
Hardy-Weinberg weights, so the estimate ends at the lowest value; `call` then
warns and uses 9. These are simulations. The estimate has not yet been checked
on real panels.

**Contamination as the only error allowance.** With the default contamination
rate of 0.001, a diploid with 10 Ref reads and 1 Alt read is called
heterozygous: a read of an allele the sample does not carry can only come from
contamination. This is polyRAD's behaviour.

**Other checks.** On the bundled example all 96 target and 96 haplotype
genotypes equal the simulated truth. `bcftools view` reads the three VCF types
without warnings, including mixed-ploidy output. Output is identical for any
`--jobs`. The simulated table above took 63 s on one core and 17 s on four for
all three VCF types.

## Downstream fixed allele IDs

This pipeline's responsibility ends with generation of the raw MADC-like
microhaplotype count table — and that boundary is certified end-to-end:
MADCATS output has been carried through the full production stack
(HapApp_utils fixed allele IDs -> `BIGr::madc2vcf_all()` -> updog
dosage) on two platforms and scored against external truth at the far
end. On GT-seq the whole downstream stack costs 0.25 points (98.42% at
the end of the chain vs 98.67% at the MADC), and the off-target SNPs
recovered from Match haplotypes score *higher* than the targets (98.79%
on 38,244 additional cells). On the 192-sample capture cohort the chain
delivers 95.13% on 295,068 target cells and reproduces the source
publication's population-genetic findings from MADCATS calls alone. Fixed allele-ID assignment and maintenance of the
cross-project panel allele database are downstream HapApp responsibilities;
the former `--allele-catalogue` reprocessing mode was removed accordingly.

## Input contracts

The DArT primer CSV must contain:
`MarkerName, LeftPrimerSequence, RightPrimerSequence, VariantRef, VariantAlt,
VariantPositionInTagSequence`. Indel markers use dash notation (`-` on the
empty side, the inserted or deleted bases on the other) with the position
naming where the affected span begins; the pipeline's stored variant
position becomes the first base at which the two tags actually differ,
which can sit beyond the declared position when the flank repeats the indel
content. Primers may contain IUPAC codes (degenerate oligos over flanking
SNPs); tags are reconstructed as the primer sequence verbatim plus the
genomic continuation, matching the vendor's representation.

The probe-design TSV/CSV must contain `MarkerName` and `TargetSequence`,
with a single bracketed variant — a SNP such as `[C/T]` or an
anchored-base indel such as `[A/ACTTG]` or `[TTTTTTG/T]`. Panel
construction rejects designs whose variant lies inside the left primer
(Alt reads could never carry the primer), identical alleles, and SNP
designs whose constructed Alt tag does not appear among the bracket
expansions (a primer-table/probe-design disagreement). An optional
`LocusGroup` column declares markers as tiles of one physical locus
(see the 0.9.0 section); members of a group must present the same
target region — identical Ref expansions, differing only in which
variant is bracketed — and violations are load-time errors. The
reference-tags CSV accepts and emits the same optional column, and it
only appears in output when some locus actually uses it, keeping
ungrouped panels byte-identical to earlier versions.

The platform-neutral reference-tag CSV must contain:
`MarkerName, LeftPrimerSequence, RefSequence, AltSequence`. Add
`VariantPositionInTagSequence` when Ref and Alt differ at more than one base;
it is validated as an in-range number pointing at an actual Ref/Alt
difference. Identical Ref/Alt sequences or alleles are rejected.
`VariantRef`, `VariantAlt` (validated for IUPAC compatibility with the tag
bases), and `RightPrimerSequence` are optional. Curated tags may contain
IUPAC codes; see the classification note above for what that means for
canonical counts.

Left primers that are strict prefixes of other panel primers are supported,
and under the default `--primer-mismatches 1` every primer within the
mismatch budget is considered regardless of depth or mismatch count — an
exact short primer cannot suppress a longer one-mismatch primer. Reads
matching several primers are resolved by full-tag identity, and
unresolvable ties are dropped as ambiguous. The default changed from 0 to
1 in 0.9.0: every validated benchmark on every platform ran with 1, and
exact-only matching silently drops Match-haplotype reads whose
distinguishing SNP falls inside the assignment anchor, which biases depth
by haplotype. Pass `--primer-mismatches 0` to restore exact-only
assignment.

Process a cohort in ONE `generate` invocation whenever novel-allele
discovery matters. The evidence model's reproducibility route counts
supporting samples within the invocation, so splitting a cohort into
several smaller `generate` runs weakens discovery in proportion to the
split: a real haplotype spread thinly across sample subsets can fail
the `--reproducible-min-samples` gate in every subset while easily
passing it cohort-wide (measured on a 192-sample capture cohort
processed as 8 batches: a materially thinner Match-row catalogue than
one full-cohort pass). Batch the expensive `normalize-reads` stage
freely — its outputs are per-sample and windows stay frozen — but feed
all samples to a single `generate`. The per-sample cache
(`--cache-dir`) makes re-running one large invocation cheap after the
first pass.

FASTQ input to `generate` is one primer-started single-end file per
sample; `normalize-reads` accepts any read that contains a panel window.
Paired-end data enters through one of the two 0.9.0 routes — `merge-pairs`
for overlapping pairs, or pooling both mates into one per-sample FASTQ for
non-overlapping dual-anchored designs — chosen explicitly from the assay's
measured geometry. Target manifests (`targets_*.csv`) are optional and used
only for sample names/metadata; by default sample columns use the unique
FASTQ filename stem.

## Tests

```bash
python3 -m unittest discover -s tests -v
```

The suite contains 59 unit/end-to-end tests, including evidence-route
protection and emission, single-sample depth discovery, the supporting-
frequency floor, exact and fuzzy prefix-nested primer resolution,
design-guard rejections, vendor-preamble sample-ID detection,
read-accounting invariants, infer-design (window enumeration, anchor
selection, exact and fuzzy end-to-end inference, ambiguous-window
exclusion), version 0.7.0's indel and degenerate-primer support (general
window enumeration, indel panel construction with round-trip, end-to-end
indel inference, primer overlay), version 0.8.0's normalize-reads
stage (window-placement invariants, orientation and extraction, and
end-to-end normalization into generation), and version 0.8.1's window
freezing and repeat-aware placement (verbatim `--reuse-tags` round-trip,
shifted-window honoring, design-mismatch refusal, and the tandem-repeat
penalty flipping a near-tied placement off a repeat), and version 0.8.2's
output framing (preamble structure with first-cell provenance and starred
cells, dotted statistics columns, sample-lookup renaming with
incomplete-table refusal, and `read_madc` round-trip of preambled
output), version 0.8.3's per-locus anchors (growth past the minimum
only for colliding prefixes, capped sharing for variant-identical
windows), and version 0.9.0's paired-end and tiling support (pair
merging across extended and engulfed/read-through geometries with
quality-picked overlap bases, locus-group tile routing end to end in
both `nearest` and `covered` modes including the winner-take-all
full-span case, and the distinct-window guard nudging a colliding
tile to its best distinct start and refusing an indistinguishable
homopolymer frame), and the `window-scan` diagnostic (a linked-SNP
fixture where a longer window earns strictly more information at equal
retention, and the scan recommends it), and version 0.10.0's gapped
rescue and review fixes (the bit-vector prescreen against a reference DP
on random indel cases, assignment and through-the-gap classification of
early, two-event and target-deleting indels, the composition guard
refusing a substitution-heavy read, the verified-suffix rule keeping the
substitution reading for a tail shift, closest-parent denoising, the
same-allele pooled pass on a private heterozygote, N-tag crediting and
dropping, intra-target k-mer deletion, R2 extension on unequal-length
merges, mate-ID errors, and the seeded scan subset), plus the second-review
fixes (N reads never credited to an ambiguity-coded canonical, design-side
IUPAC compatibility inside the gapped prescreen and alignment, dropped N
reads counted as window-scan churn, and length validation in the named
FASTQ parser). The diagnostic scripts in `scripts/` refuse
cache directories holding multiple fingerprints for one sample (pass
`--cache-fingerprint` to select one) and are runnable directly without
`PYTHONPATH` tricks. The 0.5.0 refactor was additionally gated on
byte-for-byte regressions: fresh per-sample processing reproduces every
0.4.0 cache payload exactly on all 80 samples, and the new quality-only
configuration reproduces the 0.3.x `full_s5` output file byte-for-byte.

## Interpretation and references

Breeding Insight describes MADC rows as 81-bp amplicon microhaplotypes with
per-sample counts. Ref and Alt contain only the designed target difference;
RefMatch and AltMatch have the respective target base plus other differences.
The official format overview is at
<https://breeding-insight.github.io/BIGr_documentation/posts/madc-file-format/>.

Breeding Insight's processing guidance discusses minimum per-sample support
for Other alleles and depth-based duplicate-site flags:
<https://breeding-insight.github.io/BIGr_documentation/posts/madc-data-processing/>.

The open `dartag_alleles` repository documents the reference-BLAST/paralog
audit pattern used for the historical diagnostics here:
<https://github.com/Breeding-Insight/dartag_alleles>.

Breeding Insight's public panel metadata and design Ref/Alt database are in
PanelHub:
<https://github.com/Breeding-Insight/BIGapp-PanelHub/tree/long_seq/sweetpotato>.

The sweetpotato panel publication gives the same classification definitions
and states that DArT supplies the MADC from its pipeline:
<https://pmc.ncbi.nlm.nih.gov/articles/PMC11354173/>.
