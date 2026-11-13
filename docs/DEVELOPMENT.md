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

Review of the module. Two readers went through `madcats_paralogs.py`, one for
structure and one for names, comments, sources and tests, as was done for
`madcats_call.py`. Three defects were fixed:

| Defect | Effect before | Now |
|---|---|---|
| A minimap2 index in several parts was refused only as `--genome-index` | Given as `--genome` it was loaded, mappy searched its first part, and markers on later contigs were judged without their own site | Refused whichever option names it |
| A processed MADC's `readCountSum` column was taken for a sample by `validate` and `filter-paralogs` (only `call` skipped it) | QC read totals of `filter-paralogs` doubled; `validate` listed it among the samples | Skipped by the helper all three use |
| A run that failed while writing (a non-numeric count in a late row) | A truncated MADC and labels file stayed on disk | Removed |

The rest changed no output: `filter_paralogs`, 211 lines doing five jobs, is
now five functions; both passes read their rows through one generator; the
MADC's header is checked before the index is loaded; names and comments were
corrected against the code (the genome-wide search returns up to 51 hits, not
50; a marker is placed by its CloneID position as well as its Ref tag; the
score is that of the reported alignment and not minimap2's `AS`); and the
module says what is minimap2's (Li 2018; Li 2021) and what is MADCATS's own.
The outputs of eight runs on the data below were the same before and after,
byte for byte, and the same on one thread and on four.

Rows with an empty `AlleleID`, or too short to hold the identity columns, are
passed over by both passes and are not counted in the QC totals. That is now
documented and tested; counting them would add a key to the QC JSON and was
left for a decision.

Delta smelt GT-seq: placement against the primers. The public smelt data
(see "Delta smelt GT-seq against alignment-based genotypes" under genotype
calling) give a second species and platform, and an independent answer for
where each marker lies: the 81 bases from its forward primer, found by exact
search of both primers in the genome (GCA_021917145.1, 437 Mb). The MADC has
368 markers, 1,527 rows and 10.2 million reads.

| Run | Markers at the expected site | Unplaced | Rows labelled |
|---|---|---|---|
| CloneIDs renamed to CHROM_POS of the target | 368 | 0 | 1,527 on_target |
| The same, `--clone-positions ignore` | 368 | 0 | 1,527 on_target |
| 29 markers' positions moved 5,000 bases and 15 put on another chromosome | 324 | 44, the moved ones | 1,350 on_target, 177 unplaced |
| The same, `--clone-positions ignore` | 368 | 0 | 1,527 on_target |
| CloneIDs as generated (locus names such as `NC_061060_1_U618`) | 354 | 14 | 1,494 on_target, 33 unplaced |

No row was labelled `paralog` or `ambiguous` and none was dropped: the panel
has no paralogs in this genome, so this checks placement and false positives.
A marker with a wrong position is left unplaced, not misplaced. In the last
run 100 locus names end in an underscore and digits and are read as
CHROM_POS; for 14 the number is not within the tag and the marker is left
unplaced. A CloneID that looks like a position is taken for one, which
`--clone-positions ignore` avoids. Eight runs took 38 s on four threads,
building the index once.

Planted paralogs. The record above says the agreement with bwa mem is
concordance and that cases with a known answer were still to be built. They
were built in the smelt genome: for each of 105 markers drawn at random, a
copy of the genome around its tag (the tag and 100 bases either side) was
changed and added to the genome as a contig of its own, and the MADC was
given the rows a co-amplified copy would give (the copy's tag, and that tag
with one more substitution, with half the Ref row's reads each). Fifteen
markers per class and three random draws; CloneIDs are CHROM_POS; default
options:

| Planted copy | Rows from the copy labelled `paralog` | Rows of the marker itself labelled `on_target` |
|---|---|---|
| None (263 markers per draw) | | 3,348 of 3,348 |
| 1 substitution from the Ref tag | 90 of 90 | 173 of 179 (4 `ambiguous`, 2 `paralog`) |
| 2 substitutions | 90 of 90 | 168 of 170 (2 `ambiguous`) |
| 4 substitutions | 90 of 90 | 165 of 165 |
| 8 substitutions | 90 of 90 | 183 of 183 |
| 1 substitution and a 2-base deletion | 90 of 90 | 167 of 167 |
| Identical copy | | 1 of 192; 191 `ambiguous`, of which 90 canonical rows were kept |
| 4 substitutions, copy not added to the genome | 0 of 90: all `on_target` | 177 of 177 |

- Every row from a copy that is in the genome was labelled `paralog` and
  dropped, 450 of 450, and no row of an unplanted marker was mislabelled.
- Of the 864 rows of markers with a copy in the genome, 8 were lost, all at
  copies one or two bases from the Ref tag: an allele of the marker that is
  as near the copy as its own site is `ambiguous`, and one that carries the
  copy's base is the copy's sequence.
- An identical copy makes the marker's rows `ambiguous`; the default drops
  them and keeps the canonical Ref and Alt rows.
- A copy that is not in the assembly is not seen: its rows are kept as
  alleles of the marker.

The same design, twenty markers on a random genome, is in the test suite
(`PlantedCopyTests`).

Alfalfa: two aligners and two families. The alfalfa 3K panel is a second real
panel with paralogs, in an autotetraploid, and its public data allow three
checks that the lettuce record lacks: where each marker lies, from more
sequence than its tag; each row's label from two other aligners; and where
each row maps genetically.

- MADCs: the two F1 families of the G3 article jkag131, with fixed allele
  IDs. I195 x J432 has 186 samples, 2,999 markers and 7,517 rows of 54-base
  tags; E115 x F628 has 374 samples, 2,989 markers and 8,731 rows of 81-base
  tags. The panel's haplotype database (23,469 haplotypes, Zenodo
  doi:10.5281/zenodo.10681524) was run as a third table.
- Genome: the monoploid assembly of cultivar Bolivia (GCA_051527215.1, eight
  gapless chromosomes, 776 Mb). The panel was designed on XinJiangDaYe, whose
  assembly is on figshare and was not reachable, so this is a genome of
  another cultivar: the panel's 361-base design sequences differ from it at a
  median of 3 bases, and at 17 or more for a tenth of the markers.
- A marker's home on this genome is where its design sequence (180 bases
  either side of the target, from the panel's probe table) aligns. bwa mem
  placed 2,745 of the 3,000 uniquely and minimap2 2,780; the 2,735 that both
  place uniquely are at the same base in both, and are the markers judged
  below. Another 22 have a unique place on a chromosome with another number
  than the panel names, and the rest have no unique place.
- Two runs of each table: with `--clone-positions ignore`, as a user must run
  it on a genome the CloneIDs do not refer to, and with the default after the
  CloneIDs were renamed to the chromosome and target position on this genome,
  as it runs on a design genome. The first run took 47 s and 5.2 GB with the
  index built, the others 7 to 8 s.

Placement, for the markers with a home:

| | E115 x F628, 81-base tags (2,730 markers) | I195 x J432, 54-base tags (2,734) |
|---|---|---|
| `ignore`: at home | 2,669 | 2,487 |
| `ignore`: several equal-best sites, one at home | 17 | 28 |
| `ignore`: somewhere else | 15 | 25 |
| `ignore`: not placed | 29 | 194 |
| Positions on this genome: at home | 2,730 | 2,727 (7 not placed) |

Without positions, a Ref tag that differs from this genome is placed where it
matches best, or nowhere: 0.5 to 0.9% of the markers are put in another copy,
and 1% (81 bases) to 7% (54 bases) are left unplaced with all their rows
kept. With positions, the tag is aligned again at its own site and every
marker but seven is placed.

Rows and reads by label, with positions: in E115 x F628 7,928 rows are
`on_target` (90.9% of reads), 357 `paralog` (4.3%), 297 `ambiguous` (3.8%),
9 `unmapped` and 140 `unplaced` (0.9%); 437 rows with 4.1% of the reads are
dropped, and 217 canonical rows with 4.0% are kept despite their label. In
I195 x J432 the shares are 88.7%, 4.7%, 4.5%, 0.3% and 1.8%, and 4.2% of the
reads are dropped. Of the database's 23,469 haplotypes 1,145 are `paralog`
and 930 `ambiguous`. Lettuce lost 63% of its reads; alfalfa has little to
remove.

Labels against bwa mem and blastn. Every row was aligned to the same genome
with bwa mem (`-a -k 12 -T 20`) and with blastn (word size 11, no masking),
and labelled from the aligner's own score by the rule the filter uses: the
best score is at home, is equalled elsewhere, or is elsewhere. Rows of the
markers placed at home; "same decision" joins the labels by what the default
does with the row (kept: `on_target`, `unmapped`; dropped: `paralog`,
`ambiguous`):

| Table, run | Rows | bwa mem: same label, same decision | blastn: same label, same decision |
|---|---|---|---|
| E115 x F628, positions | 7,916 | 99.42%, 99.81% | 99.42%, 99.68% |
| E115 x F628, `ignore` | 7,716 | 99.49%, 99.86% | 99.49%, 99.74% |
| I195 x J432, positions | 6,768 | 98.63%, 99.44% | 98.73%, 99.56% |
| I195 x J432, `ignore` | 6,145 | 99.27%, 99.56% | 99.35%, 99.63% |
| Haplotype database, positions | 21,155 | 99.14%, 99.53% | 98.90%, 99.28% |
| Haplotype database, `ignore` | 20,551 | 99.17%, 99.54% | 98.99%, 99.35% |

bwa mem and blastn give each other the same decision for 99.82% and 99.65% of
the two families' rows and 99.37% of the database's, so the filter is as near
to each of them as they are to one another. For RefMatch and AltMatch rows alone the filter and bwa mem
agree on 99.5% (81 bases) and 98.8% (54 bases). All 297 rows the filter calls
`paralog` in E115 x F628 are `paralog` for bwa mem.

Where the rows map genetically. In a full-sib family a haplotype is inherited
with the haplotypes near it on the same chromosome copy, so across the
progeny its share of the marker's reads rises and falls with theirs, and a
row from a paralogous copy follows the haplotypes around that copy instead.
For each row whose share varies in the family, the ten rows of other markers
that correlate with it most strongly (0.4 or more in absolute value, three or
more needed) were looked up at their markers' homes. The row is linked at
home when 60% or more of them lie within 15 Mb of its marker, and linked
elsewhere when 60% or more lie on one chromosome and fewer than 20% near the
marker. No alignment of the row enters this.

| Family, run | Rows judged | Kept rows linked at home | Rows labelled `paralog`: linked elsewhere, at home with the other site far away, at home with it near |
|---|---|---|---|
| E115 x F628, positions | 2,877 | 2,701 of 2,768 (97.6%) | 38, 7, 29 |
| E115 x F628, `ignore` | 2,816 | 2,671 of 2,736 (97.6%) | 36, 5, 21 |
| I195 x J432, positions | 2,585 | 2,425 of 2,485 (97.6%) | 22, 4, 31 |
| I195 x J432, `ignore` | 2,408 | 2,286 of 2,340 (97.7%) | 20, 1, 21 |

- Of the rows the filter keeps, 97.6% map to their marker in both families.
  More than half of the others are at some 30 markers where every kept row
  maps elsewhere: there the marker itself is somewhere else in the family
  than its design sequence is on this genome, which no row label can mend.
  The rest, 22 to 28 rows or about 1% of the kept rows judged, are rows from
  another place that the filter did not catch; a copy that this genome lacks
  or holds in another form is not seen.
- Of the `paralog` rows that linkage can decide, 84 to 95% map elsewhere, and
  for 25 of 38 and 18 of 22 the region they map to holds the site the filter
  named.
- For 21 to 31 `paralog` rows the better site is on the marker's own
  chromosome within 15 Mb, for most of them 1 kb to 1 Mb away. Such a copy is
  inherited with the marker, so linkage says nothing about them.
- RefMatch and AltMatch rows that map elsewhere are 9.2% of those judged
  before the filter and 5.5% after in the larger family, and 8.3% and 5.2% in
  the smaller.
- A RefMatch or AltMatch row with reads in every sample of the family is two
  to three times as common among `paralog` rows as among `on_target` rows:
  23% against 8%, and 43% against 13%.

Genotypes. `call --prior norm` was run on each family before and after the
filter and compared with the article's published dosages, which come from
the same reads. Over the 986 and 974 published targets the share of equal
dosages moved from 89.80% to 89.87% and from 90.11% to 90.33%. Of the 211
markers that lost a row in the larger family only 29 are among the
published ones (2.9%, against 7.1% of the panel), and there the share rose
from 86.3% to 86.9% with 11% of the genotypes changed. The estimated
overdispersion rose from 40 to 50 and from 30 to 40: with those rows gone the
counts scatter less.

A count of the progeny genotypes that the two parents' genotypes cannot give
was tried as a test and is not one. Reads of a paralog make a parent look
heterozygous, and a heterozygous parent allows every progeny genotype, so the
count rewards the rows the filter removes: it was 1.86% before and 2.20%
after at the markers that lost a row. At chr4.1_040870396 it went from 0 to
189 of 372 when one row was dropped. That row has a quarter to three quarters
of the marker's reads in every one of the progeny, which no allele of a
tetraploid cross has (two copies in each parent leave 1 progeny in 36 without
it); it differs from the Ref tag at five bases and matches a copy 10 kb away
better.

Not covered: the design genome; FASTQs of this panel, and so a MADC from
`generate` (no public alfalfa DArTag reads were found); planted copies in
this genome; and a haplotype-resolved reference, which has a copy of every
tag on each haplotype.

## Genotype calling (`call`, unreleased)

`madcats call` (in `madcats_call.py`) calls genotypes from a MADC and writes
VCFs. The model is adapted from polyRAD 2.0.1: a beta-binomial likelihood of
each allele's read count against the rest of the record, a uniform or
Hardy-Weinberg prior, and (for `hwe`) iteration of the allele frequencies.
Only the parts needed for this were adapted. polyRAD's population-structure,
linkage and mapping-population priors, its ploidy tests and its import and
export code were not.

Three parts of updog 2.1.7's model were added afterwards: the `norm` genotype
prior, the allelic bias, and a prior with a free weight for every dosage
(updog's `flex`), which `call` uses only while estimating the read model. Four
parts are MADCATS's own and have no published counterpart: the
maximum-likelihood estimate of the overdispersion, the table-wide
background-read rate estimated with it, the search that fits the bias, and
the `structure` prior. `madcats_call.py` names the polyRAD or updog function
each of its functions follows and says where it differs ("Audit against both
packages' sources" below). The paragraphs below are in the order the work was
done, so the earlier tables were measured before the later parts existed. The
figures for the finished model are under "Finished model on the alfalfa
populations", "Tables from updog's and polyRAD's simulators" and "Time and
memory".

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
records, other settings at their defaults):

| `--overdispersion` | 2x called / correct | 4x called / correct | 6x called / correct |
|---|---|---|---|
| 9 | 93.1% / 99.83% | 11.4% / 83.71% | 2.5% / 65.46% |
| 50 | 97.6% / 99.98% | 76.7% / 99.63% | 11.1% / 99.25% |
| 200 | 97.7% / 99.98% | 86.2% / 99.39% | 51.5% / 98.70% |
| `auto` (chose 300) | 97.7% / 99.98% | 86.9% / 99.33% | 55.2% / 98.34% |

`--overdispersion auto`, the default since the real-data comparison below, is
a maximum-likelihood estimate over a grid
(`estimate_overdispersion`). polyRAD's `TestOverdispersion` was ported first
and reproduced polyRAD's scores to four decimal places on three simulated
tables, but both chose 10 for a tetraploid table simulated with 50 and 4 for a
hexaploid table simulated with 20 (and 25 for a diploid table simulated with
20). It tests the read counts of samples that a rough estimate assigns one
allele copy, and in polyploids that estimate is wrong often enough to make the
counts look more scattered than they are. The port was dropped. The
maximum-likelihood estimate, on 200 simulated markers and 150 samples per
case (20–150 reads unless stated, no stray reads), with the simulated value in
the first column. The table is from the estimate as it now stands, with free
weights for the dosages and the background-read rate estimated alongside
(both described further down):

| Simulated | 2x | 4x | 6x | mixed 2x/4x/6x | 4x, 8–40 reads | 4x F1 | 6x F1 | 2x inbred |
|---|---|---|---|---|---|---|---|---|
| 5 | 5 | 6 | 8 | 6 | 6 | 5 | 7 | 2 |
| 10 | 10 | 10 | 10 | 10 | 10 | 10 | 10 | 2 |
| 20 | 20 | 20 | 20 | 20 | 20 | 20 | 20 | 2 |
| 50 | 50 | 50 | 50 | 50 | 50 | 50 | 50 | 2 |
| none (binomial) | 500 | 500 | 500 | 500 | 500 | 500 | 500 | 2 |

From 10 upward every value is recovered. The first version of the estimate,
which assumed Hardy-Weinberg proportions, read 12–30 for the F1 columns at
10 and 20, and 75 for the hexaploid F1 at 50. At a simulated 5, scatter wider
than anything seen in the real tables, polyploid estimates are 6 to 8,
because part of the scatter is taken for background reads. The
background-read rate came out at the lowest value tried, 0.0002, in 36 of the
40 cases, and at 0.001 to 0.01 in the other four: random-mating tetraploids
and hexaploids at 5, and the low-depth tetraploids at 5 and 10.

The inbred diploids here have no heterozygotes and no stray reads, so there is
no scatter to measure and the estimate ends at the lowest value; `call` then
warns and uses 9. On the two real alfalfa F1 populations below the first
version chose 40 and 50, and the finished one 30 and 40 with a background-read
rate of 0.003.

**Real alfalfa DArTag data (tetraploid, 3K panel).** Three processed MADCs
were called as tetraploids with the `call` model and with polyRAD 2.0.1
(`IterateHWE`, overdispersion 9, then `GetProbableGenotypes` with
`multiallelic = "correct"`) on the same counts: the two F1 populations of the
G3 article jkag131 (I195 × J432, 186 samples; E115 × F628, 374 samples) and
project DAl22-7249 (615 samples from 20 crosses and their parents).

| MADC | Target genotype values identical | Haplotype genotype values identical |
|---|---|---|
| I195 × J432 (tolerance 1e-7) | 1,103,352 of 1,103,352 | 1,385,886 of 1,385,886 |
| E115 × F628 (tolerance 1e-6) | 2,222,352 of 2,223,056 | 3,251,974 of 3,252,678 |
| DAl22-7249 (tolerance 1e-6) | 3,485,814 of 3,485,820 | 5,306,214 of 5,306,220 |

Every difference is at one marker (E115 × F628) or two (DAl22-7249) where
polyRAD had stopped while that marker's allele frequency was still moving: its
stopping rule is the mean change over all loci, and those frequencies were
0.01–0.03 from where they settle. Elsewhere the two tools' allele frequencies
differ by 1e-7 at the median.

The article publishes the dosages used for its linkage maps. Of its markers,
974 and 986 are panel targets, and the rest are off-target SNPs. Target calls
against the published dosages, over 166 and 338 progeny (the published tables
count Ref copies, and count the other allele at 38 and 43 markers, which were
re-oriented). These rows were measured with a contamination rate of 0.001 and
no background-read rate, which is how `call` ran at that point:

| Settings | I195 × J432: called / same dosage | E115 × F628: called / same dosage |
|---|---|---|
| Overdispersion 9, no filters (what polyRAD calls) | 100% / 88.81% | 100% / 86.33% |
| `--overdispersion auto`, no filters | 100% / 90.80% | 100% / 89.93% |
| Overdispersion 9, default filters | 40.9% / 93.17% | 41.3% / 91.78% |
| `--overdispersion auto`, default filters | 89.0% / 93.15% | 91.8% / 92.19% |
| `--overdispersion auto --prior uniform`, default filters | 86.6% / 93.01% | 90.2% / 91.98% |

In every row 98.5–99.4% of calls are within one copy of the published dosage.
The published dosages come from another caller, so a difference is not
necessarily an error of either. Calling from the Ref and Alt rows alone,
without the Match rows, lowered agreement to 80%, so the published calls also
counted Match rows.

polyRAD's other priors were run on the same target counts at the published
markers, with the overdispersion `call` estimated (40 and 50), to see what a
port would gain on these families:

| polyRAD prior | I195 × J432: called / same dosage | E115 × F628: called / same dosage |
|---|---|---|
| `IterateHWE` | 100% / 90.80% | 100% / 89.93% |
| `IteratePopStruct` | 100% / 90.74% | 100% / 89.91% |
| `PipelineMapping2Parents` (F1, no linkage) | 93.9% / 91.08% | 95.0% / 90.42% |

A single F1 family has no population structure, so the population-structure
prior changes nothing here. At these depths (medians of about 180 and 430
reads per sample and marker) the prior matters little.

**Against updog's norm model.** updog 2.1.7 (`flexdog`, `model = "norm"`) was
run on the same target counts at the published markers. Unfiltered agreement
with the published dosages was 92.98% and 92.68%, against 90.80% and 89.93%
for `call` with the HWE prior; updog took 235 s and 453 s on four cores. The
two callers give the same dosage for 96.99% and 96.54% of genotypes, and where
they differ the published dosage sides with updog about nine times in ten, so
the published dosages probably come from updog or a similar caller. Split by
the allelic bias updog fitted for each marker:

| Fitted bias of the marker | Share of genotypes | updog = `call` | updog = published | `call` = published |
|---|---|---|---|---|
| 0.8–1.25 | 73%, 69% | 99.19%, 99.28% | 93.59%, 92.91% | 93.09%, 92.50% |
| outside that, within 0.67–1.5 | 17%, 21% | 94.96%, 95.36% | 93.63%, 92.85% | 89.31%, 89.05% |
| outside 0.67–1.5 | 10%, 11% | 83.77%, 81.16% | 87.09%, 90.83% | 75.54%, 75.11% |

updog's per-marker fits had a median bias of 0.99 and 0.98 (10th–90th
percentile 0.72–1.24 and 0.72–1.30), a median sequencing-error rate of 0.09%
and 0.24%, and a median overdispersion of 395 and 183 on polyRAD's scale. The
single values `call` estimated for these tables (40 and 50) are lower because
they absorb the bias, which `call` did not model at that point.

On simulated tetraploid tables with known genotypes (100 markers, 150 samples,
20–60 reads, overdispersion 40, 0.5% read errors; one table per case),
accuracy of the most probable dosage:

| Population | updog norm | `call` defaults | `call`, contamination 0.01 and overdispersion 40 |
|---|---|---|---|
| Random mating | 87.5% | 84.8% | 89.2% |
| F1 | 88.7% | 88.2% | 89.4% |
| Two subpopulations | 88.8% | 80.3% | 89.5% |
| Inbred lines | 83.4% | 86.6% | 96.8% |
| Random mating, biased markers | 86.0% | 79.7% | 82.3% |
| F1, biased markers | 86.9% | 82.4% | 84.9% |
| Two subpopulations, biased markers | 89.0% | 77.9% | 85.3% |
| Inbred lines, biased markers | 84.9% | 86.2% | 96.6% |

Biased markers had a bias of 2^N(0, 0.5). Three things separate the columns.
(1) Read errors above the contamination rate pull the overdispersion estimate
down: with a true value of 40 it was 30, 16, 12 and 8 at error rates of 0.1%,
0.3%, 0.5% and 1%, and 40–50 throughout once the contamination rate was 0.01.
(2) The estimate assumed Hardy-Weinberg proportions and failed on the
two-subpopulation table (7) and on inbred lines (fell back to 9). (3) With
the right overdispersion and stray-read rate, `call` matches or beats updog on
unbiased markers and trails it by 2–4 points on biased ones. A prototype of a
normal-shaped prior in the `call` engine (not in the repository) changed these
figures by at most 1.5 points. The background-read rate, the free-weights
estimate and the allelic bias, all described below, answer these three points
in turn; the same eight tables are called again with the finished model under
"Tables from updog's and polyRAD's simulators".

**The norm prior (`--prior norm`).** updog's `norm` genotype distribution was
added as a prior: dosage probabilities proportional to a normal density, with
a mean and standard deviation fitted to the summed posteriors in every
iteration. The fit is exact, since the prior depends on the samples only
through those sums. It uses Newton's method on the two natural parameters,
where the function has a single peak, and searches the edges of updog's
parameter limits when the peak lies outside them; it takes 0.2 ms on average.
On 130 sets of dosage weights no point of a fine grid fitted better.

A Newton step that does not improve the fit is shortened by adding to the
diagonal of the curvature (Levenberg-Marquardt). The first version halved the
step instead, and stopped at its starting values when they put almost all of
the prior on one dosage: the curvature there is next to nothing and the step
is too long by tens of orders of magnitude. The audit found it in 1 of 40
single fits (weights 0.0003, 0.27 and 14.0 on the dosages of a diploid, start
at mean 1.98 and standard deviation 0.137), 0.39 below updog's optimum in mean
log prior. On 2,954 fits with many near-zero weights and arbitrary starts the
damped fit was never worse than the halved one and better in 337. Calls
changed little: the potato comparison below is the same to the last digit
(92.10% with `norm`, 96.57% with `norm` and the bias), and of twelve stored
runs the three that use `norm` changed in 8 of 257,400 target and haplotype
genotypes and in the off-target records of three markers, each with a higher
penalised likelihood than before.

`scripts/make_updog_fixture.R` stores updog 2.1.7 values in
`tests/fixtures/updog_norm.json`. For 14 sets of weights (ploidies 2 to 8;
one peak, rare and common alleles, two peaks, flat, a single class) the fitted
priors agree to 2e-5. For four simulated markers (random-mating tetraploid,
one hexaploid family, inbred diploid, rare-allele tetraploid), fitted by
updog with its bias, sequencing error and overdispersion held fixed, the
priors agree to 1e-5 and every genotype is the same. Posteriors agree to
1e-6 except for samples with a few reads of an allele, where updog's minimum
sequencing-error rate of 1e-6 leaves dosage 0 a probability of up to 5e-4.

Each ploidy among the samples is fitted separately, and a ploidy with fewer
than 10 samples uses Hardy-Weinberg proportions. On the two alfalfa F1
populations the norm prior agreed with the published dosages as often as HWE
did (90.38% against 90.80%, and 89.78% against 89.93%, unfiltered) and took
about 35% longer (29 s against 22 s, and 37 s against 27 s, for all 3,000
markers on four workers).

The processed layout of these files has a `readCountSum` column before the
samples, which `call` skips.

**Background-read rate (`--background-rate`).** In polyRAD a read of an allele
the sample does not carry can only come from contamination, at the allele's
frequency among the other samples. That cannot account for stray reads of an
allele that is rare in the table, and read errors beyond the contamination
rate were taken for scatter and pulled the overdispersion estimate down (point
1 above). `call` therefore has updog's sequencing-error term, under the name
background-read rate: a read shows the other allele with probability *e*, so
an expected fraction *f* becomes *f*(1 − *e*) + (1 − *f*)*e*. It is one value
for the table, where updog fits one per marker, and the default contamination
rate became 0. With polyRAD's settings a diploid with 10 Ref reads and 1 Alt
read is called heterozygous (posterior 0.95 in a table of 21 samples); with a
background rate of 0.002 and an overdispersion of 40 the same sample is
homozygous with a posterior of 0.69.

The rate is estimated together with the overdispersion
(`estimate_read_model`). For a given rate the overdispersion is the
maximum-likelihood value on its grid, and the rate kept is the one at which
the counts are most probable, found by starting at 0.002 and moving to a
neighbouring rate (eight values from 0.0002 to 0.02) for as long as that
raises the likelihood. The two can be told apart because they show in
different samples: stray reads in samples that do not carry an allele, scatter
in samples that carry both.

Eighty simulated tables of 50 markers and 120 samples were estimated: random
mating, one F1 family, two subpopulations and inbred lines, each with a
background rate of 0.001, 0.003 and 0.01 at an overdispersion of 40, and with
0.003 at 150. Estimates of the overdispersion, and the number of tables whose
background rate was the simulated one:

| Samples | Reads per sample | Simulated 40 (12 tables) | Simulated 150 (4 tables) | Background rate |
|---|---|---|---|---|
| 4x | 20–120 | 40 in all 12 | 150 in all 4 | 16 of 16 |
| 2x | 20–120 | 40 in all 12 | 150 in 3, 200 in 1 | 16 of 16 |
| 6x | 30–150 | 40 in 8; 14–50 | 150 in all 4 | 15 of 16 |
| 4x | 6–25 | 40 in 6; 30–75 | 150 in 2; 100–500 | 15 of 16 |
| 2x | 6–25 | 40 in 10; 30–50 | 100–300 | 16 of 16 |

The background rate is the simulated one in 78 of the 80 tables and one grid
step away in the other two. The overdispersion is the simulated one in 31 of
the 32 diploid and tetraploid tables with 20 or more reads. The two
hexaploid estimates furthest out (14 and 25) are at a background rate of 1%.
With 6 to 25 reads an overdispersion of 150 cannot be measured: at that depth
scatter that small is hidden by the sampling of reads, and for the same
reason it changes the calls little. Inbred lines are estimable here because
their homozygotes carry stray reads.

On the alfalfa populations the estimate is a background rate of 0.003 with an
overdispersion of 30 (I195 × J432) and 40 (E115 × F628). With the bias fitted
it is 0.0005 with 150, and 0.001 with 150, nearer updog's per-marker medians
for the same tables (0.09% and 0.24%; 395 and 183).

**Free weights for the estimate.** The estimate sums the genotypes out, so it
needs genotype proportions for each marker. With Hardy-Weinberg proportions it
failed on two subpopulations and on inbred lines (point 2 above). With the
`norm` shape, which cannot describe an excess of both homozygotes, it still
failed on inbred and partly inbred tables. The estimate now gives every dosage
its own probability at each marker (`PRIOR_FREE`: the summed posteriors,
normalised), whatever `--prior` the run uses. This is the form of updog's
`flex` model. updog's default is `norm`, and it suggests `norm` in place of
`flex` when there are fewer than 10 × (ploidy + 1) samples, since `flex` has
one free weight per dosage at every marker. Here it is used only for the
estimate, which pools the markers, and never for the calls.
Estimates without bias fitting,
as overdispersion and background rate, on eight of the simulated tables
described under "Tables from updog's and polyRAD's simulators":

| Table | Simulated | `norm` shape | Free weights |
|---|---|---|---|
| polyRAD 4x, inbreeding 0.5 | 20; 0.001 contamination and 0.001 errors | 2 (fails), 0.0002 | 20, 0.001 |
| polyRAD 4x outbred | same | 20, 0.001 | 20, 0.001 |
| polyRAD 6x outbred | same | 20, 0.001 | 20, 0.001 |
| updog 4x random mating, biased | 49–499, 0.005 | 25, 0.01 | 25, 0.01 |
| updog 4x F1, biased | 49–499, 0.005 | 30, 0.01 | 30, 0.01 |
| Two subpopulations | 40, 0.005 | 40, 0.005 | 40, 0.005 |
| Inbred lines | 40, 0.005 | 2 (fails), 0.0002 | 40, 0.005 |
| F1, biased | 40, 0.005 | 30, 0.005 | 30, 0.005 |

Free weights give the same values as the `norm` shape where that worked, and
the simulated values on both tables where it failed. On biased tables both
read low unless the bias is fitted, when both gave 100 and 0.005 on the two
updog tables. The inbred table is estimable because its homozygotes carry
stray reads; a table with neither heterozygotes nor stray reads still ends at
the lowest overdispersion tried, and `call` then warns and uses 9.

**Allelic bias (`--allele-bias fit`).** The comparison with updog put the gap
at markers where a copy of one allele yields more reads than a copy of the
other (point 3 above). updog's bias was added for records with two alleles:
with bias *h*, an expected reference fraction *f* becomes
*f* / (*h*(1 − *f*) + *f*). As in updog, the log of the bias has a normal
prior (mean 0, standard deviation 0.7), the search starts from five values
(exp(−1) to exp(1)), and the fit with the highest log-likelihood after that
prior is kept. Two things differ from updog:

- The bias is applied to the sample's own allele copies, and background reads
  are added afterwards. updog applies its sequencing error first. In updog's
  order, samples that do not carry an allele inform its bias through their
  stray reads; with one background rate for the whole table that would turn
  markers with more or fewer stray reads than average into biased ones.
- The search is over the bias alone. Each trial value is a complete fit of the
  record at that bias, scored by its log-likelihood with the genotypes summed
  out. The five starting values are fitted from scratch, and every one that
  scores at least as high as its neighbours is refined by a golden-section
  search within half a log unit of it, to 0.01 in the log of the bias. updog
  optimises bias, sequencing error and overdispersion together from each
  start. The range searched is exp(−2) to exp(2), 0.14 to 7.4.

`tests/fixtures/updog_norm.json` holds ten updog fits with the bias estimated
and the other parameters fixed: five simulated markers (tetraploids simulated
with a bias of 0.6, 1.9 and 1; a diploid with 0.5; a hexaploid family with
1.5), each under updog's `norm` and `hw` models. The fitted biases agree to
0.07%, the priors to 1e-4 and the posteriors to 1e-3, and all 1,140 genotypes
are the same. For the hexaploid family both tools fit a bias of 2.1 (`norm`)
and 2.7 (`hw`) where 1.5 was simulated: with few dosage classes a bias and a
shift of every dosage explain the counts about equally well.

The first version of the search refined only the best starting value. On
updog-simulated F1 markers with one or two dosage classes it then kept the
wrong one (90.7% of genotypes correct against updog's 94.0%); refining every
start that scores at least as high as its neighbours gave 93.3%.

Haplotype records and sites with three or more bases are called without a
bias: updog's bias is defined for two alleles, and a bias per allele of a
multi-allele record has not been worked out or tested.

**Finished model on the alfalfa populations.** Target calls against the
published dosages, without filters, with the estimates above:

| Settings | I195 × J432: same dosage | E115 × F628: same dosage |
|---|---|---|
| polyRAD's settings (overdispersion 9, contamination 0.001) | 88.81% | 86.33% |
| `--prior hwe` (the defaults) | 90.87% | 90.04% |
| `--prior norm` | 90.11% | 89.80% |
| `--prior hwe --allele-bias fit` | 92.17% | 91.79% |
| `--prior norm --allele-bias fit` | 92.26% | 91.94% |
| updog 2.1.7, `norm` | 92.98% | 92.68% |

In the four rows of the finished model 98.3–98.7% of calls are within one copy
of the published dosage. The published dosages probably come from updog
(above), so the last row is partly a comparison of updog with itself. With
`--prior norm`, of the 2,966 and 2,972 target records fitted, 38% and 40%
have a bias outside 0.8–1.25 and 26% and 27% outside 0.67–1.5 (median 1.0),
and 3% end at a limit of the search.

**Potato DArTag against SNP-array genotypes.** The alfalfa dosages above
probably come from updog, so they cannot say whether updog or `call` is nearer
the truth. The potato data of Endelman et al. (2024, The Plant Genome,
doi:10.1002/tpg2.20484) have genotypes from a second technology. The
polyBreedR repository (`inst/vignette_data`) carries the part of that dataset
from a five-parent half-diallel: 85 tetraploid clones in ten F1 families, with
DArTag read counts (`diallel_DArTag.vcf.gz`, the `AD` field) and Infinium
array genotypes called with fitPoly (`diallel_array.vcf.gz`). The two files
share 1,544 markers with the same REF and ALT; the median depth is 304 reads
per sample and marker. The VCF has no tag sequences, so `call` was given the
reference and alternate counts as a two-row-per-marker MADC with stand-in
sequences. Most probable dosage against the array dosage, without filters,
over the 131,224 genotypes with reads:

| Caller | Same dosage | RMSE | Markers with under 10% wrong | Bias 0.8–1.25 (1,135 markers) | Bias outside (409) |
|---|---|---|---|---|---|
| `call` with polyRAD's settings | 87.92% | 0.353 | 905 | 91.9% | 76.9% |
| polyRAD `IterateHWE`, overdispersion 9 | 87.92% | 0.353 | 905 | 91.9% | 76.9% |
| polyRAD `IteratePopStruct`, overdispersion 9 | 90.53% | 0.314 | 1,095 | 94.5% | 79.5% |
| `call --prior hwe` (the defaults) | 92.26% | 0.290 | 1,199 | 96.7% | 79.9% |
| `call --prior norm` | 92.10% | 0.294 | 1,191 | 96.7% | 79.4% |
| `call --prior hwe --allele-bias fit` | 96.60% | 0.192 | 1,396 | 97.1% | 95.2% |
| `call --prior norm --allele-bias fit` | 96.57% | 0.192 | 1,394 | 97.1% | 95.2% |
| updog `norm` on the same 85 samples | 96.69% | 0.188 | 1,403 | 97.2% | 95.4% |
| updog `norm` as stored in the VCF (fitted on 703 samples) | 96.86% | 0.181 | 1,399 | 97.3% | 95.7% |

The last two columns split the markers by the bias updog stored in the VCF.
Every row has 99.8–99.9% of calls within one copy.

- `call` with polyRAD's settings and polyRAD give the same share to two
  decimals, as on the alfalfa data.
- The estimated read model (overdispersion 75, background rate 0.003) is
  worth 4.3 points over polyRAD's defaults.
- Fitting the bias is worth another 4.3 points, 15 at the 409 markers updog
  found biased. The estimate with the bias fitted is 300 and 0.001, and 31% of
  the fitted biases are outside 0.8–1.25.
- With the bias fitted, `call` is 0.1 points below updog on the same samples
  and took 52–56 s on four workers against updog's 150 s (CPU 197–214 s
  against 394 s; memory 150 MB against 693 MB). Without the bias it took 6–7 s.
- polyRAD's population-structure prior adds 2.6 points at an overdispersion
  of 9, where the prior carries more weight, and still ends below `call` with
  the HWE prior and estimated values. With an overdispersion of 75 polyRAD's
  HWE prior gives 92.41%, leaving 210 genotypes uncalled, and its structure
  prior stops with an error (a regression with no usable cases).

**The whole potato deposit, and the paper's own error table.** The deposit
of Endelman et al. (2024) is on Dryad (doi:10.5061/dryad.8pk0p2nw4, CC0) with
a copy on Zenodo (record 11286799). File S1 has the DArTag read counts of 703
tetraploid samples at the 2,503 markers of the first panel version, File S5
the array genotypes of 298 of them, and File S12 the paper's error of the
DArTag genotypes (updog's `norm`, stored in File S1) against the array for
each of the 1,862 markers that are in both files with the same alleles. All
703 samples were called together, without filters, and scored on the 298 with
array genotypes. The median depth is 201 reads per sample and marker.

The first three columns are File S12's measures, per marker over the genotypes
with an array dosage, as mean and median over the 1,862 markers: the root mean
squared error, the share with another dosage, and the share in another class
of absent, mixed and fixed. The last column is the share of all 550,052
genotypes with reads that have the array's dosage.

| Caller | RMSE: mean, median | Another dosage: mean, median | Another class: mean, median | Same dosage, all genotypes |
|---|---|---|---|---|
| File S12 as published | 0.363, 0.224 | 13.4%, 5.0% | 6.5%, 0.7% | |
| updog `norm` as stored in File S1 | 0.363, 0.224 | 13.4%, 5.0% | 6.5%, 0.7% | 87.15% |
| `call --prior hwe --allele-bias fit` | 0.347, 0.224 | 12.9%, 5.0% | 5.8%, 0.7% | 87.46% |
| `call --prior norm --allele-bias fit` | 0.349, 0.232 | 12.9%, 5.0% | 5.9%, 0.7% | 87.45% |
| `call --prior hwe` (the defaults) | 0.415, 0.278 | 17.5%, 7.7% | 7.1%, 0.7% | 82.81% |
| `call --prior norm` | 0.418, 0.278 | 17.6%, 7.7% | 7.2%, 0.7% | 82.74% |
| `call --prior structure` | 0.412, 0.278 | 17.3%, 7.4% | 7.1%, 0.7% | 83.02% |
| polyRAD `IterateHWE`, overdispersion 75 | 0.414, 0.284 | 17.4%, 7.7% | 6.7%, 1.0% | 83.02% |
| polyRAD `IteratePopStruct`, overdispersion 75 | 0.413, 0.278 | 17.2%, 7.7% | 6.7%, 1.0% | 83.23% |
| `call` with polyRAD's settings | 0.498, 0.393 | 22.6%, 15.4% | 7.9%, 4.0% | 77.73% |
| polyRAD `IterateHWE`, overdispersion 9 | 0.502, 0.393 | 22.6%, 15.4% | 7.9%, 4.0% | 77.72% |
| polyRAD `IteratePopStruct`, overdispersion 9 | 0.474, 0.362 | 20.4%, 12.8% | 7.4%, 3.4% | 80.03% |

- The stored genotypes scored here give File S12's figures to the last
  decimal, so these columns measure what the paper measured.
- With the bias fitted `call` has the paper's median error and a lower mean:
  12.9% of genotypes with another dosage against 13.4%. The estimate was an
  overdispersion of 150 and a background rate of 0.003.
- Without the bias (estimate 75 and 0.005) `call` is 4.6 points behind, level
  with polyRAD given the same overdispersion and 5 points ahead of polyRAD's
  defaults. Fitting the bias is worth as much here as on the 85 clones.
- The structure priors add 0.2 points in both tools. polyRAD's finished at an
  overdispersion of 75 with 703 samples; with 85 it stopped with an error.
- Agreement is lower than on the 85 clones for every caller because the
  markers are not selected: at 331 of the 1,862 markers fewer than 80% of the
  stored genotypes have the array's dosage, and at 153 fewer than half. No
  sample is under 60%, and the median marker is at 95%. The disagreement is
  between the two assays at particular markers, which is why the medians say
  more than the means.
- On four workers that other jobs were also using, `call` took 19 s with the
  HWE prior, 38 s with the structure prior, and 233 s (HWE) and 432 s (norm)
  with the bias fitted; polyRAD's structure model took 90 s.

**A vendor MADC with diploid and tetraploid samples (potato, second panel
version).** File S4 of the same deposit is the MADC as DArT delivered it for
the second panel version: 375 samples, 52 of them diploid, and 3,915 markers.
File S3 has the paper's genotypes of those samples, and File S13 its error
against the array per marker, as the share of genotypes in another class of
absent, mixed and fixed (the measure that suits both ploidies). `call` read
File S4 as it is, with `--ploidy 4` and a `--ploidy-file` naming the diploids
(each sample's ploidy was read from its genotypes in File S3). It wrote 3,906
target records, 5,879 records in the `all` VCF and 3,915 haplotype records;
nine markers have no target record because their Ref and Alt tags differ at
more than one position. The median depth is 273 reads. 78 samples are on the
array (38 diploid, 40 tetraploid), and 2,603 markers are in both files with
the same alleles and have a target record: 202,143 genotypes with reads.

| Caller | Same dosage: all | Diploids (98,339) | Tetraploids (103,804) | Another class, per marker: mean, median |
|---|---|---|---|---|
| File S13 as published | | | | 5.3%, 3% |
| Genotypes stored in File S3 | 91.83% | 92.61% | 91.09% | 5.3%, 2.6% |
| `call --prior hwe --allele-bias fit` | 91.91% | 92.61% | 91.25% | 5.1%, 2.6% |
| `call --prior norm --allele-bias fit` | 91.89% | 92.58% | 91.23% | 5.1%, 2.6% |
| `call --prior hwe` (the defaults, no filters) | 89.14% | 92.02% | 86.41% | 5.7%, 2.6% |
| `call --prior norm` | 88.80% | 92.04% | 85.74% | 5.7%, 2.6% |
| `call --prior structure` | 89.31% | 92.03% | 86.74% | 5.7%, 2.6% |
| polyRAD `IterateHWE`, overdispersion 30 | 89.38% | 92.28% | 86.62% | 5.4%, 2.6% |
| polyRAD `IteratePopStruct`, overdispersion 30 | 89.60% | 92.31% | 87.03% | 5.4%, 2.6% |
| `call` with polyRAD's settings | 86.25% | 92.62% | 80.21% | 6.6%, 3.8% |
| polyRAD `IterateHWE`, overdispersion 9 | 86.23% | 92.60% | 80.20% | 6.6%, 3.8% |
| polyRAD `IteratePopStruct`, overdispersion 9 | 87.80% | 92.64% | 83.22% | 6.2%, 3.8% |

File S13 is rounded to two decimals; the stored genotypes scored here give a
mean of 5.30% against its 5.28%.

- A vendor MADC of mixed ploidy is called in one run, and with the bias
  fitted the result is level with the paper's genotypes for the diploids and
  0.2 points above them for the tetraploids. The estimate was an
  overdispersion of 100 and a background rate of 0.003; 1,507 of the 3,904
  fitted biases are outside 0.8 to 1.25.
- Without the bias (estimate 30 and 0.005) the tetraploids are 4.8 points
  behind, level with polyRAD at the same overdispersion and 6.2 points ahead
  of polyRAD's defaults.
- For the diploids the defaults are 0.6 points behind polyRAD's settings. The
  two differ in 1,512 diploid genotypes at 528 markers. In 1,034 the array has
  a heterozygote that the defaults call as fixed and polyRAD's settings as
  mixed; in 445 the array has a fixed genotype that the defaults get right.
  These are markers where one allele gives few reads: a diploid heterozygote
  with 210 reference and 47 alternate reads, at a marker with a fitted bias of
  0.27, is called without the alternate allele under an overdispersion of 30,
  because with a background rate the homozygote's distribution has a long
  tail, and as mixed under an overdispersion of 9, which lets a heterozygote
  stray further from half. Fitting the bias settles it: 92.61%.
- With the default filters 87.5% of the genotypes are called (98.4% of the
  diploid, 77.1% of the tetraploid) and 92.2% of those have the array's
  dosage.
- On shared workers the run that wrote all three VCFs took 88 s, the norm
  prior 34 s, the structure prior 63 s, and the runs with the bias 159 to
  169 s; polyRAD's structure model took 113 to 200 s.

**Which rows a target genotype counts.** `call` counts every row of a marker
by the base it has at the target, so a RefMatch or AltMatch haplotype adds
its reads to the allele it carries. A workflow that takes the exact Ref and
Alt rows only, as the sweetpotato study below did, leaves those reads out.
The potato vendor MADC and the array decide between the two:
File S4 was cut down to its 7,830 Ref and Alt rows (of 13,187) and called
again.

| Counts | `--prior hwe`: all, diploids, tetraploids | With the bias fitted |
|---|---|---|
| Every row, by its base at the target | 89.14%, 92.02%, 86.41% | 91.91%, 92.61%, 91.25% |
| Ref and Alt rows only | 82.06%, 85.25%, 79.04% | 84.36%, 85.77%, 83.03% |

With the Ref and Alt rows alone 7 points fewer genotypes have the array's
dosage, and 2.1% of the genotypes have no reads at all. The reads on the
other rows are reads of the target's alleles.

**Hexaploid sweetpotato: replicates of eight clones.** The deposit of a study
of DNA from sweetpotato leaves and roots (Zenodo record 22288064, CC-BY) has
a vendor MADC of 80 samples, ten of each of eight clones, and the study's
genotypes from updog's `norm` model (`multidog`, ploidy 6, written by BIGr's
`updog2vcf`). There is no second assay, so the measure is how often two
samples of one clone get the same dosage. The MADC as delivered names each
sample twice, once per extraction method, and `call` refuses repeated sample
names; the second of each pair was given a suffix.

The study's genotypes use the exact Ref and Alt rows, so `call` was first run
on those rows alone (6,240 of 19,319), which gives both tools the same
counts. 2,636 markers are in every output and 622,513 pairs of samples of one
clone both have reads:

| Caller | Pairs with the same dosage | Different dosages among the eight clones, per marker | Same genotype as the study's |
|---|---|---|---|
| updog `norm`, as stored in the deposit | 87.31% | 2.16 | |
| `call --prior norm --allele-bias fit` | 86.68% | 2.18 | 90.86% |
| `call --prior norm` | 86.14% | 2.12 | 80.76% |
| `call --prior hwe --allele-bias fit` | 85.66% | 2.19 | 87.48% |
| `call --prior hwe` | 84.79% | 2.17 | 80.24% |
| `call` with polyRAD's settings | 80.02% | 2.00 | 74.57% |
| `call --prior hwe`, default filters | 99.02% of the 53.0% of pairs with both called | 1.77 | 93.70% |
| `call --prior norm --allele-bias fit`, default filters | 96.90% of the 65.9% | 2.01 | 93.69% |

- On the same counts `call` with the norm prior and the bias fitted is 0.6
  points behind updog in agreement between replicates and gives the study's
  genotype in 91% of cases. The eight clones are ten samples each, which is
  not a population in Hardy-Weinberg proportions; the norm prior suits it
  better than the HWE prior by 1 to 1.4 points.
- The estimate was an overdispersion of 50 and a background rate of 0.002,
  and 150 and 0.001 with the bias fitted. With the default filters two
  thirds or less of the pairs are called, and 97 to 99% of those agree.
- With every row counted, as `call` does by default, the picture is another
  one: 72.8% (HWE) to 75.9% (norm) of pairs agree, the clones show 3.2
  different dosages per marker instead of 2.2, and only 36% of the genotypes
  are the study's. At Chr01_000084128 the Alt row has 7 reads in all 80
  samples and the study's genotypes have no alternate allele; AltMatch rows
  with the alternate base carry 7 to 29% of each sample's reads, alike within
  a clone and different between clones, and `call` gives dosages of 1 and 2.
- Agreement between replicates cannot rank the two count sets, because a
  marker that looks monomorphic agrees with itself. The potato comparison
  above can: there the rows left out are reads of the target's alleles. The
  lower agreement with every row counted is what more heterozygous
  genotypes cost in a hexaploid at this overdispersion, and no truth for
  these clones is at hand to say how many of them are right.

**Tables from updog's and polyRAD's simulators.** To have truth that does not
come from MADCATS code, read counts were also drawn with the simulators the
two packages ship (`rgeno` and `rflexdog` in updog; `SimGenotypes` and
`SimAlleleDepth` in polyRAD), 150 markers by 150 samples per table, with
negative-binomial depths of mean 150, or 25 where the table is marked low
depth. updog's tables have a sequencing-error rate of 0.5%, an overdispersion
of 49–499 on polyRAD's scale drawn per marker and, where marked biased, a bias
of exp(N(0, 0.5)) per marker. polyRAD's have an overdispersion of 20, a
contamination rate of 0.001 and an error rate of 0.001. The first eight rows
are the tetraploid tables of the earlier comparison (100 markers, 150 samples,
20–60 reads, overdispersion 40, 0.5% read errors, bias 2^N(0, 0.5) where
marked). Share of genotypes whose most probable dosage is the true one,
without filters; polyRAD is `IterateHWE` with its defaults:

| Table | polyRAD | updog `norm` | `call` hwe | `call` norm | hwe + bias | norm + bias |
|---|---|---|---|---|---|---|
| 4x random mating | – | 87.5% | 89.5% | 89.5% | 89.1% | 88.9% |
| 4x F1 | – | 88.7% | 89.5% | 90.7% | 85.6% | 89.9% |
| 4x two subpopulations | – | 88.8% | 89.8% | 90.2% | 87.2% | 88.7% |
| 4x inbred lines | – | 83.4% | 97.4% | 98.6% | 71.7% | 90.3% |
| 4x random mating, biased | – | 86.0% | 81.8% | 81.3% | 87.4% | 87.4% |
| 4x F1, biased | – | 86.9% | 85.5% | 85.2% | 85.9% | 88.9% |
| 4x two subpopulations, biased | – | 89.0% | 85.4% | 85.4% | 88.7% | 89.6% |
| 4x inbred lines, biased | – | 84.9% | 97.3% | 98.4% | 72.8% | 89.8% |
| updog 4x random mating | 84.3% | 96.7% | 96.8% | 96.7% | 96.7% | 96.7% |
| updog 4x random mating, biased | 72.3% | 96.3% | 82.3% | 81.5% | 96.2% | 96.2% |
| updog 4x F1, biased | 76.1% | 94.0% | 81.6% | 81.1% | 87.6% | 93.3% |
| updog 4x `norm` genotypes, biased | 74.6% | 95.9% | 84.2% | 83.6% | 95.8% | 95.9% |
| updog 2x random mating, biased | 99.5% | 100.0% | 99.6% | 99.6% | 100.0% | 100.0% |
| updog 6x random mating, biased | 58.3% | 86.4% | 67.3% | 66.6% | 86.5% | 86.3% |
| updog 4x random mating, biased, low depth | 74.5% | 84.1% | 76.3% | 75.6% | 84.6% | 84.4% |
| polyRAD 4x outbred | 87.0% | 87.3% | 89.0% | 88.8% | 88.5% | 88.5% |
| polyRAD 4x, inbreeding 0.5 | 88.7% | 93.6% | 93.2% | 93.9% | 90.7% | 93.1% |
| polyRAD 2x outbred | 99.7% | 99.9% | 99.9% | 99.9% | 99.9% | 99.9% |
| polyRAD 6x outbred | 73.1% | 65.7% | 75.3% | 75.2% | 71.7% | 65.3% |
| polyRAD 4x outbred, low depth | 81.5% | 79.2% | 82.4% | 82.4% | 81.8% | 81.7% |

`call` leaves a sample without reads uncalled and the other two call it from
the prior; such samples are 0.1% of the low-depth tables and are not counted
for `call`. What the table shows:

- Without bias in the data, `call` without bias fitting is within half a point
  of updog or ahead of it on every table: ahead by 0 to 3 points on outbred
  tetraploids, 10 on the polyRAD hexaploid and 14 to 15 on inbred lines. On
  polyRAD's partly inbred table `hwe` is 0.4 behind and `norm` 0.3 ahead.
- With bias in the data, `--prior norm --allele-bias fit` is within 0.7 points
  of updog on all seven of updog's tables and ahead of it on the four biased
  tables of the first group. Without bias fitting `call` is 8 to 20 points
  behind updog on updog's biased polyploid tables.
- Fitting a bias that is not there costs 0 to 1.5 points on outbred
  tetraploids with `--prior norm`. Fitting one costs 8 points on inbred lines,
  whether or not the data are biased, and 10 on the polyRAD hexaploid. updog,
  which always fits a bias, is at 83–85% and 66% on those tables. Without
  heterozygotes, or with seven dosage classes of scattered counts, a bias and
  a shift of the dosages are hard to tell apart.
- `--prior hwe --allele-bias fit` fails where the samples are not in
  Hardy-Weinberg proportions: 72–73% on inbred lines, and 86–88% on the three
  F1 tables against 89–93% with `norm`. The bias needs a prior that can follow
  the genotype proportions.
- polyRAD with its default overdispersion of 9 is 2 to 12.5 points below
  `call` hwe, the same model with estimated values, on updog's polyploid
  tables, whose scatter is far smaller than 9 allows. On polyRAD's own tables,
  simulated with 20, it is 0 to 4.5 points below.
- With no bias fitted, the mean probability `call` gives its called dosage is
  within 1.5 points of its accuracy on every table without bias in the data.
  With `--prior norm --allele-bias fit` it is within 3 points on the outbred
  tetraploid tables, apart from updog's F1 table (4 above). updog's is 4.5 to
  8 points above its accuracy on the outbred tetraploid tables that are not
  its own and 0 to 4 above on its own, 21 above on the polyRAD hexaploid and
  22 below on inbred lines.
- With the bias fitted, the overdispersion estimate was one grid step high on
  several unbiased tables (50 for a simulated 40, 25 for 20), and wrong on the
  polyRAD hexaploid: 10 with a background rate of 0.0002, where 20 and 0.001
  were simulated and the estimate without bias gave 20 and 0.001.

**Time and memory.** One table was called by all three tools: the first 1,000
markers of I195 × J432 (186 tetraploid samples, 2,486 MADC rows), target
records, default filters, on a four-core Intel Xeon at 2.1 GHz with 15 GB of
memory. Each command ran in its own session, and every process in that session
was sampled ten times a second: wall time, the CPU time of all its processes
together, and the peak of their summed proportional memory, which counts a
page shared by several processes once. polyRAD is `IterateHWE` with an
overdispersion of 9 followed by `GetProbableGenotypes`, and updog is
`multidog` with `model = "norm"`; both read the reference and alternate counts
`call` derives from the MADC rows. `IterateHWE` has no option to use more than
one core.

| Tool and settings | Cores | Wall | CPU | Peak memory |
|---|---|---|---|---|
| polyRAD 2.0.1, `IterateHWE` | 1 | 5 s | 5 s | 284 MB |
| `call` with polyRAD's settings (nothing estimated) | 1 | 6.5 s | 6 s | 56 MB |
| `call` with polyRAD's settings | 4 | 2 s | 7 s | 109 MB |
| `call --prior hwe` (the defaults) | 1 | 23 s | 22 s | 59 MB |
| `call --prior hwe` | 4 | 6.5 s | 23 s | 159 MB |
| `call --prior norm` | 1 | 24 s | 24 s | 58 MB |
| `call --prior norm` | 4 | 7 s | 25 s | 159 MB |
| `call --prior structure` | 1 | 30 s | 29 s | 67 MB |
| `call --prior structure` | 4 | 9 s | 32 s | 204 MB |
| `call --prior hwe --allele-bias fit` | 4 | 43 s | 160 s | 182 MB |
| `call --prior norm --allele-bias fit` | 1 | 170 s | 168 s | 64 MB |
| `call --prior norm --allele-bias fit` | 4 | 47 s | 170 s | 181 MB |
| updog 2.1.7, `norm` | 1 | 1,115 s | 1,108 s | 281 MB |
| updog 2.1.7, `norm` | 4 | 325 s | 1,128 s | 765 MB |

The estimate of the overdispersion and the background-read rate is a fixed
cost, because it works on at most 20,000 sample-by-marker cells however large
the table is (107 of the 1,000 markers here). The same runs on one core with
the two values given as numbers, so that nothing is estimated:

| Settings, values given | Wall | Share of the full run |
|---|---|---|
| `--prior hwe` | 4.5 s | 4.5 of 23 s |
| `--prior norm` | 7 s | 7 of 24 s |
| `--prior norm --allele-bias fit` | 57 s | 57 of 170 s |

- With polyRAD's settings `call` takes 6.5 s on one core where polyRAD takes
  5 s, in a fifth of the memory, and 2 s on four cores.
- A default run spends 18 of its 23 s on the estimate. Calling alone takes
  4.5 s with `hwe` and 7 s with `norm`.
- With the bias fitted the estimate takes about 115 s, since the bias of every
  sampled marker is fitted again for each pair of values tried. Calling with
  the bias takes 57 s, eight times as long as without.
- The whole run with the bias is 6.5 times as fast as updog on one core and 7
  times on four, with the accuracy shown above. Without the bias,
  `call --prior norm` is 46 times as fast as updog on one core and 45 times on
  four.
- Four workers give 3.1 to 3.6 times the speed of one. Memory with four
  workers is that of five Python processes, 110–205 MB.

The `call` rows were measured again after a review of the code for work done
more than once; polyRAD's and updog's rows are from the first measurement.
Before the review the default run took 33 s on one core and 10 s on four,
`norm` 38 s and 10 s, and `norm` with the bias 301 s and 80 s; calling alone
took 6 s, 10 s and 97 s. `call` with polyRAD's settings, whose code changed
least, took 7 s then and 6.5 s now, so the two measurements are of comparable
machines. Three changes account for the gain, and none changes a result (51
output files of twelve runs were the same byte for byte before and after):

- The likelihoods of a record's read counts are computed once and not in
  every pass. They depend on the allele frequencies only through polyRAD's
  contamination term, so they are computed in every pass only when
  `--contamination-rate` is above 0.
- The estimate chose its markers again for every background-read rate it
  tried, and now chooses them once.
- The binomial coefficient of a sample's counts and the terms that depend
  only on the overdispersion tried are computed once per sample and once per
  value, not once per dosage.

**Shortening the estimate when the bias is fitted: tried, and not done.** With
the bias fitted, two thirds of a run is the read-model estimate (113 of 170 s
above). On the benchmark table, in one process, the estimate takes 108 s
against 17 s without the bias: it works on 107 markers, tries five background
rates and passes over the markers nine times, and each pass fits every
marker's bias again, 963 fits of 0.1 s. A fit calls the record at about 21
bias values (17 to 43), and each call takes 3 passes with `hwe`, 4 with `norm`
and 6 with the free weights the estimate uses. A profile of the fit found no
step done twice: the time is in the passes themselves, not in the read-count
likelihoods, which are already computed once per bias value.

That leaves fitting the bias less often, which was tried in three ways:

| Change | Estimate on the benchmark table | What it cost |
|---|---|---|
| None | 108 s; 150 and 0.0005 | |
| Fit each marker's bias once, under the read model the search starts from | not run | The bias depends on the overdispersion it is fitted under: against the fit under the final read model, the log bias differed by a median of 0.14 and by more than 0.05 at 64 of the 107 markers (9 of 107 when fitted under an overdispersion of 100 where 150 is final) |
| Keep each marker's bias for an overdispersion when the search moves to another background rate | 86 s; 150 and 0.001 | A different background rate on both real tables (potato: 0.003 where it is 0.001), for a whole run that is 19% shorter on the benchmark table and 5% longer on potato. Potato accuracy was unchanged (96.59% against 96.57%) |
| Start the search from an estimate made on a quarter of the markers | 53 s; 150 and 0.0005 | No gain without the bias (17.3 s against 16.9 s); with an eighth of the markers, 91 s with the bias and 24 s without. The gain depends on the pilot landing near the answer |

None was kept. The first two change what is estimated, because two
neighbouring background rates often fit almost equally well and any change in
the biases decides between them; the third leaves the estimate alone but does
not reliably save time. A run with the bias is shortened by more workers, by
giving `--overdispersion` and `--background-rate` as numbers once a first run
has estimated them (57 s where the full run takes 170 s), or by a compiled
kernel, which is not planned.

**Keeping a fitted bias only where it earns its place: measured, not built.**
Fitting the bias raises agreement with the array by 4 points on the potato
data and lowers accuracy on tables without bias, by up to 8 points for
hexaploids pooled from real diploid reads (below). A rule that could be the
default was tried in a scratch copy of the module: a record keeps its fitted
bias only when the fit raises the log-likelihood, after the prior on the
bias, by a set margin over calling without one. Each data set was called with
`--prior hwe` under the read model its own run with the bias had estimated,
so that only the rule differs. Share of genotypes right, without filters:

| Data set | No bias (its own read model) | Bias always | Margin 2 | Margin 5 | Margin 10 | Margin 20 |
|---|---|---|---|---|---|---|
| Potato, 85 clones, against the array | 92.26% | 96.60% | 96.59% | 96.56% | 96.47% | 96.23% |
| Potato vendor MADC, diploids | 92.02% | 92.61% | 92.60% | 92.60% | 92.60% | 92.58% |
| Potato vendor MADC, tetraploids | 86.41% | 91.25% | 91.25% | 91.23% | 91.19% | 91.15% |
| Delta smelt targets (diploid) | 99.41% | 99.28% | 99.27% | 99.27% | 99.27% | 99.27% |
| Delta smelt off-target | 99.70% | 99.39% | 99.39% | 99.39% | 99.39% | 99.39% |
| Pooled tetraploids, 50 reads, targets | 90.56% | 89.65% | 90.17% | 90.46% | 90.53% | 90.64% |
| Pooled hexaploids, 50 reads, targets | 77.60% | 69.83% | 75.53% | 77.17% | 77.60% | 77.60% |
| Pooled hexaploids, off-target | 89.03% | 85.92% | 88.03% | 88.56% | 88.99% | 89.03% |

Records that keep a bias at margins 0, 2, 5, 10 and 20: potato 1,513, 1,153,
965, 789 and 590 of 1,544; the vendor MADC 3,851, 3,266, 2,949, 2,591 and
2,155 of 3,904; smelt 905, 330, 195, 137 and 84 of 946; pooled tetraploids
920, 123, 53, 29 and 19 of 946; pooled hexaploids 926, 89, 40, 26 and 18 of
944.

- A margin of 10 keeps the gain where there is bias (96.47% of 96.60%;
  91.19% of 91.25%) and removes the loss where there is none (90.53% against
  90.56%; 77.60% against 77.60%). On real unbiased reads nine records in ten
  give the bias up; on potato half to two thirds keep it.
- Smelt stays 0.14 points (targets) and 0.31 (off-target) under the run
  without bias at every margin. There the loss does not come with the many
  small biases; it comes with the read model the bias run estimated
  (overdispersion 200 where the run without bias has 150) or with the few
  records that keep a large bias.
- This is evidence for a third setting of `--allele-bias`, not a change:
  the repository's code is as it was, and whether to add the setting, and
  whether it should be the default, is open.

**The structure prior (`--prior structure`).** polyRAD's `IteratePopStruct`
gives each sample its own allele frequencies: it runs a principal component
analysis of the genotype estimates of all markers, regresses every allele on
the components, and repeats both with the calling until they settle. A port
would need a matrix library and would replace the marker-by-marker design
with an iteration over the whole table. `call` instead finds, once, the
samples most like each sample (`find_neighbours`) and then calls each record
on its own with Hardy-Weinberg proportions at the allele frequency among a
sample's neighbours: a second run of the record's iteration
(`_iterate_record`), started from the result of the first. Likeness is the
sum of squared
differences between two samples' expected allele shares over the markers,
from a first pass with the `hwe` prior, which ranks pairs as a genomic
relationship matrix would. The number of neighbours (0, 5, 10, 20, 40 or 80)
is the one whose neighbours, found on one half of the markers, best predict
the samples' shares at the other half; 0 is the `hwe` prior.

Test tables: tetraploids, 300 samples, 500 markers, read counts from
polyRAD's `SimAlleleDepth` (overdispersion 20, contamination and error rates
0.001), negative-binomial depths of mean 8, 25 or 150. The populations are
three groups of 100 with Balding-Nichols allele frequencies at three values
of FST; a cline in which every sample is a different mixture of two sources
(FST 0.3 between them); fifteen full-sib families of twenty; and one
random-mating population. The last two rows are the 85 potato clones above
with every read kept with probability 10/366 or 25/366, scored against the
array genotypes. polyRAD was given the simulated overdispersion (75 for the
potato rows) and `call` that value and a background rate of 0.001 (0.003).
Share of genotypes with reads whose most probable dosage is the true one:

| Population | Mean depth | polyRAD `IterateHWE` | polyRAD `IteratePopStruct` | `call` hwe | `call` structure | Neighbours chosen |
|---|---|---|---|---|---|---|
| Three groups, FST 0.05 | 8 | 71.2% | 71.3% | 71.1% | 71.9% | 80 |
| Three groups, FST 0.15 | 8 | 72.3% | 74.7% | 72.3% | 75.1% | 80 |
| Three groups, FST 0.15 | 25 | 83.2% | 84.6% | 83.4% | 84.8% | 80 |
| Three groups, FST 0.15 | 150 | 89.5% | 90.3% | 89.7% | 90.5% | 80 |
| Three groups, FST 0.30 | 8 | 73.9% | 79.0% | 73.9% | 79.4% | 80 |
| Three groups, FST 0.30 | 25 | 84.7% | 87.0% | 84.8% | 87.3% | 80 |
| Cline | 8 | 75.4% | 76.7% | 75.4% | 76.7% | 80 |
| Cline | 25 | 85.1% | 85.7% | 85.3% | 85.8% | 80 |
| Fifteen families | 8 | 70.8% | 73.1% | 70.8% | 73.8% | 20 |
| Fifteen families | 25 | 81.9% | 83.4% | 82.0% | 83.9% | 20 |
| One population | 8 | 70.5% | 67.1% | 70.5% | 70.5% | 0 |
| Potato families, thinned | 10 | 71.0% | 71.1% | 70.9% | 71.8% | 10 |
| Potato families, thinned | 25 | 80.8% | 81.4% | 80.7% | 81.4% | 10 |

- The structure prior calls as many genotypes correctly as `IteratePopStruct`
  on every table or more, by up to 0.7 points. The cline is the case
  principal components describe best, and the two are level there.
- On the six tables with three groups it equals what calling each group on
  its own with the `hwe` prior gives (72.0%, 75.1%, 84.8%, 90.5%, 79.4% and
  87.3%), which is the most a prior can gain from knowing the groups.
- With one population `IteratePopStruct` kept 51 components and lost 3.4
  points. The structure prior chose 0 neighbours.
- The gain over `hwe` falls with depth: 0.8 to 5.5 points at 8 reads, 0.5 to
  2.5 at 25, and 0.8 at 150. At full depth the potato clones go from 92.3% to
  92.3% (10 neighbours), and the two alfalfa populations from 90.87% and
  90.04% agreement with the published dosages to 90.89% and 89.98% (20
  neighbours each).
- With the read model estimated instead of given, three of the tables that
  were run again gave the same shares (79.4%, 73.8% and 71.8%).

Three ways of using the same relationships were compared on a test bench
with one likelihood and one first pass (not in the repository; it used
numpy): the nearest samples, as above; principal components, their number
chosen on held-out markers; and a ridge predictor from the whole relationship
matrix, its shrinkage chosen the same way. On the group and cline tables the
three were within 0.3 points of each other. On the families the nearest
samples gave 73.8% and 83.9%, the components 73.1% and 83.8%, and the ridge
predictor 72.7% and 83.4%. On the thinned potato counts the ridge predictor
failed (33.9% and 56.0%, with the least shrinkage offered), and the other two
were level (71.9% and 81.3% against 71.7% and 81.4%). A sample's own value
was left out of its own prediction in all three.

The number of neighbours, components or the shrinkage has to be chosen on
markers that took no part in measuring the relationships. Chosen on the same
markers, the bench took the least shrinkage and the most components it was
offered (44.1% and 72.9% on the first table, against 72.3% for `hwe`),
because each marker is then one of the features its own prediction is made
from.

On the benchmark table of "Time and memory" the structure prior took 30 s on
one worker against 23 s for `hwe`, and 9 s against 6.5 s on four, with 204 MB
against 159 MB; 18 s of both is the read-model estimate. With the read model
given it took 10.6 s against 4.0 s on one worker: the first pass over the
markers that the samples are compared on, the comparison, and the second run
of every record. (An earlier version of this paragraph gave 34.6 s against
33.4 s for the whole run, which understated the difference: that code,
measured again with the read model given, took 12.5 s against 5.1 s.) On the
simulated tables, 500 markers and 300 samples, it took 4 to 5 s against 2 s
on four workers with the read model given.

Not tested for accuracy: ploidies other than four, samples of mixed ploidy,
and read counts from other simulators. A run with diploid and tetraploid
samples and all three VCF types completes, and `bcftools view` reads its
output. The prior is not combined with the allelic bias or the selfing rate;
`call` refuses both combinations.

**Family priors, measured and not built.** updog's `f1` model and polyRAD's
`PipelineMapping2Parents` give the offspring of a cross the dosages their
parents' genotypes can segregate into. Whether `call` needs such a prior was
measured with a test bench (not in the repository) that uses the `call`
likelihood and, per family and marker, the pair of parental dosages that
makes the reads most probable, under polysomic inheritance without double
reduction and with 1% of the prior spread evenly over the dosages for
off-types.

One biparental population with both parents among the samples is the common
case. The two alfalfa populations have that; their reads were also thinned,
parents included, and the calls of the 974 and 986 published markers compared
with the published dosages, which come from the full depth:

| Population | Mean depth | hwe | norm | structure | F1, both parents' reads | F1, one parent's reads | F1, offspring only |
|---|---|---|---|---|---|---|---|
| I195 × J432 | 174 | 90.87% | 90.11% | 90.85% | 91.01% | 90.82% | 90.64% |
| I195 × J432 | 25 | 84.53% | 84.08% | 84.60% | 85.14% | 84.95% | 84.79% |
| I195 × J432 | 10 | 75.72% | 76.04% | 75.74% | 76.73% | 76.63% | 76.58% |
| E115 × F628 | 424 | 90.04% | 89.80% | 89.98% | 90.33% | 90.26% | 90.15% |
| E115 × F628 | 25 | 82.23% | 82.25% | 82.25% | 82.84% | 82.78% | 82.71% |
| E115 × F628 | 10 | 74.15% | 74.68% | 74.33% | 75.05% | 74.99% | 74.94% |

With both parents' reads the F1 prior is 0.1 to 0.3 points above the best of
the three priors `call` has at full depth, 0.5 to 0.6 at 25 reads and 0.4 to
0.7 at 10. The parents' reads account for 0.1 to 0.4 of that; the rest comes
from the segregation shape, which a family of this size shows by itself.

Many small families are the other case. On the potato half-diallel (ten
families of 1 to 20) and on fifteen simulated full-sib families of 20:

| Data | Mean depth | hwe | structure | F1, parents estimated from offspring | F1, parental genotypes given |
|---|---|---|---|---|---|
| Simulated families | 8 | 70.8% | 73.8% | 75.0% | |
| Simulated families | 25 | 82.0% | 83.9% | 84.9% | |
| Potato, thinned | 10 | 70.9% | 71.8% | 71.8% | 76.1% |
| Potato, thinned | 25 | 80.7% | 81.4% | 81.1% | 84.1% |
| Potato, full depth | about 300 | 92.3% | 92.3% | 92.3% | 93.4% |

Only parental genotypes known from elsewhere help much here; those given were
the article's phased parents, from array data and linkage analysis, which are
better than a parent's own reads would give. That needs a pedigree, which is
often not available, and is not built. A single biparental population gains
under a point from an F1 prior, so that is not built either.

The biparental case was also simulated, with the parent models of the other
two tools: polyRAD's `PipelineMapping2Parents` (F1, no linkage, given the
overdispersion `call` estimated) and updog's `f1` model with both parents
named, against updog's `norm`. Tetraploid, 200 progeny and both parents, 300
markers; genotypes and reads from polyRAD's simulators (`SimGenotypesMapping`
and `SimAlleleDepth`, overdispersion 20, no bias) or updog's (`rgeno` and
`rflexdog`, sequencing error 0.5%, overdispersion 49–499 and a bias of
exp(N(0, 0.5)) per marker). The parents were sequenced like the progeny.
Share of the progeny's genotypes called correctly:

| Reads from | Mean depth | `call` hwe | `call` norm | `call` norm + bias | F1 prior on `call`'s likelihood | polyRAD mapping | updog `norm` | updog `f1` |
|---|---|---|---|---|---|---|---|---|
| polyRAD | 8 | 74.9% | 77.1% | | 77.5% | 77.3% | 74.0% | 75.4% |
| polyRAD | 25 | 84.5% | 85.7% | | 86.0% | 85.8% | 80.6% | 82.3% |
| polyRAD | 150 | 90.3% | 91.2% | | 91.3% | 91.2% | 88.9% | 89.5% |
| updog | 8 | 69.3% | 70.1% | 75.3% | 70.3% | 70.6% | 75.5% | 74.4% |
| updog | 25 | 79.4% | 79.1% | 87.8% | 80.1% | 80.0% | 87.6% | 87.6% |
| updog | 150 | 85.2% | 84.5% | 97.1% | 85.5% | 84.9% | 97.1% | 97.3% |

On unbiased reads the `norm` prior is within 0.4 points of the F1 prior and
of polyRAD's mapping prior: it fits the narrow spread of one family's
dosages. On biased reads `call` with `norm` and the bias fitted is within 0.3
points of updog's `norm` and `f1` at 25 and 150 reads. The structure prior
gave the same shares as `hwe` on all six tables. polyRAD's mapping prior left
13–16% of the genotypes of the biased tables uncalled.

With 10% of the progeny replaced by unrelated samples (25 reads), true
progeny and off-types apart:

| Reads from | Caller | True progeny | Off-types |
|---|---|---|---|
| polyRAD | `call` hwe | 84.6% | 81.3% |
| polyRAD | `call` norm | 85.8% | 79.4% |
| polyRAD | F1 prior on `call`'s likelihood | 86.2% | 76.2% |
| polyRAD | polyRAD `IterateHWE` | 84.4% | 81.4% |
| polyRAD | polyRAD mapping | 86.1% | 70.8% |
| polyRAD | updog `norm` | 79.9% | 78.9% |
| polyRAD | updog `f1` | 82.4% | 75.6% |
| updog | `call` norm + bias | 87.6% | 82.0% |
| updog | updog `norm` | 87.6% | 82.6% |
| updog | updog `f1` | 87.3% | 79.0% |
| updog | polyRAD mapping | 78.1% | 65.6% |

A family prior costs 3 to 11 points on the samples that are not of the
family, because it leaves them only the family's dosage classes, and gains
at most 2.5 points on those that are.

The same parent models were run on the alfalfa populations, with the parents
named and polyRAD given the overdispersion `call` estimated. Agreement with
the published dosages, which were themselves called by updog from the full
depth; for polyRAD's mapping prior the share of genotypes it called is in
brackets, and its agreement is among those:

| Population | Mean depth | `call` hwe | `call` norm | `call` norm + bias | polyRAD `IterateHWE` | polyRAD mapping | updog `norm` | updog `f1` |
|---|---|---|---|---|---|---|---|---|
| I195 × J432 | 174 | 90.87% | 90.11% | 92.27% | 90.78% | 90.83% (92.7%) | 92.98% | 93.29% |
| I195 × J432 | 25 | 84.53% | 84.08% | 85.77% | 84.44% | 85.12% (90.0%) | 86.26% | 86.61% |
| I195 × J432 | 10 | 75.72% | 76.04% | 76.03% | 75.52% | 77.02% (87.5%) | 75.87% | 75.56% |
| E115 × F628 | 424 | 90.04% | 89.80% | 92.00% | 89.93% | 90.42% (93.7%) | 92.68% | 93.06% |

E115 × F628 at 25 and 10 reads was not run: the job reached its time limit
after the four rows above. updog's `f1` model is 0.3 to 0.4 points above its
`norm` at full depth and at 25 reads, and 0.3 below at 10. polyRAD's mapping
prior is 0.05 to 1.5 points above its `IterateHWE` on the genotypes it calls
and leaves 6 to 13% uncalled. `call` with `norm` and the bias is 0.5 to 0.7
points below updog's `norm` at full depth and at 25 reads and level with it at
10. The truth here favours updog, whose calls it is.

**Priors on haplotype records.** Neither polyRAD's paper nor its package
gives accuracy figures for multi-allele (haplotype) calls: the paper scores
two-allele SNPs with array truth, and the package has import examples and
small demonstration datasets without expected genotypes. On polyRAD's own
DArTag example (`inst/extdata/DArTag_example.csv`, two markers with three and
five haplotypes, four tetraploid samples) `call --vcf haplotype` with
polyRAD's settings gives polyRAD's genotype for all eight sample-by-marker
genotypes.

Accuracy was measured on microhaplotype tables from polyRAD's simulators
(`SimGenotypes`, `SimGenotypesMapping`, `SimAlleleDepth`; overdispersion 20,
contamination and error rates 0.001): 300 markers of 3 to 6 alleles, 300
samples, as one population, three groups (FST 0.15) or fifteen full-sib
families of 20. polyRAD was given the overdispersion, and `call` that and a
background rate of 0.001. Share of genotypes with reads in which the dosage
of every allele is right, without filters:

| Ploidy | Population | Mean depth | polyRAD `IterateHWE` | polyRAD `IteratePopStruct` | `call` hwe | `call` norm | `call` structure | `call` uniform |
|---|---|---|---|---|---|---|---|---|
| 4 | One population | 25 | 75.0% | stopped with an error | 74.6% | 74.6% | 74.6% (0 neighbours) | 70.9% |
| 4 | Three groups | 8 | 55.2% | 58.9% | 55.1% | 56.8% | 59.3% | 53.6% |
| 4 | Three groups | 25 | 75.1% | 77.0% | 75.0% | 75.5% | 77.0% | 73.5% |
| 4 | Three groups | 150 | 85.8% | 86.8% | 85.7% | 85.8% | 86.7% | 84.5% |
| 4 | Families | 8 | 54.7% | 56.8% | 54.7% | 54.8% | 58.6% | 48.5% |
| 4 | Families | 25 | 75.0% | 76.8% | 74.6% | 74.9% | 77.1% | 71.0% |
| 2 | Three groups | 8 | 92.5% | 93.0% | 92.6% | 92.7% | 93.0% | 92.2% |
| 2 | Families | 8 | 91.8% | 92.3% | 91.8% | 92.1% | 92.6% | 90.9% |

- `hwe` is within 0.4 points of polyRAD's HWE prior, as it should be.
- `structure` is never below `hwe` and gains 1 to 4 points in tetraploids
  with groups or families; it is within 0.1 points of `IteratePopStruct` or
  above it on every table, by 1.8 points on the low-depth families. With one
  population it chooses 0 neighbours.
- `norm` is never below `hwe` and gains at most 1.7 points. It takes two to
  five times as long, and 2 and 11 of the 300 records of the two low-depth
  tetraploid tables reached the iteration limit.
- `uniform` is the lowest on every table, by up to 6 points.

With the default filter (`--min-marginal-posterior 0.9`), 94 to 99% of the
genotypes kept are right under each of `hwe`, `norm` and `structure`. The
share kept at this overdispersion is 4 to 37% in tetraploids and 74 to 83% in
diploids, and is 3 to 8 points higher with `structure` than with `hwe` where
there are groups or families.

Not covered: read counts from other simulators, real haplotype data with
known genotypes, and hexaploids.

**Off-target records (`--vcf all`).** The same eight tables give truth for
the off-target records, which no earlier comparison scored: every haplotype
other than Ref and Alt carries one substitution of its own, so the true
dosage of a substitution is that of its haplotype. `call --vcf all` was run
with its defaults, so the overdispersion and background rate were estimated
(18 to 20 where 20 was simulated, and 30 on the tetraploid families at 8
reads; 0.0002 to 0.001). Share of genotypes with reads called correctly, for
target and off-target records, without filters; and for `hwe` with the
default filters (`--min-dp 10`, `--min-gq 10`), the share of the off-target
genotypes kept that are right and the share kept:

| Ploidy | Population | Mean depth | Off-target records | hwe: target, off-target | structure: target, off-target | norm: target, off-target | hwe, default filters: right, kept |
|---|---|---|---|---|---|---|---|
| 4 | One population | 25 | 760 | 86.8%, 88.2% | 86.8%, 88.2% | 86.7%, 88.2% | 97.8%, 61% |
| 4 | Three groups | 8 | 758 | 78.4%, 80.5% | 80.6%, 82.2% | 79.2%, 81.0% | 96.8%, 22% |
| 4 | Three groups | 25 | 759 | 87.0%, 88.6% | 88.2%, 89.6% | 87.4%, 88.9% | 98.1%, 62% |
| 4 | Three groups | 150 | 760 | 92.2%, 93.1% | 92.9%, 93.7% | 92.4%, 93.2% | 98.8%, 79% |
| 4 | Families | 8 | 758 | 77.8%, 79.7% | 79.7%, 81.3% | 77.8%, 79.8% | 96.2%, 22% |
| 4 | Families | 25 | 760 | 86.9%, 88.2% | 88.0%, 89.2% | 86.9%, 88.3% | 97.8%, 61% |
| 2 | Three groups | 8 | 747 | 96.5%, 96.9% | 96.7%, 97.1% | 96.6%, 97.0% | 99.8%, 33% |
| 2 | Families | 8 | 750 | 96.4%, 96.6% | 96.8%, 96.9% | 96.5%, 96.7% | 99.8%, 33% |

- Off-target records are called as well as target records: 0.2 to 2.1 points
  more of their genotypes are right. Their alternate alleles are a little
  rarer in these tables (mean frequency 0.14 to 0.23, against 0.25 for the
  target's).
- The priors rank as they do on target records: `structure` gains up to 2.2
  points on targets and 1.7 on off-target records where there are groups or
  families and nothing with one population, and `norm` is within 0.1 points
  below and 0.8 above `hwe`.
- With the default filters 96 to 99.8% of the off-target genotypes kept are
  right.

**Off-target substitutions shared by haplotypes.** In the tables above a
substitution is on one haplotype only and never on the Alt tag. In real tags
several haplotypes can carry the same one, and their reads are then added
together for its record. The same reads and true haplotype dosages were given
other stand-in tags: one substitution is on two of the haplotypes with the
reference base at the target, a second is on two as well, and a third is on a
haplotype that has the target's alternate base, so that the target's record
adds two haplotypes too. The true dosage of a substitution is the sum over the
haplotypes that carry it. Of the 760 off-target records 311 add two
haplotypes, and 149 of the 300 target records do. Share of genotypes with
reads called correctly, as in the table above:

| Ploidy | Population | Mean depth | hwe: target, off-target | structure: target, off-target | norm: target, off-target | hwe, default filters: target right, kept | hwe, default filters: off-target right, kept |
|---|---|---|---|---|---|---|---|
| 4 | One population | 25 | 82.8%, 85.4% | 82.8%, 85.4% | 82.8%, 85.4% | 97.0%, 47% | 97.4%, 53% |
| 4 | Three groups | 8 | 73.4%, 76.5% | 75.9%, 78.8% | 74.2%, 77.3% | 96.2%, 16% | 96.6%, 19% |
| 4 | Three groups | 25 | 83.9%, 86.1% | 85.2%, 87.3% | 84.3%, 86.5% | 97.6%, 50% | 97.8%, 56% |
| 4 | Three groups | 150 | 90.2%, 91.6% | 91.0%, 92.2% | 90.3%, 91.7% | 98.3%, 69% | 98.6%, 73% |
| 4 | Families | 8 | 71.9%, 75.7% | 74.9%, 77.9% | 72.1%, 75.8% | 95.9%, 15% | 96.0%, 18% |
| 4 | Families | 25 | 83.0%, 85.5% | 84.6%, 86.8% | 83.0%, 85.6% | 97.0%, 47% | 97.4%, 54% |
| 2 | Three groups | 8 | 95.8%, 96.3% | 96.0%, 96.5% | 95.8%, 96.3% | 99.8%, 33% | 99.8%, 33% |
| 2 | Families | 8 | 95.4%, 96.0% | 96.0%, 96.4% | 95.5%, 96.0% | 99.7%, 33% | 99.7%, 33% |

- A record that adds the reads of several haplotypes is the same call as a
  two-allele marker with those reads. The reads of the carriers and of the
  other haplotypes were added up outside `call` and called as two-row markers
  under the same overdispersion and background rate: all 2,538,900 genotypes
  of the eight tables are the same as in the `all` VCF.
- Fewer genotypes are right than in the tables with one haplotype per
  substitution (by 1 to 6 points for tetraploids, under 1 for diploids)
  because the added alleles are more frequent: a mean of 0.33 against 0.25 at
  targets and 0.29 against 0.21 at off-target records, so more genotypes have
  a dosage in the middle, where neighbouring dosages are harder to tell
  apart.
- Off-target records are still called as well as target records, the priors
  rank as before, and with the default filters 96 to 99.8% of the genotypes
  kept are right.

Real off-target genotypes are tested in the delta smelt comparison below.

**The polyRAD paper's own benchmark.** Clark et al. (2019) tested polyRAD on
read depths simulated from real genotypes: SNP-array and GoldenGate genotypes
of diversity panels and mapping populations, with reads drawn binomially at a
mean depth of about 10 (15 to 21% of the cells have no reads). The data, the
scripts and the results are in the Illinois Data Bank
(doi:10.13012/B2IDB-9729830_V2), with the error of every marker under each
model of polyRAD 0.4-0 and under EBG and GATK-style calls. `call` was run on
the same read depths and scored as the paper did: the root mean squared error
between true and estimated copy number at each marker, over the cells with
reads, for the posterior mean ("continuous") and for the most probable
genotype ("discrete"), averaged over the markers. Lower is better. `call`
estimated an overdispersion of 500, the highest it tries, as it should for
binomial reads (100 on the soybean set). Continuous and discrete error:

| Data set | Samples × markers | polyRAD HWE (archive) | `call` hwe | polyRAD PopStruct (archive) | `call` structure | `call` norm |
|---|---|---|---|---|---|---|
| *Miscanthus sinensis* diversity panel, 2x | 565 × 395 | 0.154, 0.175 | 0.154, 0.175 | 0.140, 0.161 | 0.140, 0.162 | 0.150, 0.170 |
| The same, simulated tetraploids | 565 × 395 | 0.420, 0.464 | 0.421, 0.464 | 0.338, 0.384 | 0.326, 0.372 | 0.398, 0.457 |
| Apple diversity panel, 2x | 96 × 3,704 | 0.173, 0.197 | 0.173, 0.197 | 0.176, 0.200 | 0.172, 0.196 | 0.173, 0.196 |
| Potato diversity panel, 4x | 221 × 3,763 | 0.439, 0.499 | 0.439, 0.499 | 0.427, 0.491 | 0.428, 0.488 | 0.438, 0.502 |
| *Glycine soja*, 2x, selfing rate 0.95 | 1,179 × 2,957 | 0.037, 0.034 | 0.037, 0.034 | 0.038, 0.037 | not with a selfing rate | |

Two sets are F1 populations, which polyRAD called with its mapping prior and
the parents named; `call` has no family prior:

| Data set | Samples × markers scored | polyRAD mapping (archive) | `call` hwe | `call` norm |
|---|---|---|---|---|
| *Miscanthus* F1, 2x | 85 × 185 | 0.186, 0.205 | 0.200, 0.226 | 0.186, 0.213 |
| Potato F1, 4x | 238 × 2,306 | 0.391, 0.447 | 0.409, 0.478 | 0.385, 0.446 |

- `call` with the `hwe` prior has polyRAD's published HWE error to the third
  decimal on all five diversity sets, soybean with `--selfing-rate 0.95`.
- `structure` equals polyRAD's population-structure model on *Miscanthus* and
  potato and is better on the simulated tetraploids (0.326 against 0.338) and
  on apple, where polyRAD's model was worse than its HWE.
- `norm`, told nothing about the cross, matches polyRAD's mapping prior on
  the *Miscanthus* F1 and is 0.006 below it on the potato F1 in continuous
  error; in discrete error it is within 0.008.
- polyRAD's models that use linkage between markers are better than `call` on
  five of the seven sets, by 0.003 to 0.044 in continuous error: 0.137 on
  *Miscanthus* (`call` at its best 0.140), 0.412 on potato (0.428), 0.164 on
  apple (0.172), and 0.170 and 0.340 on the two F1 populations (0.186 and
  0.385). `call` treats markers independently. On the simulated tetraploids
  `structure` is below polyRAD's structure-and-linkage model (0.326 against
  0.335), and on soybean linkage made polyRAD worse.
- Without the selfing rate the soybean errors are 0.157 with `hwe`, 0.124 with
  `norm` and 0.110 with `structure`, against 0.037 with it: for inbred
  samples the selfing rate has to be given.
- EBG's HWE model, the other caller in the archive with a population prior,
  has the discrete error of `call` with `hwe` (0.464 against 0.464 on the
  simulated tetraploids, 0.500 against 0.499 on potato), as two HWE callers
  should on binomial reads. The GATK-style calls, which have no population
  prior, are far above both on the polyploid sets (0.523 and 0.696).
- With polyRAD's present defaults (overdispersion 9) `call` is a little
  worse on these binomial reads, 0.458 against 0.421 on the simulated
  tetraploids.

The archive's eighth set is a simulated allohexaploid (soybean genotypes on
one of three subgenomes: 0 to 2 copies of the allele out of 6, mean depth 30).
polyRAD modelled it as three diploid genomes and has 0.194 and 0.219 with HWE
and 0.178 and 0.199 with population structure. `call` has no such model and
called it as an autohexaploid: 0.360 and 0.429 with `hwe` (82.1% of genotypes
right), 0.312 and 0.365 with `norm` (87.3%), 0.336 and 0.402 with `structure`.
Markers that sit on one subgenome of an allopolyploid are outside what `call`
models, because it spreads the prior over all 0 to 6 copies.

**Delta smelt GT-seq against alignment-based genotypes (diploid).** The
example data of gtseq_microhap (Zenodo, doi:10.5281/zenodo.19069550,
CC-BY-4.0) are 96 delta smelt, 410 loci, 150,000 read pairs per fish, with the
panel's primer file and a VCF of 2,787 biallelic SNPs called from the same
reads by BWA and bcftools on the reference genome. `generate` has been
measured on them (see "External benchmarks on public datasets"); here they
test `call` on real diploid reads, at targets and at off-target SNPs, against
genotypes from an alignment-based caller.

Every locus was placed on the genome (GCA_021917145.1) by its two primers,
all 410 once. For the 368 loci with a polymorphic SNP of the comparison VCF in
the first 81 bases, the tag is the genome sequence from the forward primer and
the target is the SNP with the highest minor allele frequency; `generate` on
the R1 reads gave 1,527 allele rows and assigned about 100,000 reads per
fish, some 280 per marker. `call --vcf all --ploidy 2` wrote 368 target records and 615
off-target records; 467 of the off-target records are a SNP of the comparison
VCF with the same two bases (111 are not in it, 36 have more than two alleles
and 1 has other bases). Agreement with the comparison genotype:

| Caller | Targets (35,073 genotypes) | Off-target SNPs (44,777) |
|---|---|---|
| `call --prior hwe` (the defaults; overdispersion 150, background 0.003) | 99.57% of the 98.8% called | 99.80% of the 99.1% called |
| `call --prior norm` | 99.58% of 98.8% | 99.81% of 99.1% |
| `call --prior structure` (5 neighbours) | 99.58% of 98.8% | 99.80% of 99.1% |
| `call --prior hwe --allele-bias fit` | 99.42% of 98.8% | 99.48% of 99.1% |
| `call` with polyRAD's settings | 99.20% of 98.1% | 99.43% of 98.7% |
| `call --prior hwe`, no filters | 99.41% | 99.70% |
| polyRAD `IterateHWE`, overdispersion 9 | 98.71% | 99.08% |
| polyRAD `IteratePopStruct`, overdispersion 9 | 98.85% | 99.13% |
| polyRAD `IterateHWE`, overdispersion 150 | 99.38% of the 99.7% called | 99.66% |
| updog `norm` | 99.27% | 99.38% |

Every row is scored on the same genotypes: those where the sample has reads
at the record and the comparison genotype rests on 10 or more reads. polyRAD
and updog also give a genotype where a sample has no reads; those are left
out.

polyRAD and updog were given the reference and alternate counts of the 947
two-allele records and call every genotype, so the row to compare them with is
`call` without filters.

- On real diploid reads `call` with its defaults agrees with the
  alignment-based genotypes more often than polyRAD or updog do on the same
  counts: by 0.14 points on targets and 0.32 on off-target SNPs over updog,
  and by 0.70 and 0.62 over polyRAD's defaults.
- Off-target SNPs are called as well as targets, here on real reads; the
  simulated tables under "Off-target records" said the same.
- Fitting the bias lowers the agreement, by 0.15 points on targets and 0.3 on
  off-target SNPs. updog, which always fits one (median 0.98 here; 176 of 946
  outside 0.8 to 1.25), is also below `call` without it. On the potato data
  fitting the bias raised agreement with the array by 4.3 points. Whether it
  helps depends on the data, which is a reason to leave it off by default.
- polyRAD's population-structure model stopped with an error at an
  overdispersion of 150, as it did on the potato data at 75.
- The comparison genotypes come from the same reads, so this measures
  agreement between callers and not error against an independent assay.

**Audit against both packages' sources.** Every function of `madcats_call.py`
was read against polyRAD 2.0.1's and updog 2.1.7's source and run against the
packages on the same inputs, to settle what is adopted, what is adapted and
what is MADCATS's own. The module's opening comment and the functions'
comments now say so, with the full references.

| Compared | Result |
|---|---|
| Six simulated records (ploidy 2, 4, 6 and mixed; 2 to 4 alleles; one with selfing) against `IterateHWE` | likelihoods, posteriors, frequencies and corrected genotypes agree to 1e-13 |
| 366 random corrections of dosages to the ploidy against `CorrectGenos` | 366 the same |
| The same where no step has any probability | `call` gives no call; polyRAD returned all-zero dosages in 33 cases and did not return in one |
| Ties between dosages | `BestGenos` compares in single precision and takes the first: posteriors 0.4, 0.4, 0.2 give dosage 0 in polyRAD and no call in `call` |
| Prior on the bias, bias formula at error rate 0, limits of `norm` | identical to updog's `pen_bias`, `xi_fun` and bounds |
| 40 single `norm` fits against updog's | 39 agree to 4e-5; the fortieth is the stalled fit described under "The norm prior", since fixed |

Three statements were wrong and are corrected: this record listed the
free-weights prior as MADCATS's own (it is updog's `flex`; only its use for
the estimate is new), the command's description said "The model is
polyRAD's", and a comment described the frequency floor of 1 / (2 × samples)
as polyRAD's for a diploid population (polyRAD divides by the largest locus
ploidy it is given, which is 2 in the stored comparisons, whatever the
samples' ploidy).

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

The genotype model of `call` is adapted from polyRAD and updog (see "Genotype
calling" for which part comes from where):

- Clark LV, Lipka AE, Sacks EJ (2019) polyRAD: genotype calling with
  uncertainty from sequencing data in polyploids and diploids. G3
  9(3):663-673. <https://doi.org/10.1534/g3.118.200913>
- Gerard D, Ferrão LFV, Garcia AAF, Stephens M (2018) Genotyping polyploids
  from messy sequencing data. Genetics 210(3):789-807.
  <https://doi.org/10.1534/genetics.118.301468>
- Gerard D, Ferrão LFV (2020) Priors for genotyping polyploids.
  Bioinformatics 36(6):1795-1800.
  <https://doi.org/10.1093/bioinformatics/btz852>
- de Silva HN, Hall AJ, Rikkerink E, McNeilage MA, Fraser LG (2005)
  Estimation of allele frequencies in polyploids under certain patterns of
  inheritance. Heredity 95:327-334. <https://doi.org/10.1038/sj.hdy.6800728>
