# Phasing Heterozygous SNPs from Single-Cell Bisulfite Sequencing — Project Report

**Period:** 2026-07-12 → 2026-09-02   **Donors:** H1930001, H1930002, H1930004 (human cortex, CX45)
**Data:** Tian W, Zhou J, … Ecker JR, *Science* 2023 (DOI 10.1126/science.adf5357), GEO GSE215353,
per-cell BAMs from the NeMO Archive. No re-mapping — the authors' alignments are used as published.

---

## 1. Objective

Recover **haplotypes** — which of a person's two chromosome copies each allele sits on — from
single-cell **bisulfite** sequencing, where the C→T conversion breaks ordinary SNP callers, and
where each cell alone is far too shallow. The approach: merge many single cells into a pseudo-bulk
BAM, call SNPs with a bisulfite-aware caller, and phase the heterozygous ones with **HapCUT2**,
which links two het sites whenever a single read or read pair covers both.

The project asked three questions in sequence, and answered all three:

1. Does statistical phasing work here? **No.**
2. Does Hi-C (snM3C-seq) fix the read-span bottleneck that limits snMC-seq? **Yes, decisively.**
3. What actually limits phasing once Hi-C is in place? **Per-site read depth** — not homozygosity,
   and not block gaps.

## 2. Pipeline

`download per-cell BAMs → merge to pseudo-bulk → (snM3C only) repair read pairing →
bisulfite-aware SNP calling → HapCUT2 preprocessing → extractHAIRS → HapCUT2`

Two SNP callers were run on every dataset: **bsgenova** (posterior model) and a **naive** caller
with fixed thresholds. Both apply a hard **DP ≥ 10** cutoff. Reference: hg38 + lambda. Full
pipeline runtime at the largest scale is ~26 h per donor.

Two snM3C-specific fixes were required and are the reason the experiment succeeded:
`repair_m3c_pairs.py` (the merged 3C BAMs contain **no paired reads at all** as delivered), and
`extractHAIRS --hic 1` together with `HAPCUT2 --hic 1`.

## 3. Experiments

| assay | tiers | phasing mode | pseudo-bulk input |
|---|---|---|---|
| snMC-seq | 10-cell ×3 donors; 200-cell (d2, d4); 1000-cell (d1) | standard | d1 1000c: 1,454,425,153 reads / 123.5 GB · d2 200c: 225,256,293 / 22.9 GB · d4 200c: 357,187,360 / 33.4 GB |
| snMC-seq downsample | donor-1 1000-cell BAM at 0.02–0.70 | standard | — |
| snM3C-seq | 10-cell ×3 donors | `--hic 1` | — |
| snM3C-seq | 100-cell ×3 donors | `--hic 1` | d1: 843,150,557 reads / 45.5 GB · d2: 868,532,873 / 44.5 GB · d4: 780,730,691 / 40.9 GB |

36 caller × dataset combinations in total.

## 4. Genotyping and phasing results

### snM3C-seq (`--hic 1`)

| donor | cells | caller | raw calls | het calls | het offered | **phased** | rate |
|---|---|---|---|---|---|---|---|
| H1930001 | 10 | bsgenova | 186,988 | 165,413 | 163,894 | 48,239 | 0.294 |
| H1930001 | 10 | naive | 73,089 | 61,426 | 60,099 | 20,192 | 0.336 |
| H1930002 | 10 | bsgenova | 302,301 | 263,936 | 262,266 | 81,033 | 0.309 |
| H1930002 | 10 | naive | 123,590 | 101,713 | 99,905 | 31,374 | 0.314 |
| H1930004 | 10 | bsgenova | 230,426 | 199,739 | 198,106 | 61,774 | 0.312 |
| H1930004 | 10 | naive | 93,657 | 77,093 | 75,502 | 23,809 | 0.315 |
| H1930001 | 100 | bsgenova | 4,939,635 | 3,919,539 | 3,897,860 | 2,239,973 | 0.575 |
| H1930001 | 100 | naive | 5,242,025 | 4,431,872 | 4,389,230 | **2,939,059** | 0.670 |
| H1930002 | 100 | bsgenova | 5,278,393 | 4,227,826 | 4,203,619 | 2,459,141 | 0.585 |
| H1930002 | 100 | naive | 5,771,712 | 4,930,636 | 4,885,323 | **3,338,060** | 0.683 |
| H1930004 | 100 | bsgenova | 4,995,524 | 4,029,953 | 4,013,332 | 2,351,363 | 0.586 |
| H1930004 | 100 | naive | 4,884,362 | 4,134,601 | 4,093,725 | **2,708,901** | 0.662 |

### snMC-seq (standard mode)

| donor | cells | caller | raw calls | het calls | het offered | phased | rate |
|---|---|---|---|---|---|---|---|
| H1930001 | 10 | bsgenova | 1,166 | 869 | 792 | 361 | 0.456 |
| H1930001 | 10 | naive | 660 | 440 | 419 | 216 | 0.516 |
| H1930002 | 10 | bsgenova | 421 | 292 | 271 | 130 | 0.480 |
| H1930002 | 10 | naive | 306 | 190 | 180 | 106 | 0.589 |
| H1930004 | 10 | bsgenova | 1,712 | 1,378 | 1,211 | 384 | 0.317 |
| H1930004 | 10 | naive | 688 | 485 | 457 | 250 | 0.547 |
| H1930002 | 200 | bsgenova | 1,255,862 | 844,504 | 847,372 | 137,155 | 0.162 |
| H1930002 | 200 | naive | 410,865 | 210,915 | 200,326 | 21,057 | 0.105 |
| H1930004 | 200 | bsgenova | 2,268,560 | 1,426,975 | 1,432,952 | 305,197 | 0.213 |
| H1930004 | 200 | naive | 1,260,616 | 711,759 | 708,448 | 117,414 | 0.166 |
| H1930001 | 1000 | bsgenova | 3,749,921 | 2,319,448 | 2,314,396 | 976,750 | 0.422 |
| H1930001 | 1000 | naive | 3,423,882 | 2,040,577 | 2,012,384 | 694,795 | 0.345 |

### Donor-1 downsample series (fractions of the 1000-cell BAM, bsgenova / naive phased)

| fraction | 0.02 | 0.05 | 0.10 | 0.20 | 0.40 | 0.70 | 1.00 |
|---|---|---|---|---|---|---|---|
| bsgenova phased | 319 | 2,372 | 37,331 | 239,048 | 498,146 | 788,324 | 976,750 |
| bsgenova rate | 0.466 | 0.181 | 0.127 | 0.192 | 0.282 | 0.374 | 0.422 |
| naive phased | 209 | 894 | 4,826 | 60,888 | 368,204 | 606,465 | 694,795 |
| naive rate | 0.571 | 0.378 | 0.137 | 0.138 | 0.259 | 0.322 | 0.345 |

(`frac0.02` is degenerate — 684 het sites, all in collapsed repeats — and is not a data point.)

**"Het calls" vs "het offered"** differ in both directions: HapCUT2 preprocessing drops 0.2–1.0 % of
het sites *and* splits multi-allelic sites across records. Phased counts join on `(chrom, pos)` and
run ~0.16 % above the per-record totals in the runner logs.

## 5. Haplotype block structure

| dataset | caller | largest block span | blocks ≥ 1 Mb |
|---|---|---|---|
| snM3C 10-cell (d1 / d2 / d4) | bsgenova | 83.8 / 71.1 / 69.7 Mb | 284 / 576 / 504 |
| snM3C 10-cell (d1 / d2 / d4) | naive | 37.0 / 42.2 / 56.6 Mb | 98 / 212 / 125 |
| snM3C 100-cell (d1 / d2 / d4) | bsgenova | 127.2 / 141.1 / 150.1 Mb | 85,866 / 93,133 / 104,017 |
| snM3C 100-cell (d1 / d2 / d4) | naive | 155.4 / 163.9 / 174.8 Mb | 123,566 / 139,511 / 133,633 |
| snM3C 10-cell, **without** `--hic 1` | bsgenova | **897 bp** | 0 |
| snMC 1000-cell | bsgenova | **2,742 bp** | 0 |

The last two rows are the headline: the same data phased without Hi-C mode reaches ~900 bp,
and 1000 snMC cells at three times the per-site depth reach 2.7 kb. Span N50 at 100 cells is
~22–24 Mb.

**Blocks are a sparse overlapping skeleton, not a tiling.** The ≥1 Mb blocks sum to ~1,400 Gb of
span against a 3.1 Gb genome — they overlap 300–450 deep, and most phased SNPs still sit in small
local blocks nested inside them. More depth buys *more* Mb-scale blocks, not longer ones.

## 6. Validation (no truth haplotype exists for these donors)

| check | pairs | discordance |
|---|---|---|
| same caller, split-half, ≥ 1 Mb | 10,182 | **18.48 %** |
| same caller, split-half, ≥ 10 Mb | 8,292 | 18.52 % |
| cross-caller, ≥ 1 Mb (one pair per anchor) | 53,539 | **28.01 %** |
| null control (alleles randomised) | 10,182 | 50.18 % |

Short-range phasing is in good shape (6.73 % switch rate below 1 kb at 100 cells). Long-range
agreement is real signal but noisy, and flat from 1 Mb to 10 Mb — the Hi-C signature. These are
**internal consistency bounds, not switch-error rates**; no absolute accuracy can be quoted.

## 7. Coverage and homozygosity analysis

Measured from the published callsets and `.blocks` files alone (both carry per-site `DP`), across
all 36 runs:

1. **Phasing rate is a monotonic function of a site's own depth** — ~46 % at DP 10–14, ~100 % by
   DP 150. The single `coverage-per-variant` figure in `HAPCUT2.log` averages this away.
2. **Naive's win at high depth is positional, not qualitative.** At equal depth the callers are
   within 2–6 points; ~77 % of naive's lead is having more calls at greater depth. Replicated
   across three donors (22.6 / 22.9 / 18.8 % attributable to per-site efficiency).
3. **The depth curve transfers within a regime only:** across donors 0.4–2.7 % error, across assays
   200 %, across depths without Hi-C 28–67 %. Without long-range links, a site's phasing depends on
   its *neighbours'* depth, not its own.
4. **There are no runs of homozygosity.** Longest in any callset is 500 kb; real ROH is multi-Mb.
   The hypothesis that homozygous stretches cause block boundaries is dead.
5. **The two assays fail in opposite ways.** With `--hic 1`, 99.8 % of unphased het sites lie
   *inside* a block's span (inter-block gaps total ~35 Mb). Without it, 99.8 % lie *between* blocks,
   in gaps totalling ~2.95 Gb — effectively the whole genome.
6. **Cross-donor differences in snMC-seq are coverage, not biology.** Donor-2 reaches DP ≥ 20 at
   only 5.0 % of sites versus donor-4's 31.7 %; just above a hard cutoff, callset size is a steep
   function of depth. An earlier reading of donor-4 as intrinsically richer was withdrawn.
7. **Caveat:** extreme-depth pileups (DP p99 up to 13,605) are collapsed satellite repeats. A
   `DP > 200` ceiling would be a cheap, defensible filter on future runs.

Per-site depth medians for orientation: snM3C 100-cell 21–28; snMC 1000-cell 65; snMC 200-cell
**13–20** — i.e. 200 snMC cells are much *shallower* per site than 100 snM3C cells. Cell count is
not depth.

### Depth vs. linkage — the two regimes, mechanistically

The two assays are limited by different things, and the whole difference reduces to one
distinction: **what a site's phasing depends on.**

- **snMC-seq (standard) is linkage-limited.** A het site phases only if a read or read pair
  spans it *and* an adjacent het, so its fate depends on its **neighbours'** depth. An isolated
  deep site with shallow neighbours strands; 99.8 % of unphased het sites lie *between* blocks.
- **snM3C-seq (`--hic 1`) is depth-limited.** A long-range Hi-C contact can link a site from
  anywhere on the chromosome, so its fate depends on its **own** depth; 99.8 % of unphased het
  sites lie *inside* block spans.

Hi-C linkage therefore does not merely add information — it **decouples** a site's phasing from
its local neighbourhood. This is why, at equal per-site depth, snM3C-seq phases a higher
*fraction* of het sites (compare rate, not raw counts — absolute totals also ride on callset
size). The evidence is indirect but strong: the depth→phasing-rate curve transfers across donors
(0.4–2.7 % error) yet breaks across assays (~200 %, §7 item 3). That ~200 % gap *is* the linkage
premium, quantified.

**Not yet directly measured.** The equal-depth claim rests on the depth curve failing to transfer,
not on a controlled comparison. Binning het sites by DP in both assays and comparing phasing rate
bin-for-bin would turn "strongly inferred" into "directly measured" — cheap, since every callset
and `.blocks` file already carries per-site DP.

## 8. Conclusions

1. **Statistical phasing (SHAPEIT5) does not work** for this data — the calls are low-coverage
   single-cell variants largely absent from population reference panels. Read-backed phasing is the
   only usable route.
2. **Hi-C linkage is worth more than depth.** 100 snM3C cells phase ~67 % of het sites into
   chromosome-scale blocks; 1000 snMC cells at 3× the per-site depth phase 42 % into kb-scale blocks.
3. **Caller ranking tracks the assay/depth regime, never the donor** — bsgenova wins at low depth
   (6.5× and 2.6× more phased SNPs at 200 cells), naive at high depth (+15–36 % at 100 cells). Three
   donors agree within each regime. Always run both; never carry a ranking to a new experiment.
4. **What limits phasing is per-site depth in the Hi-C regime, and linkage in the standard regime.**
   Adding cells helps by lifting individual sites over the DP threshold.
5. Best current result: **3,338,060 phased het SNPs** (H1930002, 100-cell snM3C, naive), largest
   block 163.9 Mb.

## 9. Limitations

- No truth haplotype → no absolute switch-error rate, only the bounds in §6.
- ~1 in 5 to 1 in 4 long-range relationships flips between independent halves or callers; treat
  megabase-scale phase as provisional.
- 10-cell tiers are tiny in snMC-seq (190–1,378 het sites) — a pipeline smoke test, not evidence
  about callers.
- Three donors, one brain region, one study. Nothing here is a population claim.

## 10. Where things live

| what | where |
|---|---|
| All published callsets, phasings, BAMs | `workspace/data/` (see its `README.md` for the full layout) |
| Coverage / homozygosity outputs and scripts | `workspace/coverage-homozygosity-computation/` |
| Analysis and pipeline scripts | `workspace/Python-scripts/`, `workspace/shell-scripts/` |
| Full session-by-session narrative | `snMC-seq_SNP_phasing_experiment.md`, `snM3C-seq_SNP_phasing_experiment.md` |

## 11. Validation phase (planned — 2026-09-14, not yet started)

The CX45 results above are to be validated two ways (full plans in the two session logs' 2026-09-14
entries).

1. **Cross-region replication.** Rerun the full pipeline on **two additional brain regions** (beyond
   CX45) to test whether the trends here reproduce region-to-region: per region, snMC-seq pseudo-bulk
   at **10 and 200 nuclei** (standard mode), and snM3C-seq pseudo-bulk at **10 nuclei** plus a large
   tier whose **read depth is matched to the 200-nuclei snMC-seq** run (rather than the fixed 100
   nuclei used at CX45) — i.e. a deliberately **depth-matched** snMC-vs-snM3C comparison, which the
   CX45 work lacked (§7). Both callers throughout; specific regions TBD.
2. **External genotype truth set.** The Science-paper authors provided genotyping data for these
   donors (details TBD). It will be characterised, then used to validate the phasing — potentially
   yielding genotype-concordance for the callers and, if the authors' data is phased, the first
   **absolute switch-error rate** for our haplotypes, which §6 and §9 note has never been available.

## 12. Validation against WGS ground-truth genotypes (2026-09-17)

Author whole-genome-sequencing genotypes (Salk/BICAN, `neomorph.salk.edu/ftp/bican/WGS`, one
single-sample `HBAgenomics` VCF per donor, GATK `GT:AD:DP:GQ:PL`, **unphased**) give the first
external truth for this project. Each donor's WGS holds **~1.94 M true heterozygous SNPs**
(H1930001 1,942,674 · H1930002 1,946,484 · H1930004 1,922,030; ~3.37 M total SNPs).

**Method.** Our preprocessed het SNP calls vs the WGS via `bcftools norm -m -any` + `isec` on
`(chrom,pos,ref,alt)`. **precision** = our-het-also-het-in-WGS / our-total-het (treats WGS-absent
sites as hom-ref, i.e. a *lower bound* — a minority of "absent" could be WGS no-calls). **recall** =
our-confirmed-het / WGS-het. "Concord% (shared)" = of our het calls that coincide with a WGS *variant*
site, the fraction WGS also calls het.

### 12a. Genotype concordance — our het calls vs WGS truth

| assay/tier | donor | caller | our het | shared w/ WGS variants | confirmed het (TP) | **precision** | **recall** | concord%(shared) |
|---|---|---|---|---|---|---|---|---|
| snMC 1000-cell | H1930001 | bsgenova | 2,314,396 | 1,749,111 | 1,727,147 | **74.6%** | **88.9%** | 98.7% |
| snMC 1000-cell | H1930001 | naive | 2,012,384 | 1,735,969 | 1,725,607 | **85.7%** | **88.8%** | 99.4% |
| snMC 200-cell | H1930002 | bsgenova | 847,372 | 791,708 | 505,553 | 59.7% | 26.0% | 63.9% |
| snMC 200-cell | H1930002 | naive | 200,326 | 171,248 | 161,941 | 80.8% | 8.3% | 94.6% |
| snMC 200-cell | H1930004 | bsgenova | 1,432,952 | 1,317,033 | 1,010,267 | 70.5% | 52.6% | 76.7% |
| snMC 200-cell | H1930004 | naive | 708,448 | 629,164 | 589,677 | 83.2% | 30.7% | 93.7% |
| snM3C 100-cell | H1930001 | bsgenova | 3,897,860 | 1,369,689 | 1,166,617 | 29.9% | 60.1% | 85.2% |
| snM3C 100-cell | H1930001 | naive | 4,389,230 | 944,928 | 911,787 | 20.8% | 46.9% | 96.5% |
| snM3C 100-cell | H1930002 | bsgenova | 4,203,619 | 1,369,558 | 1,181,735 | 28.1% | 60.7% | 86.3% |
| snM3C 100-cell | H1930002 | naive | 4,885,323 | 946,417 | 925,283 | 18.9% | 47.5% | 97.8% |
| snM3C 100-cell | H1930004 | bsgenova | 4,013,332 | 1,288,213 | 1,084,016 | 27.0% | 56.4% | 84.1% |
| snM3C 100-cell | H1930004 | naive | 4,093,725 | 847,643 | 818,100 | 20.0% | 42.6% | 96.5% |

### 12b. Findings

1. **Depth dominates genotype accuracy.** The deep **snMC 1000-cell** pseudo-bulk (~58×) recovers
   **~89 % of the donor's true het SNPs at 75–86 % precision** — the pipeline is *positively
   validated at depth*. Shallow pseudo-bulks are far worse: snMC 200-cell (~11–17×) reaches only
   8–53 % recall, and snM3C 100-cell (~3×) sits at 20–30 % precision. Both callers **over-call at low
   depth** — snM3C reports 3.9–4.9 M hets against ~1.94 M real, so most shallow-callset hets are
   false positives.

2. **Precision/recall is a caller trade-off whose winner is depth/regime-dependent** — the accuracy
   analogue of the phasing-yield rule (§8.3). **naive is consistently the more precise caller in
   snMC** (81–86 % vs 60–75 %) and bsgenova the more sensitive (always higher recall). But in the
   shallow **snM3C** regime this inverts: bsgenova is *both* more precise (27–30 % vs 19–21 %) and
   more sensitive. Best single callset: **snMC 1000-cell naive — 85.7 % precision, 88.8 % recall.**

3. **What the false positives are — and it is NOT what was first assumed.** The WGS-absent ("private")
   hets are mostly false, but their substitution spectrum splits by depth. At high depth
   (snMC 1000-cell) bsgenova's private hets are **70.7 % C>T/G>A** — classic **bisulfite C→T
   conversion artifacts**; naive's conversion-masking (the hickit rule) strips these to **24.8 %**,
   which is exactly why naive out-precises bsgenova here. At low depth (snM3C 100-cell) private hets
   are *depleted* in C>T/G>A (7–17 %) — the FPs there are **non-conversion mapping/low-depth errors**,
   which masking cannot catch. So bisulfite conversion FP is a *high-depth* problem the masking solves;
   shallow FP is a different, unsolved problem.

4. **Consequence for phasing.** At 100-cell snM3C, naive *phases the most* SNPs (§8.3) but has the
   *lowest* genotype precision (~20 %) — a large share of its phased hets are false. The deep snMC
   1000-cell callset is genotype-accurate (~85 %) but phases only into kb-scale blocks (no Hi-C).
   Neither current operating point is both accurate and long-range; the ideal is deep **snM3C**
   (accurate calls + Hi-C linkage), not yet run.

### 12c. Caveats

- WGS is variant-only, so precision treats every WGS-absent site as hom-ref: it is a **lower bound**;
  a definitive FP-vs-no-call split needs the WGS BAM/gVCF (not held).
- Cross-donor snMC rows are at unequal depth (d1 1000-cell vs d2/d4 200-cell) — compare within tier.
- This validates **genotype calling**. The **switch-error** validation of the *haplotypes* (statistically
  phasing the WGS into a truth and scoring our HapCUT2 blocks against it) is the next step, underway.

## 13. Follow-up validation (2026-09-21): switch-error, FP mechanism, CX47 replication

### 13a. Switch-error vs a WGS truth — the first *absolute* haplotype-accuracy number

We built a phased truth by statistically phasing each donor's WGS genotypes with **SHAPEIT5
`phase_common`** against the 1000G high-coverage reference panel (streamed) + b38 genetic map, then
scored our HapCUT2 haplotypes against it. Switch error = within-block adjacent het pairs (het in both,
phased in ours) whose cis/trans orientation disagrees with the truth (orientation-invariant, so it is
insensitive to each block's arbitrary global flip). First pass: **donor H1930001, chr20** (truth =
43,465 phased het sites).

| phasing (CX45, d1, chr20) | shared phased het | comparable pairs | switches | **switch error** |
|---|---|---|---|---|
| snMC 1000-cell bsgenova | 13,645 | 7,458 | 810 | **10.86 %** |
| snMC 1000-cell naive | 13,601 | 7,516 | 816 | **10.86 %** |
| snM3C 100-cell bsgenova (`--hic`) | 17,629 | 8,977 | 2,105 | **23.45 %** |
| snM3C 100-cell naive (`--hic`) | 16,172 | 7,600 | 2,003 | **26.36 %** |

**Findings.** (1) This is the **first switch error against an external truth** — previously only
internal bounds existed (§6: ~18.5 % same-caller split-half, ~28 % cross-caller at ≥1 Mb). The snM3C
numbers (23–26 %) land squarely inside those bounds, **externally validating them.** (2) **Deep
short-range snMC (≈11 %) is ~2× more accurate than long-range snM3C Hi-C (23–26 %)** — the long-range
contacts that give snM3C its Mb-scale blocks also carry substantially more switch noise, confirming
the "Mb blocks are real but provisional" caveat (§9). (3) Even the best case (~11 %) is far above a
production WGS phasing (<1–2 %), so these pseudo-bulk haplotypes remain research-grade.
**Caveats:** chr20 / donor-1 only; the truth is itself *statistically* phased so it carries its own
error (our numbers are a combined upper bound); extension to more chromosomes/donors pending.

### 13b. What the false positives actually are — refined

§12b(3) is sharpened by measuring the private-FP substitution spectrum **and depth**:
- **High depth (snMC 1000-cell):** private FPs are **70.7 % C>T/G>A** → bisulfite **C→T conversion
  artifacts**; naive's masking strips them to 24.8 % (why naive out-precises bsgenova at depth). Solid.
- **Low depth (snM3C 100-cell):** the FPs are **neither** conversion **nor** shallow-depth noise —
  their median DP is **20** (≈ the true hets' 23), and their spectrum is strongly skewed to
  **T>C/A>G (41 %)** with C>T *depleted* (16 %), vs a balanced ~33 %/33 % in true hets. So they are a
  **structured, systematic artifact**, not random error — most plausibly the **bisulfite
  strand-masking asymmetry** (the caller masks Watson-T and Crick-A, which can distort the A/G–T/C
  balance). Confirming the mechanism needs the per-strand Watson/Crick counts (`DPW`/`DPC`) — open.

### 13c. CX47 regional replication (validation Task 1, snMC 200-cell)

CX47 snMC phasing is complete for all 3 donors (10 & 200-cell, both callers). 200-cell phased SNPs:

| donor | bsgenova | naive |
|---|---|---|
| H1930001 | 501,284 | 343,299 |
| H1930002 | 284,533 | 94,919 |
| H1930004 | 397,788 | 237,093 |

**The CX45 trends replicate:** bsgenova > naive at 200-cell snMC in **every** donor (1.5–3.0×),
phasing rate ~25–28 %, blocks short/read-length-scale (~2.4 SNPs/block, no long-range) — as expected
for standard-mode snMC. CX47's absolute yields run higher than CX45's, but this is **depth**, not
region biology: at the same 200-cell count the CX47 pool carries more reads (d2 321 M vs 225 M; d4
465 M vs 357 M) and higher per-site depth (d4 median DP **22 vs 17**), because the largest-200 CX47
cells are deeper-sequenced libraries. So absolute yields are depth-confounded across regions — the
**robust replication is the caller-ordering trend**, not the raw counts. (CB63 and CX46 200-tiers
were ~5 % short of the cell threshold and are completing after a tolerance relaxation, 2026-09-21.)

## 14. Cross-region snMC phasing complete (2026-09-22) — trend replicates in 3 regions incl. non-cortex

All snMC regional phasing (validation Task 1) is finished: **CX47** (3 donors), **CX46** (donors 1–2),
**CB63** (donor-4, *cerebellum* — non-cortex), each 10 + 200-cell × both callers. 200-cell phased-SNP
counts:

| region | donor | bsgenova | naive |
|---|---|---|---|
| CX47 (cortex) | H1930001 | 501,284 | 343,299 |
| CX47 | H1930002 | 284,533 | 94,919 |
| CX47 | H1930004 | 397,788 | 237,093 |
| CX46 (cortex) | H1930001 | 303,346 | 112,931 |
| CX46 | H1930002 | 386,413 | 216,958 |
| CB63 (cerebellum) | H1930004 | 366,083 | 194,643 |

**bsgenova > naive at 200-cell snMC in every region×donor combination (6/6), including the non-cortex
cerebellar region CB63** — the CX45 caller-ordering trend (§12/§8.3) replicates robustly across brain
regions and tissue type. As at CX45, absolute yields track pseudo-bulk depth (per-region cell-library
depth), so the reproducible finding is the **caller ordering**, not the raw counts. (CX46/CB63 200-tiers
required relaxing the merge/availability tolerance to 190/200 cells — a few per-region cells were
unrecoverable mount objects or, for CB63, download-failed; negligible for the pseudo-bulk.)

**FP composition, made explicit (donor-1 CX45 snMC 1000-cell, from §12).** Of bsgenova's 2,314,396 het
calls: **75.6 % overlap the WGS variant set** (74.6 % correct het + ~1 % miscalled genotype) and
**24.4 % are private/not-in-WGS — overwhelmingly false positives** (70.7 % of them are C>T/G>A
bisulfite-conversion artifacts); a genuinely rare/de-novo real fraction cannot be separated from a
variant-only WGS but is a small minority. naive is cleaner (86.3 % overlap, 13.7 % private) because
its conversion-masking removes the C>T FPs.
