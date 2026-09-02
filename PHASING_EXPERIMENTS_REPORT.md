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
