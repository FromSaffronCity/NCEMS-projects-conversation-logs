# Phaser verification on NA12878: SHAPEIT5 → RefLinker with single-cell Hi-C
*Status report, 2026-10-20 (HapCUT2 results added; all rates re-verified with WhatsHap)*

**In short.** We phase NA12878's het SNPs the way we plan to phase the brain donors, and score the result against her pedigree-phased truth.
- SHAPEIT5 gets the local phase right but not the chromosome-scale phase.
- RefLinker (single-cell Hi-C) fixes the chromosome-scale phase but adds small local errors.
- A simple **hybrid** step keeps the best of both.
- HapCUT2 (WGS reads) adds a small further gain.

Proposed pipeline: **SHAPEIT5 → RefLinker (default) → hybrid → HapCUT2**.

## 1. Key findings
1. **SHAPEIT5 alone is locally accurate but has no chromosome-scale phase.** It has 0.4% switch error, yet ~45% of het SNPs sit on the wrong haplotype across whole chromosomes.
2. **RefLinker (Hi-C) restores chromosome-scale phase** in both independent datasets¹. Wrong-haplotype SNPs drop to 7.0% (Dip-C, 17 cells) and 2.1% (LiMCA, 220 cells), improving steadily with cell number.
3. **RefLinker leaves many tiny wrong-haplotype islands** (median 2 SNPs), which triple the switch error. Two fixes remove them and also improve whole-chromosome accuracy:
   - a **hybrid** post-processing step (keep SHAPEIT5's local phase, take block orientation from RefLinker), which reaches **0.6%** wrong-haplotype SNPs at 220 cells;
   - a **stricter RefLinker flip threshold**.
4. **Removing duplicate and low-quality Hi-C reads does not help.** The islands come from RefLinker's decisions, not from read noise.
5. **Without HapCUT2, error-free segments are no longer than SHAPEIT5's** (longest ~11.5 Mb, N50 ~1.3 Mb).
   - The breaks that remain are SHAPEIT5's own errors, which Hi-C cannot fix.
   - These breaks (non-flip switches) stay at ~0.1% of SNP pairs in every configuration except default RefLinker.
6. **Read-based patching (HapCUT2, NA12878 WGS ~50×)** adds a small gain:
   - **On its own:** HapCUT2 phases 69.2% of het SNPs, in small blocks (block-span N50 0.6 kb; switch error within blocks 0.12%).
   - **Added to the hybrid (even chromosomes):**
     - error-free N50 1.26 → 1.38 Mb (Dip-C) and 1.40 → 1.58 Mb (LiMCA 220);
     - longest error-free segment 11.55 → 13.41 Mb in both;
     - Hamming 0.0522 → 0.0521 (Dip-C) and 0.0058 → 0.0056 (LiMCA 220).

## 2. How the results are scored
All metrics compare a phasing with the truth, allowing one global swap of the two haplotypes. Only het SNPs present in both the phasing and the truth are scored. Each chromosome counts as one block.

| Metric | What it counts | Scale |
|---|---|---|
| **Hamming error** | SNPs on the wrong haplotype ÷ SNPs scored (0.5 = random) | whole chromosome |
| **Switch error** | switches ÷ adjacent SNP pairs scored. A switch is a pair of neighbouring SNPs whose relative phase is wrong. A single wrongly phased SNP (a *flip*) counts as 2 switches | local |
| **Flip error** | single wrongly phased SNPs ÷ adjacent SNP pairs (WhatsHap's "switch/flip decomposition") | local |
| **Non-flip switch error** | (switches − 2 × flips) ÷ adjacent SNP pairs. These are real block breaks: after one, a whole stretch is on the wrong haplotype (more severe) | local |
| **Longest error-free segment**, **error-free N50** | length (Mb) of the stretches between any two errors (flips included) | — |

- **Same SNPs:** within each table row group, every configuration is scored on the same SNPs (those RefLinker phased), so the rows are directly comparable.
- **Independent check:** `whatshap compare` (WhatsHap 2.8) reproduces every switch, flip, non-flip switch and Hamming rate in both results tables on the same SNPs (all 20 rows; differences ≤ 0.00001, from rounding only).
- **Fair testing:** the settings of the hybrid (window) and of tuning (threshold) are chosen on odd chromosomes. Results are reported on even chromosomes.

## 3. Pipeline and data
| Step | Tool | Input | Output |
|---|---|---|---|
| 1. Statistical phasing | SHAPEIT5 5.1.1 | unphased genotypes, reference panel, genetic map | complete phasing, one chromosome-long haplotype pair |
| 2. Hi-C read processing | `pre-meta` → `bwa mem -5SP` | raw single-cell Hi-C reads | merged BAM per dataset |
| 3. Hi-C refinement | RefLinker (`extract`, `pop`) | step-1 phasing + BAM | re-oriented phasing of Hi-C-linked SNPs |
| 4. Hybrid | our script | steps 1 + 3 | SHAPEIT5 phase with RefLinker's block orientation (majority vote over a window of SNPs) |
| 5. Read-based patching | HapCUT2 1.3.4 | WGS reads + genotypes | small read-backed blocks merged into step 4 |

- **Truth:** NA12878 pedigree-phased het SNPs (Platinum Genomes 2016-1.0, hg19).
  - 2,083,940 SNPs on chr1–22 and 63,598 on chrX.
  - Independent of any statistical or Hi-C phasing.
- **SHAPEIT5 input:** the truth het sites with the phase removed (plus hom-alt SNPs), mimicking a WGS genotype call set like the brain donors'.
  - Panel: 1000G phase3 without NA12878 and her relatives (2,501 samples; chrX: 1,269 women).
  - Genetic maps: b37.
- **Hi-C reads:** raw reads only, aligned by us to hg19. Two independent datasets¹:
  - Dip-C (Tan et al. 2018): 17 cells, 854 M read pairs.
  - LiMCA (Wu et al. 2024): 220 cells, 4,200 M read pairs, also run as nested subsets of 17, 50 and 100 cells.
- **WGS reads for HapCUT2:** NA12878 Platinum Genomes WGS (~50×). The brain donors' raw WGS reads are also available (Salk server).

## 4. Results

### 4.1 Hi-C restores chromosome-scale phase (SHAPEIT5 → RefLinker, default settings, chr1–22)
On the same SNPs, SHAPEIT5 alone has Hamming 0.446–0.447, switch error 0.38%, flip error 0.13–0.14% and non-flip switch error 0.11%.

| Hi-C data | het SNPs phased by RefLinker | Hamming | switch error | flip error | non-flip switch error | error-free N50 |
|---|---|---|---|---|---|---|
| Dip-C, 17 cells | 95.1% | 0.070 | 1.06% | 0.14% | 0.78% | 0.41 Mb |
| LiMCA, 17 / 50 / 100 / 220 cells | 87.4 / 94.5 / 94.9 / 95.2% | 0.217 / 0.057 / 0.035 / 0.021 | 0.87–1.29% | 0.13–0.28% | 0.61–0.83% | 0.36–0.53 Mb |

- Hamming falls steeply with more cells. Switch error rises about threefold because of RefLinker's tiny islands (fixed in 4.2).
- chrX shows the same pattern (Hamming 0.46 → 0.07 with Dip-C; 0.05 → 0.03 with 100–220 LiMCA cells). The exception is 17 LiMCA cells, which were too few on chrX (no gain).

![Metrics vs number of LiMCA cells](figures/fig1_cells_curve.png)
*Fig. 1. LiMCA: share of het SNPs phased, Hamming error, longest error-free segment and error-free N50 vs number of cells. Blue: + RefLinker. Orange: SHAPEIT5 alone, scored on the same SNPs. The SHAPEIT5 phasing itself is identical at every cell count; its orange values move only because the set of scored SNPs (those RefLinker phased) grows with cell number.*

### 4.2 Fixing RefLinker's islands, and adding HapCUT2 (even chromosomes, same SNPs in every row)
Each cell gives Hamming / switch / flip / non-flip switch / error-free N50 / longest error-free segment.

| Configuration | Dip-C, 17 cells | LiMCA, 220 cells |
|---|---|---|
| SHAPEIT5 alone | 0.442 / 0.35% / 0.13% / 0.10% / 1.28 / 11.5 Mb | 0.442 / 0.35% / 0.13% / 0.10% / 1.28 / 11.5 Mb |
| + RefLinker (default) | 0.064 / 1.03% / 0.13% / 0.76% / 0.41 / 11.5 Mb | 0.019 / 1.25% / 0.27% / 0.70% / 0.37 / 11.6 Mb |
| + RefLinker (tuned threshold²) | **0.028** / 0.39% / 0.13% / 0.14% / 1.24 / 11.5 Mb | 0.0075 / 0.37% / 0.14% / 0.10% / 1.35 / 11.5 Mb |
| + RefLinker + hybrid³ | 0.052 / 0.38% / 0.13% / 0.13% / 1.26 / 11.5 Mb | 0.0058 / 0.38% / 0.14% / 0.11% / 1.40 / 11.5 Mb |
| + RefLinker + hybrid + HapCUT2 | 0.052 / **0.34%** / **0.10%** / 0.14% / **1.38** / **13.4** Mb | **0.0056** / **0.34%** / **0.11%** / 0.11% / **1.58** / **13.4** Mb |

² Flip threshold `-e` −50 (Dip-C) and −100 (LiMCA 220), instead of the README default −10.
³ Window of 1,001 (Dip-C) and 101 (LiMCA 220) SNPs. Any window from 3 to 201 SNPs improves on default RefLinker in every dataset.

**Reading the table, step by step:**
- **RefLinker (default):** Hamming drops sharply. Its extra errors are mostly short wrong stretches of ≥2 SNPs, so they show up as non-flip switches (0.70–0.76%) and the N50 shrinks.
- **Hybrid or tuned threshold:** both remove those stretches, bringing switch error back to SHAPEIT5's level and lowering Hamming further.
- **HapCUT2:** lowers the flip rate (0.13–0.14% → 0.10–0.11%), which accounts for its switch-error gain. It also lengthens error-free segments.
- **Why segments stop growing:** apart from default RefLinker, non-flip switches (real block breaks) stay at ~0.1%.
  - RefLinker re-orients the long blocks between SHAPEIT5's breaks, which is what drops Hamming, but it rarely removes the breaks themselves.
  - The hybrid keeps 741 of SHAPEIT5's 995 non-flip switches (LiMCA 220, even chromosomes) and adds new ones a median of 5 SNPs (~21 kb) away.
  - Long-range errors thus become short misphased stretches.
- **Cleaner reads give no gain:**
  - About 53–54% of Hi-C alignments were duplicates; only ~37% of alignments remain after deduplication and MAPQ ≥ 20.
  - Re-running RefLinker on them leaves Hamming essentially unchanged (Dip-C 0.064 → 0.065; LiMCA 220 0.019 → 0.019).

![All configurations, Dip-C](figures/fig3_pipeline_steps_all17.png)
![All configurations, LiMCA 220](figures/fig3_pipeline_steps_limca_n220.png)
*Fig. 2. The four metrics for every configuration (even chromosomes, same SNPs per row), Dip-C (top) and LiMCA 220 cells (bottom). Blue: proposed pipeline SHAPEIT5 → RefLinker → hybrid → HapCUT2. Also shown: HapCUT2 merged onto SHAPEIT5 alone and onto default RefLinker.*

## 5. Practical notes
- **RefLinker** needed three memory-safety patches to run on sparse single-cell Hi-C (hard-clipped reads; chromosomes with few linked SNPs). They are included in our build.
- **Transfer to brain data:**
  - **Hybrid over tuning:** tuned thresholds are data-specific (the best value grew from −20 to −100 with cell number), and there is no truth to tune on for the donors. The hybrid is more robust: a fixed window of ~100–200 SNPs worked everywhere.
  - **Bisulfite reads:** the donors' snm3C-seq Hi-C reads are bisulfite-converted. Before RefLinker they need bisulfite-aware alignment (e.g. Bismark), plus masking of conversion-ambiguous bases using the donors' WGS het SNPs.
- **SNPs not phased:** 3.9% of truth het SNPs are absent from the panel, so SHAPEIT5 can't phase them. Of the rest, ~1% have no Hi-C link.
- **Cost:** ~4.2k core-hours so far. RefLinker needs up to ~52 GB of memory per large chromosome at 220 cells.

## 6. Next steps (for discussion)
1. Decide whether read-based patching (HapCUT2) is worth adding, based on the last row of the table in 4.2 and Fig. 2.
2. Choose the final pipeline: SHAPEIT5 → RefLinker (default) → hybrid (fixed window) [→ HapCUT2, if it helps].
3. Brain donors:
   - hg38 panel and maps;
   - bisulfite-aware alignment and masking of the snm3C-seq reads;
   - download the donors' WGS reads if HapCUT2 helps.

¹ *The two datasets differ in depth, read length, lab and cell clone, so they are independent reference points and are not compared with each other.*

*Code, logs and data: `/anvil/projects/x-bio260019/x-aahmed7/phaser-verification/` on Purdue Anvil.*
