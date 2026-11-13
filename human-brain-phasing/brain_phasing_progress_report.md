# Phasing the human brain donors: progress report
*2026-11-13 · data: Tian et al. 2023 snM3C-seq + donor WGS · code and logs: Anvil `x-aahmed7/human-brain-phasing/`*

## In short
- **Goal:** phase each donor's WGS het SNPs, then label every snM3C-seq contact end with the chromosome copy it came from
  (0 / 1), for ChromoGen.
- **Pipeline** (verified on NA12878, see the [pipeline guide](../phaser-pipeline-guide/phaser_pipeline_guide.md)):
  SHAPEIT5 → RefLinker → hybrid refiner, with HapCUT2 in parallel → combiner → pruner.
- **Done so far: donor 1 (H1930001)**, phased with 4, 100 and 2,370 of its own snM3C cells, with and without pruning.
  Contacts of the 4 pilot cells in x-tzahin's `pilot_v2/original.h5` are labelled.
- **Main results** (truth-free check on 20 held-out cells; on the authors' contacts, the share of same-chromosome
  contacts whose two ends get the same label at > 10 Mb):
  - SHAPEIT5 alone: 49% (coin flip).
  - 4 cells: no gain.
  - 100 cells: 86%.
  - **2,370 cells + pruning: 88%** (95–96% at 0.1–10 Mb).
- **Labels per cell stay sparse:** ~7% of contact ends and ~12% of contacts get a label, whatever the phasing. Labels
  come only from reads that cover a het SNP, so **phasing accuracy is what we can improve.**
- **In progress:** donors 1, 2 and 4 with ~5,000 cells each, then labelling two new contact files (350 cells). Currently
  blocked by the exhausted Slurm allocation.

## 1. Data and setup
- **Donor 1 WGS** (GATK, hg38): 1,927,881 het + 1,356,652 hom-alt SNPs (PASS, biallelic, chr1–22).
- **Reference panel:** 1000G 30x phased, 3,202 people; SHAPEIT4 b38 maps.
- **snM3C-seq reads:** the authors' per-cell Bismark BAMs from NeMO. Donor 1 has 49,480 snM3C cells (S5 table).
- **Cell sets** (cell type ignored: every cell of a donor carries the same SNPs):

  | Set | Cells | How chosen |
  |---|---|---|
  | n4 | 4 | the 4 donor-1 pilot cells (V1C, L4-IT) of `pilot_v2/original.h5` |
  | n100 | 100 | the 4 pilot cells + the 96 batch-1 donor-1 cells with the most contacts |
  | n2370 | 2,370 | the 4 pilot cells + 5% (2,474) random donor-1 cells; 2,370 downloaded (104 failed) |
  | held-out | 20 | random donor-1 cells, never used for phasing, only for checking |

- **Read preparation (per cell):**
  - remove duplicate molecules;
  - mask letters made ambiguous by bisulfite (e.g. T at a C/T SNP on a C→T-converted strand);
  - re-join the separately aligned pieces of each DNA molecule.

  | | n100 | n2370 |
  |---|---|---|
  | molecules (duplicates) | 220.3 M (26.7%) | 4.14 G (29.0%) |
  | molecules covering a readable het SNP | 12.0 M | 223.1 M |
  | molecules linking SNPs ≥ 1 kb apart on one chromosome | 89,585 | 1,435,492 |
  | SNP letters read: ref / alt / other | 7.73 M / 6.47 M / 0.12 M | 145.2 M / 122.6 M / 2.09 M (0.8%) |
  | letters masked for bisulfite | 5.69 M | 108.8 M |

## 2. SNP phasing of donor 1
- **SHAPEIT5** phases **1,890,745 of 1,927,881 het SNPs (98.1%)**. The rest are not in the panel.

| | n4 | n100 | n2370 |
|---|---|---|---|
| SNPs linked by RefLinker (scaffold), % of SHAPEIT5 SNPs | 7.2% | 85.1% | **97.6%** |
| SNPs phased by HapCUT2 (`--hic`) | 51,564 | 759,714 | 1,589,258 |
| SNPs kept after pruning (switch energy < 2 removed) | – | **38.8%** | **84.3%** |
| SNPs kept if the cut-off were 0 | – | 87.1% | 97.5% |

- RefLinker flips ~half of the SHAPEIT5 blocks (47% of scaffold SNPs). That is expected: SHAPEIT5's long-range
  orientation is random.
- **Pruning cut-off:** the cut-off of 2 (chosen on NA12878) prunes far more on sparse data, because switch energies grow
  with the amount of Hi-C. **From 2026-11-06 the cut-off is 0.**

## 3. Truth-free check on 20 held-out cells
Most contacts join two places on the **same** chromosome copy, so a correct phasing gives both ends the same label. A
random orientation gives ~50%. Only contacts with both ends on phased SNPs count.

**(a) On the authors' contact tables** (5.41 M same-chromosome contacts ≥ 1 kb; ~12.5 k with both ends labelled):

| Phasing | 1–10 kb | 10–100 kb | 0.1–1 Mb | 1–10 Mb | > 10 Mb | all |
|---|---|---|---|---|---|---|
| SHAPEIT5 alone | 99.1% | 96.5% | 79.7% | 57.7% | 49.1% | 79.6% |
| n100 | 98.7% | 96.5% | 93.6% | 92.3% | 85.7% | 94.1% |
| n100 + pruning | 99.2% | 98.2% | 94.8% | 91.8% | 84.7% | 95.5% (2.5 k) |
| n2370 | 97.3% | 96.1% | 95.5% | 94.5% | 87.2% | 94.8% |
| **n2370 + pruning** | 97.7% | 96.6% | **96.1%** | **95.0%** | **88.4%** | 95.5% (10.4 k) |

**(b) On contacts rebuilt from the reads** (same cells; first version of the check), which also covers n4 and the
confidence-guided (unpruned) variants:

| Phasing | 1–10 kb | 10–100 kb | 0.1–1 Mb | 1–10 Mb | > 10 Mb |
|---|---|---|---|---|---|
| SHAPEIT5 alone | 98.3% | 96.4% | 79.3% | 57.4% | 49.0% |
| n4 | 98.3% | 96.1% | 78.5% | 56.5% | 50.5% |
| n100 | 97.7% | 96.1% | 93.4% | 92.3% | 86.1% |
| n100 + pruning | 97.0% | 98.0% | 94.6% | 92.7% | 86.1% |
| n2370 | 95.8% | 95.6% | 95.3% | 94.8% | 86.9% |
| n2370, confidence-guided | 96.4% | 95.9% | 95.8% | 95.1% | 88.0% |
| n2370 + pruning | 96.4% | 96.0% | 95.9% | 95.1% | 88.1% |

- (a) and (b) agree within ~1 point in every bin.
- **4 cells don't help**: RefLinker links only 7% of SNPs.
- **100 cells restore long-range phase. 2,370 cells add 2–3 points at 0.1–10 Mb and 1–2 points at > 10 Mb**, i.e.
  diminishing returns.
- **Pruning helps only with enough data:** at 100 cells it removes 61% of SNPs for no gain; at 2,370 cells it removes
  16% and gives the best result.
- **Short range dips slightly with 2,370 cells** (1–10 kb: ~98–99% → 96–97%). Not yet explained.
- **What the ~12% at > 10 Mb is:** an upper bound on phasing error. Part of it is real contacts between the two copies.

## 4. Contact labels for the 4 pilot cells (`pilot_v2/original.h5`, entries 0, 1, 2, 4)
- **How a label is made:** each contact end is matched to its own read piece in the cell's BAM (676 k / 676 k ends
  matched). The SNP letter on that piece, or on another piece of the same molecule within 1 kb, gives 0 or 1.
- chrX ends get 1 (maternal) and chrY ends 0 (paternal): the donor is male. Entry 3 (donor 3) is left at −1.

| File (in `data/contacts/`) | phasing | ends labelled | contacts ≥ 1 end | both ends | same label, both ends SNP-labelled |
|---|---|---|---|---|---|
| `original_phased.h5` | n100 | 6.9–7.0% | 11.7–11.8% | 2.0–2.2% | 98.4–98.7% |
| `original_phased_pruned_n100.h5` | n100 + pruning | 4.7% | 7.3–7.4% | 1.9–2.1% | 98.7–99.2% |
| `original_phased_n2474.h5` | n2370 | 6.9–7.0% | 11.7–11.8% | 2.0–2.2% | 97.1–98.3% |
| **`original_phased_pruned_n2474.h5`** | n2370 + pruning | 6.5–6.6% | 10.8–11.0% | 2.0–2.2% | 97.9–98.6% |

- **Ends labelled:** per cell, ~3.8% of ends from SNPs and ~3.1% from the chrX/Y rule.
- **Why so few:** het SNPs are ~1 per 1.5 kb, read pieces are 40–130 bp, and bisulfite masking removes some letters.
  More cells don't change this; they make the labels more correct.
- **Label meaning:** autosome labels are consistent along a chromosome only (parent of origin unknown, no parental DNA).

## 5. Checks and side experiments
- **Verification:**
  - An independent re-implementation of hybrid + combiner gave 0 mismatches.
  - 600 / 600 randomly chosen contact labels were re-derived from the raw reads.
  - The output files differ from the source only in the phase columns.
- **SNPsplit** (standard bisulfite allele-sorting tool) agrees **100%** with our per-read calls (409,860 read pieces in the
  4 pilot cells). Ours assigns 8% more, because it keeps unambiguous converted letters.
- **hickit imputation** (Dip-C method, filling unlabelled ends from neighbours) recovers only 6–13% of hidden labels, at
  72–82% accuracy. It is **not used**; too few ends are labelled to start with.
- **WGS genotypes can't replace bisulfite masking:** on a converted strand, an unmethylated C and a real T both read as T.

## 6. Ongoing and next
- **~5,000 cells per donor, for donors 1, 2 and 4** (pruning cut-off 0), then labelling two new contact files: x-tzahin's
  `mtg50.h5` (200 cells: 100 donor 1, 100 donor 2) and `pilot_v5/original.h5` (150 cells: 51 / 52 / 47).

  | Donor | SHAPEIT5 | cells downloaded | cells processed | next |
  |---|---|---|---|---|
  | H1930001 | done (98.1%) | 4,828 (2,370 reused + 2,458 new; 192 failed 3×) | 2,370 | process → RefLinker / HapCUT2 → combine → check |
  | H1930002 | done (98.0% of 1,930,628) | 5,020 / 5,020 | 0 | same |
  | H1930004 | done (97.9% of 1,907,685) | 5,006 / 5,020 | 3,296 | same |

  - Blocked since ~11-09 by the exhausted Slurm allocation ("AssocGrpCPUMinutesLimit").
  - About 1 day of computing remains once jobs run.
  - The drivers (tmux on login05) resume where they stopped. They will need a restart after the Anvil downtime
    (Mon–Wed).
- **Joe's suggestion: pool reads from all cells instead of sampling whole cells.** The idea is that a contact between the
  two copies in one cell may repeat several times in that cell, while sampling a few reads from many cells dilutes it.
  - Cheap tests planned on donor 1, with no new downloads:
    - do discordant contacts cluster within a cell?
    - equal-budget comparison: 1,185 whole cells vs 50% of the molecules of all 2,370 cells.
  - Pooling across all ~49 k cells would mean ~16 TB of downloads per donor.
- **Questions for the working group:**
  - How should we judge whether phasing is good enough when only ~12% of contacts carry a label (Longzhi)?
  - Which cell-selection criteria should we use (Joe)?
- **Also open:** why short-range agreement dips with more cells; a cut-off for pruning tuned to each data depth.

## 7. Practical notes
- **Compute:** pilot ~101 core-h; pruning + 2,370-cell runs ~283 core-h; RefLinker peak 14 GB at 2,370 cells.
- **Downloads are the bottleneck:**
  - NeMO gives ~30–45 MB/s in total, ~3–5 cells per minute.
  - ~7% of downloads fail per pass; 2–4% never succeed.
  - Raw BAMs are deleted after processing, except the pilot, held-out and contact-file cells.
- **Issues fixed along the way:**
  - SHAPEIT5 needed allele-count fields I had stripped from the inputs (~6 h lost).
  - Three driver bugs (merge stats step, restart overwriting the download status, a too-strict download guard).
