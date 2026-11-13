# conversation-logs

Conversation and experimental logs from my NCEMS internship work, collaborating with a Claude agent.
2026-07-12 → present.

Project goal: haplotype-resolved 3D genome structures of the human brain. This needs the heterozygous SNPs
of each donor phased into chromosome-length haplotypes. The experiments below test ways of getting there.

## Experiments

### 1. Read-backed phasing from single-cell bisulfite data (2026-07-12 → 2026-10-05, closed)
Read-backed haplotype phasing of heterozygous SNPs called from single-cell **bisulfite** sequencing of human brain
(snMC-seq and snM3C-seq), across three donors and several pseudo-bulk depths.

| file | what it is |
|---|---|
| [`snMC-seq-snM3C-seq-phasing-experiment/PHASING_EXPERIMENTS_REPORT.md`](snMC-seq-snM3C-seq-phasing-experiment/PHASING_EXPERIMENTS_REPORT.md) | **Start here.** Objective, methods, every result table, findings and limitations |
| [`snMC-seq-snM3C-seq-phasing-experiment/snMC-seq_SNP_phasing_experiment.md`](snMC-seq-snM3C-seq-phasing-experiment/snMC-seq_SNP_phasing_experiment.md) | Session-by-session log for the snMC-seq experiments |
| [`snMC-seq-snM3C-seq-phasing-experiment/snM3C-seq_SNP_phasing_experiment.md`](snMC-seq-snM3C-seq-phasing-experiment/snM3C-seq_SNP_phasing_experiment.md) | Session-by-session log for the snM3C-seq experiments |

The data, callsets, phasings and analysis outputs of this experiment live on the CyVerse data store, not in
this repository; `workspace/data/README.md` there documents the layout.

### 2. Phaser verification on NA12878 (2026-10-07 → ongoing)
Statistical phasing (SHAPEIT5) followed by Hi-C-based refinement (RefLinker), validated against NA12878's
pedigree-phased truth with two independent single-cell Hi-C datasets (Dip-C, 17 cells; LiMCA, up to 220 cells).

| file | what it is |
|---|---|
| [`NA12878-phaser-verification/phaser_verification_report.md`](NA12878-phaser-verification/phaser_verification_report.md) | **Start here.** Status report for the PIs (updated 2026-10-23): summary, key findings, how results are scored, pipeline and data, results (Hi-C baseline; hybrid / tuning / filtered-read / HapCUT2 tests; keeping unlinked SNPs; confidence-score refinement and pruning), practical notes, next steps |
| [`NA12878-phaser-verification/figures/`](NA12878-phaser-verification/figures/) | The report's figures: `fig3_pipeline_steps_all17.png`, `fig3_pipeline_steps_limca_n220.png`, `fig1_cells_curve.png` |

The scripts, logs and data of this experiment live on Purdue Anvil
(`/anvil/projects/x-bio260019/x-aahmed7/phaser-verification/`), not in this repository.

### 3. Phaser pipeline guide (2026-11-10)
Step-by-step description of the phaser pipeline used for the human brain donors: SHAPEIT5 → RefLinker → hybrid refiner,
HapCUT2 in parallel, then combiner and pruner. For each step it gives the inputs, the command-line settings, the outputs,
and how our own steps work, with toy examples.

| file | what it is |
|---|---|
| [`phaser-pipeline-guide/phaser_pipeline_guide.md`](phaser-pipeline-guide/phaser_pipeline_guide.md) | **The guide** |
| [`phaser-pipeline-guide/figures/phaser_pipeline.png`](phaser-pipeline-guide/figures/phaser_pipeline.png) | Pipeline overview figure |

The code and data of the brain phasing live on Purdue Anvil (`/anvil/projects/x-bio260019/x-aahmed7/human-brain-phasing/`).

### 4. Human brain phasing (2026-10-22 → ongoing)
Phasing of the Tian et al. 2023 donors' WGS het SNPs with the pipeline above, using the donors' own snM3C-seq cells, and
labelling snM3C contacts with the chromosome copy for ChromoGen. Donor 1 is complete (4 / 100 / 2,370 cells, with and
without pruning); donors 1, 2 and 4 with ~5,000 cells are in progress.

| file | what it is |
|---|---|
| [`human-brain-phasing/brain_phasing_progress_report.md`](human-brain-phasing/brain_phasing_progress_report.md) | **Progress report** (2026-11-13): setup, all donor-1 numbers, checks, next steps |

Code, logs and data on Purdue Anvil: `/anvil/projects/x-bio260019/x-aahmed7/human-brain-phasing/`.
