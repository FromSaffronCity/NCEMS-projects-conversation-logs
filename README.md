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
| [`NA12878-phaser-verification/phaser_verification_report.md`](NA12878-phaser-verification/phaser_verification_report.md) | **Start here.** Status report for the PIs (updated 2026-10-19): key findings, method and data, results incl. hybrid / tuning / filtered-read tests (HapCUT2 in progress), next steps |
| [`NA12878-phaser-verification/figures/`](NA12878-phaser-verification/figures/) | The report's figure: phasing accuracy vs number of LiMCA cells (`fig1_cells_curve.png`) |

The scripts, logs and data of this experiment live on Purdue Anvil
(`/anvil/projects/x-bio260019/x-aahmed7/phaser-verification/`), not in this repository.
