# Figures for the supervisor (2026-11-19)
- `fig1_pipeline_schematic.png`: the current phaser pipeline: each step, what it does, inputs and outputs.
- `fig2_dipc17_ablation.png`: NA12878 Dip-C 17 cells, even chromosomes (October runs). What each step adds:
  SHAPEIT5 → + RefLinker → + hybrid refiner → + HapCUT2 (NA12878 WGS reads) → + pruning (RefLinker switch energy < 2
  unless HapCUT2 quality ≥ 30). Shows Hamming error, switch error and the longest error-free segment.
- `fig3_limca_cells.png`: NA12878 LiMCA with 17 / 50 / 100 / 220 cells (RefLinker gets the Hi-C reads of those N cells; HapCUT2 always gets the same bulk NA12878 WGS reads). Full pipeline with the same settings at every
  cell count (hybrid window 101, switch-energy voting >= 0, HapCUT2 quality >= 30, pruning < 2). 2x2: Hamming error, switch error, longest error-free segment, phased SNPs retained after pruning.
