# snMC-seq SNP Phasing Experiment — Consolidated Log

Consolidated summary of all prior Claude Code sessions (2026-07-12 → 2026-07-14) on
the Jingtian/Ecker human-brain snMC-seq single-cell SNP calling + phasing experiment.
Purpose: let a Claude agent catch up on everything done so far.

---

## Overall goal
Take one real single-nucleus **snmC-seq3** (bisulfite, PBAT) BAM from the Ecker human
brain methylome study, call **bisulfite-aware SNPs** with bsgenova, then **phase** the
resulting heterozygous SNPs (statistically with SHAPEIT5 and read-backed with HapCUT2).

## The dataset / cell
- Paper: Tian W, **Zhou J**, … Ecker JR, "Single-cell DNA methylation and 3D genome
  architecture in the human brain," *Science* 2023. DOI 10.1126/science.adf5357,
  PMID 37824674.
- GEO superSeries **GSE215353** (snmC-seq3, single nuclei, human).
- Cell used: `HBA_200130_H1930002_CX56_V2_1_P1-1-A18-A13` = **GSM5093977**
  (subseries GSE167119); donor H1930002.

---

## Session 1 (2026-07-12) — BAM download
- Located the data: GEO has only methylation calls (`.allc.tsv.gz`, no BAM); SRA/ENA
  have raw FASTQ only. The mapped **BAMs live at the NeMO Archive**
  (`nemo:dat-jx4eu3g` → Alignments `nemo:dat-budp3jg`, 524,004 per-cell
  `.final.bam.tar`, ~101 TB total).
- Downloaded `HBA_200130_...A13.final.bam` (94,834,923 B) + built `.bai`.
  - Source: `https://data.nemoarchive.org/biccn/grant/u01_ecker/ecker/epigenome/sncell/mCseq/human/processed/align/<cell>.final.bam.tar`
  - Integrity: tar md5 `9075a50ee859d075b9b0f6695fcc051a`; extracted BAM md5
    `ab51df927ddcbc65ef201442e06e9991`. `samtools quickcheck` PASS, 681,290 reads,
    coordinate-sorted, GRCh38.
- Authors' pipeline (from `@PG`): Bismark v0.22.3 (`--pbat` R1, normal R2, bowtie2) →
  `samtools view -q 10` → sort → Picard MarkDuplicates (REMOVE_DUPLICATES) →
  merge R1+R2 → `.final.bam`. **No re-mapping was done** — only downloaded + verified.
- Env note: `workspace/` is a FUSE iRODS (CyVerse) mount; large direct writes lag/corrupt,
  so always stage in local `/tmp`, verify, copy to mount, re-verify byte-for-byte.

## Session 2 (2026-07-12) — cleanup, alignment verify, bsgenova prep
- **Cleanup:** saved future-useful artifacts into `data/jingtian-snMCseq/`:
  `nemo_align_manifest.tsv` (524,004 rows, md5 `954359ea98c3006364cdc7f1233ad012`),
  its README, and `geo_gsm_to_cell.GSE167119.tsv` (3,071 rows). Removed /tmp intermediates.
- **Alignment confirmed:** ALIGNED (Bismark, bisulfite-aware, PBAT); 100% mapped,
  0 dups, 456 `@SQ` contigs = UCSC hg38 (chr1…chrY + randoms/alts) **+ `chrL`
  (48,502 bp) = lambda spike-in** (bisulfite conversion control). No re-alignment needed.
- **Reference built:** `data/reference/hg38_chrL.fa` (+`.fai`) = UCSC hg38.fa
  (md5 `1c9dcaddfa41027f17cd8f7a82c7293b`) + lambda NC_001416.1 renamed `>chrL`.
  Combined md5 `612ab3bbf05a3f916093e5cd11b31d5a`. All 456 contig names+lengths
  diff-identical to BAM `@SQ`.
- **Env `bsgenova`:** python 3.11.15, numpy 2.4.6, pysam 0.24.0 (numpy 2.x works despite
  pyproject pinning numpy==1.18). Smoke-tested bsgenova on bundled example + a real
  `chr1:100k-1M` pipe (exit 0; 0 SNVs in that sparse window = expected). bsgenova repo
  cloned at `workspace/bsgenova`.
- **Gotcha:** don't run two `conda run -n bsgenova` across a pipe (breaks stage 1);
  `conda activate` then call `python` directly.

## Session 3 (2026-07-14) — bsgenova SNP calling run
- **Ephemeral node:** conda envs from prior sessions were GONE (only `base` survived).
  Only `workspace/` mount data persists. Recreated env:
  `mamba create -n bsgenova -c conda-forge -c bioconda python=3.11 numpy pysam samtools`
  (→ py 3.11.15, numpy 2.4.6, pysam 0.24.0, samtools 1.24). Must
  `source /opt/conda/etc/profile.d/conda.sh` first (conda not on PATH in non-login shells).
- Staged reference (3.2 GB) + BAM to local `/tmp/bsg_run/` (iRODS reads are slow).
- **Run (SUCCESS, ~18 min, exit 0):** both stages in one `bash -c` pipe:
  ```
  python ./bsextractor.py -b /tmp/bsg_run/<cell>.final.bam -g /tmp/bsg_run/hg38_chrL.fa \
      --output-atcgmap - | python ./bsgenova.py -o /tmp/bsg_run/output/<cell> \
      --sample-name <cell> -P 16
  ```
  Default model params. Note: output files sit at 0 B until the very end (bsgenova
  writes in-order at completion) — normal, not a hang.
- **Output** → `data/jingtian-snMCseq/snv/` (md5-verified on mount):
  - `<cell>.vcf.gz` (md5 `c1444f9eb2df0a63b949ea5010a19aa5`)
  - `<cell>.snv.gz` (md5 `0c9cd690fab52ffcce5bd989ec7d6c96`)
- **11 SNVs** genome-wide (10 on chr16, 1 on chrY), all PASS, all het `0/1`. Sparse
  count expected for a single cell at ~1x.
- **Unphased confirmed** (user requirement): bsgenova is a per-site caller with no
  phasing; 11/11 GT use `/`, none `|`, no `##phasing` header.

## Session 4 (2026-07-14, session 2) — phasing: SHAPEIT5 vs HapCUT2
- Targets = the 10 chr16 het sites (chrY het excluded as haploid artifact).
- **Input prep:** bsgenova VCF lacked `##contig` headers → injected from reference `.fai`,
  bgzipped, tabix-indexed, split the one multiallelic site to biallelic, added AC/AN.
- **Envs (recreated, node ephemeral):**
  - `shapeit5`: shapeit5 5.1.1 (+ `libboost=1.85`; binaries are `SHAPEIT5_phase_common` etc.).
  - `bcftools`: bcftools **1.14** + htslib + `gsl=2.6` (bioconda's 1.23 has a broken gsl
    pin → solver hell; 1.14 + gsl 2.6 is the stable combo).
  - `hapcut2`: hapcut2 + samtools.
- **Approach 1 — SHAPEIT5: 0 sites phased.** Pulled only the two ±1 Mb windows around
  our loci from the remote 1000G 30x NYGC chr16 panel via tabix (14 MB, 68,659 biallelic
  SNPs, 3,202 samples). **None of the 10 SNPs are in the panel** (chr16:46,389,000–46,391,000,
  6 SNPs, is pericentromeric/segdup with zero panel variants). A bigger panel won't help —
  these are de-novo/rare single-cell calls, not catalogued common variants.
- **Approach 2 — HapCUT2: 2 blocks, 5 SNPs phased.** `extractHAIRS --bam --VCF` → 7
  fragments; `HAPCUT2` → 2 blocks:
  - Block 1 (46,389,911–46,389,974, 63 bp): copy1 = G,C / copy2 = A,G.
  - Block 2 (46,390,591–46,390,790, 200 bp): copy1 = T,A,G (all ALT) / copy2 = A,T,C (all REF).
  - Unphased/isolated: 34,592,801 / 46,390,374 / 46,394,853 / 46,399,426; 46,390,811 pruned
    (→ ./.); chrY.
  - **Caveat:** extractHAIRS ran in standard (non-bisulfite) mode; bisulfite C→T/G→A can
    bias C/T & C/G sites. Block 2 backbone (A/T sites) is bisulfite-neutral → solid.
- **Outputs** → `data/jingtian-snMCseq/phasing/`: `README.md` (full summary + haplotypes),
  `hapcut2/` (blocks.txt, phased.vcf.gz+tbi, extractHAIRS fragments), `shapeit5/`
  (small 1000G chr16 window panel bcf+csi).

---

## Current status
- BAM downloaded + verified; confirmed Bismark bisulfite-aligned. ✅
- Reference `hg38_chrL.fa` built + contig-matched. ✅
- bsgenova SNP calling done → 11 unphased het SNVs. ✅
- Phasing done: **HapCUT2 is the usable result** (2 phased blocks, 5 SNPs); SHAPEIT5
  not applicable to these de-novo/rare single-cell calls. ✅

## Key persistent facts / gotchas for next agent
- **Compute node is EPHEMERAL** — conda envs do NOT survive; only the `workspace/`
  iRODS mount persists. Recreate envs each session (recipes above); `source
  /opt/conda/etc/profile.d/conda.sh` before `conda activate`.
- iRODS FUSE mount is slow/lossy for large writes → always stage in `/tmp`, md5-verify
  both ways.
- Data locations under `workspace/data/jingtian-snMCseq/`: BAM+bai, `PROVENANCE.md`,
  manifests, `snv/` (bsgenova calls), `phasing/` (results); reference at
  `workspace/data/reference/hg38_chrL.fa`; bsgenova code at `workspace/bsgenova`.

## Possible next steps
- Re-run extractHAIRS with a bisulfite-aware mode/tool for higher confidence at C/T & C/G
  sites, or accept HapCUT2 blocks as final.

---

## Appendix — summaries of README markdown files from `workspace/data/`
Summaries of the three markdown files that lived under `workspace/data/jingtian-snMCseq/`
(consolidated here, then the originals removed from `data/`).

### `PROVENANCE.md` (BAM provenance)
- Describes one snmC-seq3 mapped BAM from the Tian/Zhou/Ecker human-brain methylome
  paper (*Science* 2023, DOI 10.1126/science.adf5357, PMID 37824674).
- Files: `HBA_200130_H1930002_CX56_V2_1_P1-1-A18-A13.final.bam` (94,834,923 B) + local `.bai`.
- Cell = **GSM5093977** (GSE167119, superSeries GSE215353); donor H1930002, region CX56,
  plate P1-1, well A18-A13; *Homo sapiens*, GRCh38/hg38; NovaSeq 6000.
- **Why no BAM on GEO/SRA:** GEO hosts only ALLC methylation calls, SRA only raw FASTQ;
  the mapped BAMs live at the **NeMO Archive** (`nemo:dat-jx4eu3g` → Alignments
  `nemo:dat-budp3jg`), downloaded as a per-cell tarball and extracted.
- Integrity: tar md5 `9075a50ee859d075b9b0f6695fcc051a`; extracted BAM md5
  `ab51df927ddcbc65ef201442e06e9991`.
- Header confirms **Bismark v0.22.3 `--pbat`** (PBAT = snmC-seq) pipeline → `samtools view
  -q 10` → sort → Picard MarkDuplicates → merge R1+R2 → `.final.bam`; ref path names first
  author Wei Tian's account. flagstat: 681,290 reads, 100% mapped, 0 dups.
- Note: `workspace/` is a FUSE iRODS (CyVerse) mount — large direct writes can corrupt, so
  the file was staged in `/tmp`, md5-verified, copied in, and re-verified byte-for-byte.

### `phasing/README.md` (SHAPEIT5 vs HapCUT2 phasing)
- Input: 11 unphased het SNVs from bsgenova (10 chr16, 1 chrY); chrY excluded as haploid →
  10 chr16 targets. Date 2026-07-14.
- **Approach 1 — SHAPEIT5 (statistical): 0 phased.** None of the 10 SNPs exist in the
  1000G 30x NYGC chr16 panel (dense locus chr16:46,389,000–46,391,000 is
  pericentromeric/segdup with zero panel variants). Statistical phasing only handles
  population-catalogued variants, not de-novo/rare single-cell calls.
- **Approach 2 — HapCUT2 (read-backed): 2 blocks, 5 SNPs phased.**
  - Block 1 (chr16:46,389,911–46,389,974, 63 bp): copy1 = G,C / copy2 = A,G.
  - Block 2 (chr16:46,390,591–46,390,790, 200 bp): copy1 = T,A,G (ALT) / copy2 = A,T,C (REF).
  - Unphased/isolated: 34,592,801 / 46,390,374 / 46,394,853 / 46,399,426; 46,390,811 pruned; chrY.
- Caveat: extractHAIRS ran in standard (non-bisulfite) mode; Block 2's A/T backbone is
  bisulfite-neutral → solid, Block 1 (single-fragment link) tentative.
- Includes an env recipe (shapeit5 5.1.1 + libboost 1.85, bcftools 1.14 + gsl 2.6, hapcut2).

### `phasing/bisulfite-aware-phasing-attempt.md` (bisulfite-aware re-run investigation)
- Goal: re-run HapCUT2 `extractHAIRS` in bisulfite-aware mode. Date 2026-07-15.
- **Finding: extractHAIRS has NO bisulfite mode** — verified three ways (binary rejects
  `--bisulfite`; help lists no such flag; upstream source has zero `bisulfite`/`methyl`
  matches). The earlier README suggestion was mistaken.
- Payoff would be limited anyway: only 3 of 9 chr16 het sites are bisulfite-sensitive
  (46399426, 46389911, 46394853); Block 2 is bisulfite-neutral and unaffected.
- **Decision (user, 2026-07-15):** keep existing HapCUT2 outputs as-is; Block 1 treated as
  tentative, Block 2 solid. Output files unchanged.
- If revisited: (1) emulate bisulfite-awareness by classifying Bismark `XG:Z:` CT/GA
  conversion strand and masking ambiguous alleles before standard extractHAIRS;
  (2) cross-check with WhatsHap (validation only); (3) accept standard blocks as final.

---

## Session 5 (2026-07-21) — FRESH RESTART: wiped `data/`, downloaded 10-cell CX56 BAM set
**Drastic reset requested by user** — start the experiment over from scratch.

### Cleanup
- **Deleted everything inside `workspace/data/`** (~5.5 GB: the prior single BAM+bai,
  `reference/hg38_chrL.fa`, `jingtian-snMCseq/` snv+phasing+manifests, 1000G chr16/chr1
  panels). Temp staging (`/tmp/bsg_run`, `/tmp/hapcut2-src`) already gone — nothing to clean.
- Created new subfolder **`workspace/data/Science-snMC-seq/`** for the fresh dataset.
- Note: prior consolidated experiment (Sessions 1–4 above) is now removed from disk; this
  log is the only record of it.

### Paper / GEO / NeMO verification (web research)
- **Yes, snMC-seq3 & from the Science paper.** GEO **GSE215353** is single-nucleus
  methylC-seq (snmC-seq3); it cites **PMID 37824674** = Tian W, Zhou J, … Ecker JR,
  *Science* 2023, DOI 10.1126/science.adf5357. Same dataset, same paper.
- **GEO/SRA host no BAMs** — only raw **FASTQ** (SRA) and processed **ALLC**
  (`*_allc.tsv.gz`) methylation-call tables. Aligned **BAMs live only at the NeMO Archive**
  (`nemo:dat-jx4eu3g` → alignments `dat-budp3jg`), per-cell `<cell>.final.bam.tar`.
- Paper scope confirmed verbatim: **46 brain dissection regions** (CX 22, BF 2, BN 11,
  HIP 5, THM 2, MB 1, PN 1, CB 2) and **188 cell types**; 517k cells from 3 adult male brains.
- Aligner = **Bismark** (bisulfite, `--pbat`); reference build = **`hg38-donor-snp-substituted`**
  (a common hg38 base with each donor's SNPs substituted in, 3 donors).

### Download — 10 BAMs, region CX56
- Selected the **10 smallest BAM-tars in region CX56** (region code = 4th `_`-field of cell
  name) from 2 donors (H1930001 ×6, H1930002 ×4), ~333 MB of tar total.
- Per cell: `curl` tar from
  `https://data.nemoarchive.org/biccn/grant/u01_ecker/ecker/epigenome/sncell/mCseq/human/processed/align/<cell>.final.bam.tar`
  → md5-verify tar vs manifest → `tar -xf` (extracts to `<cell>.final/<cell>.final.bam`)
  → copy `.bam` to `Science-snMC-seq/`. **All 10 md5-OK, all pass `samtools quickcheck`.**
  One (`HBA_200124_H1930002_CX56_V1C_1_P5-2-P3-H4`) re-extracted independently → mount copy
  **byte-identical** (md5 `a154f8a1317517dec5d61d3711d16f77`).
- **Gotcha (FUSE):** immediate post-`cp` `cmp`/`md5sum` on the iRODS mount gives false
  mismatches (read-after-write lag); re-reading after the write settles shows correct md5.
  So verify large-file copies on a short delay, not instantly.

### Reference-genome finding (for downstream tasks)
- All 10 BAMs share an **identical 456-contig `@SQ`** (same fingerprint
  `fa19c2b842380540589a2082c9de98b4` = same hg38 coordinate system, incl. **`chrL` lambda**
  spike-in), all mapped with **Bismark `--pbat`**. H1930002 cells' `@PG` names ref path
  `/gale/netapp/home/wtian/refs//hba-donor/h1930002`.
- **Same coordinates across both donors; base sequence differs only at donor-specific SNP
  sites** (hg38-donor-snp-substituted).
- **User Q — one reference enough despite 2 donors?** → For coordinate/methylation/coverage
  work, **yes: a single standard hg38** (matching the 456-contig `@SQ`, ideally +`chrL`) is
  fine. Only ref-base-sensitive steps (bisulfite-aware SNP calling, ref-vs-read allele
  comparison) might want a per-donor SNP-substituted reference.

### Current contents of `workspace/data/Science-snMC-seq/`
- 10 × `HBA_*_CX56_*.final.bam` (~34–35 MB each; no `.bai` built yet).
- `nemo_align_manifest.tsv` (524,004 rows, md5 `954359ea98c3006364cdc7f1233ad012`) +
  `nemo_align_manifest.README` — kept here for future NeMO downloads.

### Next steps (fresh experiment)
- Obtain/build the reference genome (standard hg38, add `chrL` lambda to match `@SQ`;
  consider per-donor SNP substitution only if a ref-base-sensitive step needs it).
- Build `.bai` indexes for the 10 BAMs; then proceed to SNP calling / phasing as before,
  now across a 10-cell CX56 set rather than a single cell.

## Session 6 (2026-07-22) — 3-donor × 10-cell CX45 set + per-donor merged bulk BAMs
**User goal:** 3 donor subfolders, 10 same-"cell-type" snMC-seq BAMs each, then
`samtools merge` → one bulk BAM per donor.

### Two conceptual corrections made up front (verified against manifest + its README)
- **"CX56" is a brain REGION, not a cell type.** Manifest README: *"brain-region code =
  4th underscore field."* Cell name encodes donor + region + plate/well only; the paper's
  **188 cell types** are downstream methylation-cluster annotations, NOT in the filename or
  manifest. Knowing a cell's true type would need the paper's per-cell annotation table.
  → **User chose to group by brain REGION.** (For germline SNP/phasing, cell-type
  homogeneity is irrelevant anyway — germline variants are identical across all cells of a
  donor.)
- **The 3 donors are H1930001, H1930002, H1930004** (NOT H1930003). **CX56 exists only for
  H1930001 & H1930002**; H1930004 has no CX56. Regions common to all 3 donors: CX44–CX52,
  CX55, BS91. → **User chose region CX45** for all three.
- **Manifest mixes two assays:** snMC-seq (bare cell names → `.final.bam.tar`) AND
  **snm3C-seq** (names contain `_3C_`, listed as `*.3C.sorted.bam.tar`). Also stray
  `M1C_3C_00{1,2}_*` (mouse 3C) + a `Methylome` comment row. **Must filter OUT `3C`** to get
  pure snMC-seq. (A naive 250 MB size target initially pulled mostly 3C files — caught & fixed.)

### What was done
- Created `Science-snMC-seq/{H1930001,H1930002,H1930004,merged-BAM}/`.
- **Selection:** the **10 LARGEST snMC-seq (non-3C) CX45 tars per donor** (max coverage for
  pseudo-bulk). snMC-seq CX45 pool sizes: H1930001 n=11797, H1930002 n=5839, H1930004 n=11657;
  tar sizes ~180–600 MB for the chosen cells.
- **Downloaded 30 BAMs** (curl tar → md5-verify vs manifest → `tar -xf` → cp to donor folder).
  **All 30: md5=YES, quickcheck=YES.** (script `/tmp/dl_cx45.sh`)
- **`samtools merge` (v1.9) per donor** → `merged-BAM/<donor>_CX45_snMCseq_10cells_merged.bam`
  (+`.bai`). Precheck confirmed 10 CX45-only per donor before merging; read counts cross-checked:
  - H1930001: **30,441,510** reads (3.9 GB), sum-of-inputs == merged ✅
  - H1930002: **19,211,865** reads (2.5 GB) ✅
  - H1930004: **37,517,754** reads (4.6 GB) ✅
  - all quickcheck OK, all 456-contig hg38+chrL `@SQ`, 100% mapped, 0 dups.
  - **`samtools merge` is the correct tool**: concatenates all reads into one coord-sorted
    BAM, merges headers, no cross-cell dedup (want every read for pseudo-bulk); all inputs
    share identical `@SQ` so coordinate-consistent.
- Merged-BAM filename attributes = donor + region (CX45) + assay (snMCseq) + #cells (10).
  Note field-5 sub-region varies within a donor (CaB/NAC/A25/A38/A24/Pir/A44_A45) — not a
  shared attribute, so not in the name.

### FUSE gotcha (again, bit us this session)
- Early `mv HBA_*_CX56_*.final.bam prior-CX56-unused/` reported "10 moved" off a **lagged
  listing**; the CX56 files actually landed back in the donor folders (sorted by donor ID)
  and `prior-CX56-unused/` was empty. Merge was unaffected (precheck saw only the 10 CX45).
  Cleaned by `find . -name '*_CX56_*.final.bam' -delete`. **Lesson: after any mv/cp/rm on the
  mount, re-verify with `find`, not the immediate `ls`.**

### Current contents of `Science-snMC-seq/`
- `H1930001/`, `H1930002/`, `H1930004/` — 10 CX45 snMC-seq `.final.bam` each (30 total, no `.bai`).
- `merged-BAM/` — 3 bulk BAMs + `.bai` (names above).
- `nemo_align_manifest.tsv` + `.README`. (Prior 10 CX56 BAMs deleted — reproducible from manifest.)
- A `hg38-reference/` folder is present (not created by this task).

---

## Session 6 (2026-07-22) — dataset re-layout, `bsgenova` env rebuilt, reference-genome decision
The `Science-snMC-seq/` folder has been **reorganised since Session 5**: no longer a flat
10-cell CX56 set. It is now **CX45**, split by donor, plus staging folders.

### Current contents of `workspace/data/Science-snMC-seq/` (after this session's reorg)
- `H1930001/` — **16 BAMs** (10 CX45 [NAC+CaB, ~350–570 MB] + 6 CX56 [~34–35 MB]).
- `H1930002/` — **14 BAMs** (10 CX45 [A38+A25, ~194–483 MB] + 4 CX56).
- `H1930004/` — **10 BAMs** (10 CX45 [A24 + A44_A45 + Pir]; grew from 7→10 mid-session as the
  user was still populating it; **no CX56** for this donor).
- `merged-BAM/` — empty (staging for a future merge → the intended "bulk MC-seq BAM").
- `hg38-reference/` — **NEW: reference genome for bsgenova** (see below).
- `nemo_align_manifest.tsv` (+`.README`) — kept for future NeMO pulls.
- `prior-CX56-unused/` — **REMOVED this session**; its 10 CX56 BAMs were moved into the donor
  folders above (donor = 3rd `_`-field of the cell name), sizes verified byte-for-byte post-move.
- No `.bai` built yet.

### hg38 reference — BUILT + placed on mount ✅ (`Science-snMC-seq/hg38-reference/`)
- `hg38_chrL.fa` (3,273,530,352 B, md5 `612ab3bbf05a3f916093e5cd11b31d5a`) + `.fai`
  (md5 `ec67be0edabb601e9e6484cf0c2aa4a5`) + `README.md` (full provenance).
- = **UCSC hg38** (`bigZips/hg38.fa.gz`, 455 contigs) **+ lambda `NC_001416.1`→`chrL`** = 456 contigs.
- md5 is **byte-identical to the S2 reference** (`612ab3bbf05a3f916093e5cd11b31d5a`) and its 456
  contig names+lengths **diff-match the current CX45 `@SQ` exactly** (verified this session).
- Built in `/tmp/hg38_build`, staged→mount, md5-verified (`.fai` matched immediately; `.fa` hit a
  transient FUSE "No such file" read-after-write glitch then verified on re-read — the known lag).

### `bsgenova` env — RECREATED (node ephemeral, prior envs gone; only `samtools` env survived)
Recipe (identical outcome to Sessions 2/3 — numpy 2.x works despite pyproject pinning `numpy==1.18`):
```bash
source /opt/conda/etc/profile.d/conda.sh          # conda not on PATH in non-login shells
mamba create -y -n bsgenova -c conda-forge -c bioconda python=3.11 numpy pysam samtools
conda activate bsgenova
```
Resolved versions: **python 3.11.15, numpy 2.4.6, pysam 0.24.0, samtools 1.24.**
Verified: `bsgenova.py` on the bundled `data/example.ATCGmap.gz` → valid `.snv.gz`+`.vcf.gz`
(exit 0); `bsextractor.py --help` loads (pysam OK). bsgenova code now lives at
`workspace/tools/bsgenova` (was `workspace/bsgenova` in earlier sessions).

### `samtools` env — separate standalone env (survived this session; here for reproducibility)
A dedicated `samtools` env exists alongside `bsgenova` (used for `@SQ`/`@PG` header inspection,
`quickcheck`, `.bai` indexing). It is **samtools 1.9 / htslib 1.9** (bioconda). Recipe to recreate:
```bash
source /opt/conda/etc/profile.d/conda.sh
mamba create -y -n samtools -c bioconda -c conda-forge samtools=1.9
conda activate samtools
```
(Note: the `bsgenova` env already bundles samtools **1.24**; this standalone 1.9 env is just a
lighter, independent tool env. Either works for BAM inspection/indexing.)

### Reference-genome decision (KEY — answers the recurring donor-reference question)
Inspected `@SQ` + `@PG` of one BAM per donor:
- **Coordinates identical across all 3 donors** — SN/LN-only fingerprint
  `4a4a814d679ecd4cd8d312b0db8d157f` (456 contigs incl. `chrL` lambda 48,502 bp), and it is
  **byte-identical to the old CX56 set**. (The full-`@SQ` md5 `7766f5c0…` differs from S5's
  `fa19c2b8…` only because of non-SN/LN tags/ordering — the real coordinate system is the same.)
- **`@PG` shows each donor aligned to its OWN reference:** Bismark `--pbat` against
  `/home1/07771/wtian/work/refs/hba-donor/h1930001` (resp. `h1930002`, `h1930004`) — i.e. the
  paper's per-donor **SNP-substituted (personalised) hg38**. (Path differs from S5's
  `/gale/netapp/...` note = same author, different cluster; structure `hba-donor/hXXXXXXX` is consistent.)

**Verdict: use the ORIGINAL standard hg38 (+`chrL` to match `@SQ`) — NOT the per-donor refs.** Why:
1. *Coordinate-compatible*: SNP substitution changes bases, never contig names/lengths, so one
   standard hg38+chrL matches all 3 donors' `@SQ`. (Rebuild the S2 `hg38_chrL.fa`; it was wiped in S5.)
2. *Scientifically correct for SNP calling*: bsgenova calls a site as variant when reads differ
   from the reference base. A **donor-substituted reference has the donor's own alleles baked in**,
   so at germline SNP sites reads would *agree* with the reference → those SNPs would be
   **suppressed** (you'd only see deviations from the donor's genome, i.e. somatic/mosaic). Calling
   against standard hg38 is what surfaces the donor's germline variants — the thing we want.
3. *Bonus, no downside*: because reads were aligned to a **personalised** reference, ref-mapping
   bias at het sites is already reduced → allele counts are *cleaner* than a standard-hg38 alignment
   would give. bsextractor only tallies read bases at recorded coordinates (no realignment), so
   feeding it standard hg38 is consistent.
- We **do not have and do not need** the authors' per-donor refs (they live on TACC
  `/home1/07771/wtian/...`, inaccessible). They'd only matter for calling *somatic/mosaic* variants
  relative to each donor's germline — not the current goal.

### Next steps
- Rebuild `data/reference/hg38_chrL.fa` (+`.fai`) = UCSC hg38 + lambda `NC_001416.1`→`chrL`
  (S2 combined md5 was `612ab3bbf05a3f916093e5cd11b31d5a`); diff its contig names/lengths against
  a current CX45 `@SQ` before use.
- Build `.bai` for the CX45 BAMs; run bsgenova per cell (or per merged-BAM) against `hg38_chrL.fa`.

---

## Session 7 (2026-07-22/23) — dataset re-layout + bsgenova run on 3 pseudo-bulk BAMs

### Directory re-layout (new top-level experiment folder)
Created `workspace/data/snMC-seq_SNP_phasing_experiment/` and moved everything under it:
- `hg38-reference/` — moved OUT of `Science-snMC-seq/` up to the experiment root.
- `Science-snMC-seq/` — moved wholesale into the experiment folder (now WITHOUT `hg38-reference/`;
  still holds `H1930001/ H1930002/ H1930004/ merged-BAM/` + `nemo_align_manifest.tsv`/`.README`).
- `bisulfite-aware-SNPs/from-bsgenova/` and `bisulfite-aware-SNPs/from-naive/` — new empty subdirs
  (the `from-naive/` one reserved for a later naive/non-bisulfite SNP-calling comparison).
- Final layout under `snMC-seq_SNP_phasing_experiment/`: `hg38-reference/`, `Science-snMC-seq/`,
  `bisulfite-aware-SNPs/{from-bsgenova,from-naive}/`.

### bsgenova SNP calling on the 3 per-donor pseudo-bulk merged BAMs (SUCCESS, all exit 0)
- **Inputs** = the 3 bulk BAMs in `Science-snMC-seq/merged-BAM/` (each = `samtools merge` of 10
  CX45 cells, from S6): `H1930001_…` (3.88 GB), `H1930002_…` (2.55 GB), `H1930004_…` (4.62 GB),
  all with `.bai`.
- **Env:** `bsgenova` conda env SURVIVED this session (also `samtools`; only these two exist).
  Verified: python 3.11.15, numpy 2.4.6, pysam 0.24.0, bsgenova.py 1.0.0. Tools at
  `workspace/tools/bsgenova`. Node had 256 cores, /tmp on 6.9 TB overlay (ample).
- **Reference:** standard `hg38-reference/hg38_chrL.fa` (the S6 build, md5
  `612ab3bbf05a3f916093e5cd11b31d5a`) — per the S6 decision (call germline variants vs standard
  hg38, NOT the per-donor SNP-substituted refs). Default bsgenova model params.
- **Procedure (FUSE-safe):** staged ref + all 3 BAMs(+`.bai`) to `/tmp/bsg_run/` (background cp,
  ~7 min — iRODS reads slow, foreground `cp` timed out at 2 min). Ran all 3 pipelines in PARALLEL
  (`conda activate bsgenova`, then per donor:
  `python ./bsextractor.py -b <bam> -g hg38_chrL.fa --output-atcgmap - | python ./bsgenova.py
  -o /tmp/bsg_run/output/<name> --sample-name <name> -P 16`). Wall time ~30–45 min per BAM
  (bsextractor single-threaded scan over 20–37 M reads is the bottleneck; outputs sit at 0 B
  until the very end — normal). All three printed `EXIT_0`.
- **Output naming:** bsgenova `-o <prefix>` emits `<prefix>.snv.gz` + `<prefix>.vcf.gz`; used
  `<bam_basename_minus_.bam>` as prefix → filenames match the requested
  `<input_filename_minus_ext>.<ext>` format.
- **Results (all PASS, all het GT `0/1`, ZERO phased `|` — unphased as required):**
  - H1930001: **1,166** SNVs total, **863** het.
  - H1930002: **421** SNVs total, **289** het.
  - H1930004: **1,712** SNVs total, **1,366** het.
  - (Far more than the single-cell S3 run's 11 — expected, these are 10-cell pseudo-bulk.)
- **Outputs** → `bisulfite-aware-SNPs/from-bsgenova/` (6 files, gzip-tested OK, md5 byte-identical
  local↔mount after a short settle):
  - `H1930001_CX45_snMCseq_10cells_merged.vcf.gz` (md5 `ce0fe7842a23a84ec563ae6a8790329e`),
    `.snv.gz` (`9334b3d8b10d15abe17ca7dfc25d477e`)
  - `H1930002_CX45_snMCseq_10cells_merged.vcf.gz` (`8aa364b71ff9cf20badaf583fc769f77`),
    `.snv.gz` (`b31af01ce8cb7ca0c7c3b4c2894f3a07`)
  - `H1930004_CX45_snMCseq_10cells_merged.vcf.gz` (`c56e8a528e5a8f9e8c0000ea21e6b55f`),
    `.snv.gz` (`5e89f91d78e99feee2fd35348a9c1539`)
  - `.vcf.gz` = phasing-ready VCF (`GT:GQ:GQH:DP:DPW:DPC`); `.snv.gz` = fuller per-site table
    (Watson/Crick strand base counts + posterior stats).
- **Cleanup:** `/tmp/bsg_run` removed (confirmed absent).

### Next steps
- Phase the bsgenova het SNPs (read-backed HapCUT2 as in S4; SHAPEIT5 likely N/A again for
  rare/de-novo calls). Now 3 richer per-donor callsets (289–1,366 het) instead of 1 cell.
- (Optional) populate `bisulfite-aware-SNPs/from-naive/` with a naive (non-bisulfite-aware)
  caller run for comparison against the bsgenova bisulfite-aware calls.

---

## Session 8 (2026-07-28) — naive bisulfite-aware SNP caller written; env recipes added; naive run launched

### Progress since last session (before this run)
- **Wrote a naive bisulfite-aware SNP caller** at
  `workspace/Python-scripts/bisulfite_aware_naive_SNP_caller.py` (v1.0.0). It calls SNPs
  directly from a Bismark-aligned bisulfite BAM + reference FASTA and emits **bsgenova-compatible
  output** (`<prefix>.snv.gz` + `<prefix>.vcf.gz`, same columns / VCF `GT:GQ:GQH:DP:DPW:DPC`) so
  the two callers are directly comparable (`from-naive/` vs `from-bsgenova/`).
  - **Bisulfite-awareness = hickit.js masking rule.** Based on `_hic_bs_skip_snp` in **hickit.js**
    (https://github.com/FromSaffronCity/hickit): on a forward/Watson read an observed `T` is
    ambiguous (real T or converted C) and on a reverse/Crick read an observed `A` is ambiguous
    (real A or converted G). The caller turns this into a strand-specific base mask — keep
    `n_A = W_A`, `n_T = C_T`, `n_C = W_C + C_C`, `n_G = W_G + C_G` — then runs a plain diploid
    genotype-likelihood model (flat error model, bsgenova-style prior; **no** methylation
    transition matrix). "Naive" = the mask is the only bisulfite ingredient. Consequence (expected):
    genomic C>T / G>A changes are indistinguishable from conversion → deliberately not called;
    all other substitution classes remain callable.
  - Read filtering / strand-split coverage (`count_coverage` per Watson/Crick) mirror
    `bsextractor.py`; genotype/allele-weight bookkeeping mirrors `bsgenova.py`.
- **Created environment recipes** in `workspace/environments/` so every conda env is rebuildable
  from scratch on the ephemeral node (only the `workspace/` iRODS mount persists):
  - `_common.sh` — shared helpers. **KEY LESSON (2026-07-28):** the node's base conda config lists
    only `defaults`; appending `-c bioconda -c conda-forge` alone makes mamba mix ABIs across
    channels and all four envs fail to solve (libdeflate/libnghttp2/c-ares conflicts). Fix =
    always use `--override-channels --strict-channel-priority -c conda-forge -c bioconda`.
  - `01_bsgenova.sh` (python 3.11 + numpy + pysam + samtools — also serves the naive caller),
    `02_htslib-tools.sh` (bcftools/htslib), `03_hapcut2.sh`, `04_shapeit5.sh`,
    `00_build_all.sh` (rebuild all four). `FORCE_RECREATE=1` to wipe+rebuild.

### This session's run (naive SNP caller on the 3 pseudo-bulk BAMs — IN PROGRESS)
- Cleanup done first (per user): wiped everything in `bisulfite-aware-SNPs/from-naive/` (a prior
  attempt had left 0-byte `.snv.gz`/`.vcf.gz` + a `logs/` folder — it had logged "2008 windows"
  then died, presumably from running all 3 BAMs at once / a disconnect). Also removed
  `.claude/scheduled_tasks.lock`.
- Recreated the `bsgenova` env via the recipe → python 3.11.15, numpy 2.4.6, pysam 0.24.0.
- **Strictly one-BAM-at-a-time** (user's explicit ask, to avoid memory overload from the prior
  parallel attempt): reference staged once to `/tmp`; then per BAM → stage BAM+`.bai` to `/tmp` →
  run the naive caller → copy `.snv.gz`/`.vcf.gz` to `from-naive/` (gzip-tested + md5-verified) →
  wipe the per-BAM `/tmp` staging → next BAM. Reference `hg38-reference/hg38_chrL.fa` (standard
  hg38+chrL, the S6 decision: call germline variants vs standard hg38). Results recorded below on
  completion.

### Still to do (both callers)
- Run **HapCUT2** (read-backed) and **SHAPEIT5** (statistical) on the unphased heterozygous SNPs
  in `bisulfite-aware-SNPs/` from **both** callers (`from-bsgenova/` and `from-naive/`), to phase
  and cross-compare the two callsets.

---

## Session 9 (2026-07-28) — naive SNP caller: 3 BAMs run CONCURRENTLY → `from-naive/` populated ✅
Completes the S8 naive-caller run. This time the user explicitly asked for a **concurrent**
(all-3-at-once) run, with resources capped so it would NOT disrupt the VSCode-server (the earlier
one-at-a-time framing is superseded).

### Node state (ephemeral, as always)
- Only `base` conda env survived → **rebuilt `bsgenova` via `environments/01_bsgenova.sh`**
  (python 3.11.15, numpy 2.4.6, pysam 0.24.0, samtools 1.24). The `_common.sh` channel fix
  (`--override-channels --strict-channel-priority -c conda-forge -c bioconda`) worked cleanly.
- Node had **128 cores, 476 GiB free RAM, 6.2 TB free `/tmp`** — ample.
- `from-naive/` was already empty (cleaned in S8); nothing to wipe.

### Resource-capped concurrent run (the key ask)
- Staged reference (3.27 GB) + all 3 merged BAMs(+`.bai`, 3.88/2.55/4.62 GB) to `/tmp/naive_run/`
  (background cp, ~9 min; iRODS reads slow). Sizes verified vs mount.
- Launched **all 3 callers at once**, detached with `setsid nohup` (survives disconnects — the
  suspected cause of the S8 death), each `nice -n 10` and **`-P 16`** workers:
  **48 of 128 cores + 3 mains**, leaving ~77 cores free. Peak RAM stayed tiny (459 GiB still free
  during the run). This is the "wise allocation" that kept the VSCode-server responsive.
  ```
  python bisulfite_aware_naive_SNP_caller.py -b <bam> -g hg38_chrL.fa \
      -o /tmp/naive_run/output/<name> --sample-name <name> -P 16
  ```
- 2008 genome windows each; all three **exit 0** (`ALL_FINISHED overall_rc=0`) in **~13–14 min**
  wall (concurrent): H1930002 790 s, H1930001 830 s, H1930004 850 s.

### Results (all PASS, ZERO phased `|` — unphased as required)
| donor | total SNVs | het (`0/1`+`1/2`) | hom (`0/0`+`1/1`) |
|-------|-----------:|------------------:|------------------:|
| H1930001 | 660 | 440 | 220 |
| H1930002 | 306 | 190 | 116 |
| H1930004 | 688 | 485 | 203 |

- **Fewer calls than bsgenova** (S7 het: 863 / 289 / 1366) — expected: the naive strand-mask
  discards conversion-ambiguous evidence, so it is a more conservative caller.

#### `from-bsgenova` vs `from-naive` — SNP-count difference (naive is the more conservative caller)
Het column for bsgenova is the S7-recorded het count; bsgenova total-SNV counts (S7): 1166 / 421 / 1712.
| donor | bsgenova total | naive total | Δ total (naive−bsgenova) | bsgenova het | naive het | Δ het |
|-------|---------------:|------------:|-------------------------:|-------------:|----------:|------:|
| H1930001 | 1166 | 660 | −506 | 863 | 440 | −423 |
| H1930002 |  421 | 306 | −115 | 289 | 190 |  −99 |
| H1930004 | 1712 | 688 | −1024 | 1366 | 485 | −881 |
| **sum**  | **3299** | **1654** | **−1645** | **2518** | **1115** | **−1403** |
- Across all 3 donors the naive caller yields **~half** the SNVs of bsgenova (1654 vs 3299 total;
  1115 vs 2518 het). Direction is expected: bsgenova models conversion probabilistically and keeps
  all read evidence, whereas the naive caller hard-masks conversion-ambiguous bases (Watson `T`,
  Crick `A`) and so discards evidence at C/T- and G/A-type sites → fewer, more conservative calls.
- **NOTE on C>T / G>A calls (contradicts the script docstring):** the caller *does* emit C>T /
  G>A variants (107 / 51 / 119). The docstring says these are "deliberately not called," but the
  actual mask keeps `n_T = C_T` (Crick-strand T) and `n_A = W_A` (Watson-strand A). A forward-ref
  `T` on a Crick read is not a conversion product (bottom-strand C→T shows as G→A in ref coords),
  and a Watson `A` likewise is not — so C>T/G>A stay callable from their *conversion-unambiguous*
  strand; only the ambiguous-strand evidence is dropped. The outputs are correct per the code as
  written; it is the docstring's blanket "not called" wording that overstates it. **Script left
  UNCHANGED** (not in scope this session) — flagged here for the `from-naive` vs `from-bsgenova`
  comparison.

### Outputs → `bisulfite-aware-SNPs/from-naive/` (gzip-tested OK; md5 byte-identical local↔mount)
- `H1930001_CX45_snMCseq_10cells_merged.vcf.gz` (md5 `0ce16e02f696b3a50640ca9431c91ee2`),
  `.snv.gz` (`5a855a8d46ae2285d52a5e36f2021491`)
- `H1930002_CX45_snMCseq_10cells_merged.vcf.gz` (`50c4ec0da8cf8ef9fe2ba38068537b2a`),
  `.snv.gz` (`678b8a8989ec4eda0812d1eba992fe26`)
- `H1930004_CX45_snMCseq_10cells_merged.vcf.gz` (`244dff3065a1d94cfb04aa301cde5648`),
  `.snv.gz` (`0b15d0cf352381c1e83291d7d34f65b9`)
- Same formats as bsgenova (`.vcf.gz` = `GT:GQ:GQH:DP:DPW:DPC`; `.snv.gz` = 13-col per-site table).
- **Cleanup:** `/tmp/naive_run` removed (confirmed absent).

### Next steps
- Both callsets now present (`from-bsgenova/` + `from-naive/`). Proceed to phasing: **HapCUT2**
  (read-backed) and **SHAPEIT5** (statistical, likely N/A again for rare/de-novo calls) on the
  het SNPs of both, then cross-compare bsgenova vs naive.
- (Optional) reconcile the naive caller's docstring with its actual C>T/G>A behavior.

---

## Session 10 (2026-07-28) — HapCUT2→SHAPEIT5 tooling set up (NO tool run yet)

### Naive caller results — CONFIRMED PRESENT AND VALID (see Session 9 above)
This session ran concurrently with the naive-caller session. A mid-session check
at 21:08 saw `from-naive/` empty and `/tmp/naive_run/` still staged, and I briefly
recorded that the results had not been copied — **that was a snapshot of work still
in flight, and it is superseded.** The copy landed at 21:14 and `/tmp/naive_run`
was cleaned. Independently re-verified here:

| donor | records | het (`0/1`) | phased (`\|`) | `.snv.gz` | `.vcf.gz` |
|---|---|---|---|---|---|
| H1930001 | 660 | 439 | 0 | gzip OK | gzip OK |
| H1930002 | 306 | 186 | 0 | gzip OK | gzip OK |
| H1930004 | 688 | 481 | 0 | gzip OK | gzip OK |

Record counts match the caller logs exactly (660/306/688), byte sizes match the
`/tmp` originals, and **all genotypes are unphased** — correct input for phasing.

**Concurrency lesson:** two agents were appending to this log at once (hence two
"Session 9" headings — this one renumbered to Session 10). Before reporting the
state of the mount, re-check immediately; a listing minutes old may describe work
another session has since finished.

### Environments — ALL FOUR REBUILT AND VERIFIED (node was ephemeral again; only `bsgenova` survived)
`bash environments/00_build_all.sh` → **ALL ENVS BUILT OK**. The
`--override-channels --strict-channel-priority` fix in `_common.sh` works; every
recipe solved first try.

| env | contents | smoke test |
|---|---|---|
| `bsgenova` | python 3.11.15, numpy 2.4.6, pysam 0.24.0, samtools 1.24 | pre-existing, OK |
| `htslib-tools` | bcftools 1.14, samtools 1.14, htslib 1.14, gsl 2.6 | `+fill-tags` present |
| `hapcut2` | hapcut2 1.3.4, samtools 1.23.1, python 3.11, pysam | both binaries print usage |
| `shapeit5` | shapeit5 5.1.1, libboost 1.85 | `phase_common --scaffold` present |

- **`SHAPEIT5_xcftools` is MISSING** from the bioconda 5.1.1 build (the other four
  `SHAPEIT5_*` binaries are all present). Only needed for BCF↔XCF conversion —
  not required for our plan, but don't plan around it.
- `tools/HapCUT2` repo clone is intact incl. `utilities/`
  (`calculate_haplotype_statistics.py`, `LinkFragments.py`, `prune_haplotype.py`,
  `hg38.chrom.sizes`) — these ship only in the repo, not in the conda package.

### NEW: `phasing-resources/` on the mount
- `genetic-maps/` — 25 × `chr*.b38.gmap.gz` (GRCh38 SHAPEIT maps), 22 MB, **all 25
  md5-verified**.
- **CORRECTION (same session, after user pointed it out): SHAPEIT5 has a NEW OFFICIAL
  REPO — `https://github.com/odelaneau/shapeit` — and it is LIVE.** Its README states it
  "replaces the original odelaneau/shapeit5 repo, which has been repeatedly disabled by
  GitHub." My earlier note that bioconda was "the only supply line" is superseded.
  The maps stored here were fetched from the SHAPEIT4 repo, but the new repo's
  `resources/maps/b38/genetic_maps.b38.tar.gz` is md5 `5aef9b828bb15ce15920c73988bbad07`
  — **byte-identical** to ours, so they ARE the official SHAPEIT5 maps. No re-fetch needed.
- **New repo also provides:** `resources/chunks/b38/{4cM,20cM}/` chunk coordinates (useful
  for memory-bounded chunked phasing), `test/` simulated data with a runnable
  `test/scripts/phase.array.scaffold.sh` scaffold example, and full source incl.
  `xcftools`. But `static_bins/` holds no committed blobs and there is no tagged release,
  so **bioconda is still the practical binary source**.
- `README.md` — provenance, format, caveats, and the verified **remote panel-slicing**
  recipe. `htslib-tools` htslib 1.14 is built `libcurl=yes` and a live tabix slice
  of the 1000G chr16 panel succeeded, so **the panel is streamed, never downloaded**.
- **FUSE gotcha (new variant):** the first copy pass wrote `chr7.b38.gmap.gz` as a
  **0-byte file** — a genuine silent write failure, *not* read-after-write lag
  (md5 was the empty-file `d41d8cd9…`). Re-copy fixed it. Also, tools that write
  via `tmpfile + rename()` fail on this mount with `EREMOTEIO`; write to `/tmp`
  and `cp` instead. **Lesson: always size-check AND md5-check every file after copy.**

---

## NEXT TASK — sequential phasing: unphased SNPs → HapCUT2 → SHAPEIT5

Plan as stated by the user: take the unphased het SNPs from **both** callers
(`from-bsgenova/`, `from-naive/`), run **HapCUT2** (read-backed) first, then feed
the HapCUT2 output into **SHAPEIT5** (statistical) for a sequential/hybrid phasing.

### Assessment — will it work? Short answer: HapCUT2 yes; SHAPEIT5 stage, very likely not on this data.

**0. Re-assessed against the new official repo (2026-07-28) — verdict UNCHANGED.**
The repo's own glossary (`docs/_glossary/scaffold.md`) defines a scaffold as *"a set of
**highly confident haplotypes**, typically on a subset of the data. SHAPEIT5 uses the
haplotypes derived at **common variants** as haplotype scaffolds onto which heterozygous
genotypes are phased one **rare variant** at a time."* That is precisely the
common-variant→rare-variant use case — **not** read-based blocks of arbitrary relative
orientation. The only shipped scaffold example is `phase_common --input target.unrelated.bcf
--scaffold target.scaffold.bcf --region 1 --map chr1.gmap.gz`. Nothing in the repo
documents or supports a read-backed scaffold, and there is still no `--use-PS`. So the
concerns below stand.

**1. The interface exists, and the order is right.** `SHAPEIT5_phase_common` has
`-S/--scaffold` = *"Scaffold of haplotypes in VCF/BCF/XCF format"*, and it *"can
contain data for a subset of samples at a subset of variant sites"*. HapCUT2 with
`--outvcf 1` emits `<out>.phased.vcf`. So the hand-off is mechanically real, and
read-backed→statistical is the only sensible direction. Note SHAPEIT5 has **no
`--use-PS`** (SHAPEIT4 had it); `--scaffold` is the only entry point.

**2. But HapCUT2's output violates the scaffold contract.** The docs require
scaffold data be *"non-missing and phased"*, and SHAPEIT5 treats a scaffold as
**globally consistent across the chromosome**. HapCUT2 emits **independent blocks**
whose orientations are mutually arbitrary — nothing links block *i* to block *j*.
Dumping a multi-block phased VCF in as a scaffold therefore feeds SHAPEIT5
**hard constraints that are wrong ~50% of the time** for every block after the
first. That degrades results rather than improving them. Mitigation: use **at most
one block per chromosome** as scaffold, and drop HapCUT2's pruned (`./.`) and
unphased sites first.

**3. The deeper problem — SHAPEIT5 has almost nothing to work on here.**
- **N = 3 unrelated donors**, each phased independently (HapCUT2 is one individual
  at a time). Statistical phasing is a population method; with N=3 and no panel
  there is no linkage information at all. A panel is mandatory, not optional.
- **Panel/map coverage is absent exactly where our SNPs are.** Chromosome spread of
  the het calls (checked this session):
  - H1930001 (863 het): chr16 232, **chrM 92**, **chr22_KI270733v1_random 89**, chr17 71, chr1 55, **chrUn_KI270438v1 54**, chr21 54, **chrUn_KI270744v1 49** …
  - H1930002 (289 het): chr16 146, chr17 56, **chrUn_KI270438v1 36**, **chrM 18**, chr1 18 …
  - H1930004 (1366 het): chr16 213, **chrM 158**, chr21 101, chr17 90, **chr22_KI270733v1_random 88**, chr1 71, **chrUn_KI270438v1 68**, **chr17_GL000205v2_random 63** …
  There are **no genetic maps and no panel haplotypes for `chrM`, `chrUn_*`,
  `*_random`, or alt contigs** — SHAPEIT5 cannot even be run on them. And chr16 is
  dominated by the same pericentromeric/segdup locus that S4 already showed has
  **zero** 1000G panel variants.
- **`chrM` het calls are not diploid at all** (92/18/158 per donor). They are
  heteroplasmy or NUMT mismapping. They must be dropped before *both* tools or
  they will pollute the results. Same reasoning as the S4 `chrY` exclusion.
- **This is a ~1× pseudo-bulk callset.** A real genome has ~2 M germline hets; we
  have 289–1366. So these are a tiny, coverage-biased, artifact-enriched subset,
  concentrated precisely where mapping is worst — the opposite of the well-behaved
  common-variant set SHAPEIT5 is designed for. S4 already got **0 sites phased**
  for exactly this reason, and the contig spread above says that will repeat.

**4. What to expect from HapCUT2 alone: modest but real.** At ~1× the odds of two
hets sharing a fragment are low, *but* the calls are clustered (232 on chr16 for
H1930001) because those are high-coverage pileups — so blocks will form. Expect
considerably more than S4's 2 blocks / 5 SNPs, since there are ~100× more sites.

**5. Recommendation.** Run **HapCUT2 as the primary, reportable result**. Keep
SHAPEIT5 as a clearly-labelled *exploratory* arm, restricted to (a) real autosomes
only, (b) sites actually present in the 1000G panel, (c) one HapCUT2 block per
chromosome as the scaffold. If it yields 0 again, that is a legitimate negative
finding about statistical phasing of low-coverage bisulfite single-cell callsets —
not a pipeline failure.

### Concrete prep steps required before either tool runs
1. **Drop non-diploid / unphaseable contigs**: exclude `chrM`, `chrY`, `chrL`; decide
   whether to keep `chrUn_*`/`*_random` (HapCUT2 can use them, SHAPEIT5 cannot).
2. **Keep het only** — `bcftools view -g het` (extractHAIRS defaults `--hom 0`, but
   filter explicitly so both tools see the same site list).
3. **Split multiallelics** — `bcftools norm -m -any`. Our VCFs contain sites like
   `chr1:16644262 C→G,A`; extractHAIRS defaults `--triallelic 0` and would skip them.
4. **extractHAIRS needs an UNCOMPRESSED VCF** (`--VCF <FILENAME> ... (unzipped)`).
   Our callsets are `.vcf.gz` → `bgzip -d` a working copy.
5. **Same VCF for both stages** — HAPCUT2 warns: *"use EXACT SAME file that was used
   for the extracthairs program"*.
6. **Inject `##contig` headers + AC/AN** (S4 had to do this for the bsgenova VCF;
   `--filter-maf` in SHAPEIT5 requires AC/AN).
7. Note **extractHAIRS has NO bisulfite mode** (established in S4's appendix, verified
   three ways). C/T and G/A sites stay bias-prone; A/T-backbone blocks are solid.

### Memory-constrained execution notes (user requirement)
- `extractHAIRS` streams the BAM → low, flat memory. Bound it further with
  `--region chr<N>` and loop over chromosomes rather than one whole-genome pass.
- `HAPCUT2` memory scales with fragments × SNPs-per-fragment — trivial at these
  site counts. `--long_reads` is not needed (short reads).
- `SHAPEIT5` is the only real memory risk, and it scales with **panel size**.
  → **Never download or load a whole-chromosome 1000G panel.** Slice ±1 Mb windows
  remotely via `bcftools view -r` over https (verified working this session).
- **Stage one BAM at a time** in `/tmp`, delete before the next. `/tmp` is currently
  clean (the naive run's ~14 GB of staging was removed in Session 9), and the
  overlay has ~6.2 TB free — but extractHAIRS still only needs one BAM resident.

---

## TOMORROW'S TASK (queued 2026-07-28) — 4 preprocessing steps before HapCUT2

Readiness audit done 2026-07-28. **Everything heavy is already in place:** `hapcut2` env
(extractHAIRS + HAPCUT2 1.3.4), `htslib-tools` env, all 3 merged BAMs **with `.bai`**
(3.9/2.6/4.6 GB), `hg38_chrL.fa`+`.fai`, both callsets valid and unphased, ~6.2 TB free
on `/tmp`. Verified also: **VCF contig ordering already matches the BAM/`.fai` exactly and
every VCF contig exists in the reference → NO `bcftools sort` needed.**

Only the VCF side needs work. Apply to BOTH `from-bsgenova/` and `from-naive/`, all 3 donors:

1. **Inject `##contig` headers** — both callsets have **zero** `##contig` lines. Build them
   from `hg38-reference/hg38_chrL.fa.fai` and reheader. Required for tabix indexing and for
   HAPCUT2 `--outvcf` to emit a valid phased VCF. (Same fix S4 needed.)
2. **Split multiallelics** — `bcftools norm -m -any`. extractHAIRS defaults `--triallelic 0`
   and would **silently skip** these. Matters mostly for bsgenova — het multiallelic sites:
   **57 / 24 / 77** (H1930001/2/4). Naive is barely affected: 2 / 2 / 5.
3. **Decompress for extractHAIRS** — `--VCF` requires an **unzipped** file. `bgzip -d` a
   working copy, and pass the *same* file to HAPCUT2 (it warns: "use EXACT SAME file that
   was used for the extracthairs program").
4. **Drop `chrM` and `chrY`** — neither is diploid; chrM hets are heteroplasmy/NUMT
   artifacts. Counts (het): bsgenova chrM **69/14/146**, chrY 12/9/30; naive chrM 15/3/14,
   chrY 7/5/12. (`chrL` lambda has no calls but exclude defensively.)

### Phaseable input after filtering

| callset | donor | het on primary chr (autosome+X) | additional unplaced/`_random` |
|---|---|---|---|
| bsgenova | H1930001 / H1930002 / H1930004 | 483 / 177 / 812 | +299 / +89 / +378 |
| naive    | H1930001 / H1930002 / H1930004 | 302 / 148 / 331 | +115 / +30 / +124 |

**Open decision:** unplaced/`_random` contigs are 25–35% of the het calls. **HapCUT2 CAN
phase them; SHAPEIT5 CANNOT** (no genetic map, no panel haplotypes). Proposed: keep them
for HapCUT2, drop them only at the SHAPEIT5 stage — but flag their blocks as
low-confidence, since those contigs are where mapping artifacts concentrate.

### Reference panel status (asked 2026-07-28)
**Not held locally — by design, and it would not rescue the SHAPEIT5 stage anyway.**
- Access is **verified working**: `htslib-tools`' htslib 1.14 is built `libcurl=yes`, and a
  live remote tabix slice of the 1000G 30x NYGC chr16 panel returned records. So the panel
  is streamed in ±1 Mb windows, never downloaded (multi-GB per chromosome = the memory
  blow-up we are avoiding). Recipe in `phasing-resources/README.md`.
- **But S4 already established our SNPs are essentially absent from that panel**, and the
  contig spread confirms it: the calls concentrate on chr16 pericentromeric/segdup (zero
  panel variants there), `chrM`, and unplaced/`_random` contigs that panels do not cover at
  all. Obtaining a bigger panel does not fix this — these are low-coverage de-novo/rare
  single-cell calls, not catalogued common variants.

### Run-command templates (verified flags, NOT yet executed)
```bash
source /opt/conda/etc/profile.d/conda.sh && conda activate hapcut2
# stage ONE BAM at a time to /tmp; loop chromosomes to bound memory
extractHAIRS --bam <donor>.bam --VCF <prepped.vcf> --region chr<N> --out frag.chr<N>
cat frag.chr*  > frags.all                       # fragments are concatenable
HAPCUT2 --fragments frags.all --VCF <prepped.vcf> --output <donor>.blocks --outvcf 1
```
- `--outvcf 1` is **required** — it produces `<output>.phased.vcf`, the only thing that
  could later be handed to `SHAPEIT5_phase_common --scaffold`.
- extractHAIRS defaults already suit us: `--hom 0` (het only), `--mbq 13`, `--mmq 20`,
  `--maxIS 1000`. Do NOT set `--hic/--pacbio/--ont` (short Illumina reads).
- Block stats afterwards: `tools/HapCUT2/utilities/calculate_haplotype_statistics.py`
  (repo-only, not in the conda package; needs the `hapcut2` env's python+pysam).

### Genetic-map caveat to remember
The maps' 2nd column is `16`, but our data uses `chr16`. Pass `--region chr16` explicitly;
if SHAPEIT5 rejects the map, rewrite col 2 with
`awk 'NR==1{print;next}{$2="chr"$2;print}'`. Full note in `phasing-resources/README.md`.

---

## END-OF-SESSION STATE (2026-07-28) — resume checklist for the next agent

**Nothing was executed this session beyond env builds and read-only inspection.**
No SNP caller, no HapCUT2, no SHAPEIT5 was run on the data.

READY (verified this session):
- 4 conda envs: `bsgenova`, `htslib-tools`, `hapcut2`, `shapeit5` — **but the node is
  EPHEMERAL, so expect to rebuild with `bash environments/00_build_all.sh` first.**
- 3 merged BAMs + `.bai`; `hg38-reference/hg38_chrL.fa` + `.fai`.
- Both callsets complete and unphased: `from-bsgenova/` (863/289/1366 het),
  `from-naive/` (439/186/481 het).
- `phasing-resources/genetic-maps/` — 25 official b38 maps, md5-verified.
- Remote 1000G panel slicing confirmed working (stream, never download).

TODO next session, in order:
1. The 4 preprocessing steps above (both callsets, 3 donors each).
2. Run HapCUT2 per donor per callset → blocks + `.phased.vcf`.
3. Only then consider the SHAPEIT5 exploratory arm (autosomes only, panel-present sites,
   ≤1 block per chromosome as scaffold). A 0-site result is a legitimate negative finding.

FUSE rules that bit us this session (all three are real, all cost time):
- `Write`-style `tmpfile + rename()` fails with **`EREMOTEIO`** → write to `/tmp`, then `cp`.
- A copy silently produced a **0-byte file** (`chr7.b38.gmap.gz`) — not lag, a real failed
  write. **Always size-check AND md5-check after every copy.**
- Directory listings lag: a folder that looks empty may have been populated minutes ago.
  Re-check with `find`/`stat` before concluding anything about the mount's state.
- Two agents were appending to this log concurrently — hence two "Session 9" headings
  (mine renumbered to 10). Check `stat -c%s` before and after any whole-file rewrite.
## Session 11 (2026-07-29) — TASK ASSIGNMENT (recorded so it survives disruption)

User's task:
1. Read this log to catch up.
2. Create `workspace/shell-scripts/HapCUT2_preprocessing.sh` with the 4-step HapCUT2
   preprocessing workflow queued on 2026-07-28.
3. Create `bisulfite-aware-SNPs/{from-bsgenova,from-naive}-HapCUT2_preprocessed/`.
   (User wrote `workspace/data/bisulfite-aware-SNPs/...`; real tree is under
   `workspace/data/snMC-seq_SNP_phasing_experiment/` — new folders placed there as
   siblings of the existing `from-bsgenova/`+`from-naive/`.)
4. Run the preprocessing script on the unphased het SNPs in `from-bsgenova/` +
   `from-naive/` CONCURRENTLY, resource-capped (don't disrupt the VSCode-server);
   outputs to the matching `*_preprocessed/` folder.
5. Run HapCUT2 on the preprocessed SNPs CONCURRENTLY (resource-capped); phased outputs
   to `phased-from-HapCUT2/from-bsgenova/` or `.../from-naive/` by source.

The 4 preprocessing steps: (1) inject `##contig` from `hg38_chrL.fa.fai`;
(2) split multiallelics `bcftools norm -m -any`; (3) keep het only + drop chrM/chrY/chrL;
(4) decompress to plain `.vcf` for extractHAIRS (same file to HAPCUT2).

Status: DONE — see "Session 11 RESULTS" below.
[Note: this block was rewritten after a FUSE write corrupted the original into whitespace.]
### Session 11 RESULTS (2026-07-29) — preprocessing + HapCUT2 done, both callsets ✅

**Deliverables created:**
- `workspace/shell-scripts/HapCUT2_preprocessing.sh` — the 4-step preprocessor for ONE VCF
  (STEP 1 inject `##contig` from `.fai` via `bcftools annotate --header-lines`; STEP 2
  `bcftools norm -m -any`; STEP 3 `bcftools view -g het -t ^chrM,chrY,chrL`; STEP 4 emit
  plain `.vcf` for extractHAIRS + a bgzip/tabix copy). FUSE-safe (works in `/tmp`, copies
  + verifies). Writes `<name>.preprocessed.vcf(.gz/.tbi)` + a per-step `.preprocessing.log`.
- `workspace/shell-scripts/HapCUT2_run.sh` — companion runner (extractHAIRS → HAPCUT2
  `--outvcf 1`) for ONE (callset,donor); also FUSE-safe.
- Output folders created: `bisulfite-aware-SNPs/{from-bsgenova,from-naive}-HapCUT2_preprocessed/`
  and `phased-from-HapCUT2/{from-bsgenova,from-naive}/`.
  (NOTE the real tree is under `data/snMC-seq_SNP_phasing_experiment/`; the user's
  `workspace/data/bisulfite-aware-SNPs/...` shorthand maps here.)

**Env (node ephemeral — only `base` survived):** rebuilt `htslib-tools` (bcftools 1.14) and
`hapcut2` (extractHAIRS + HAPCUT2 1.3.4) via `environments/0{2,3}_*.sh`. Node: 256 cores,
~480 GiB free RAM, 6 TB `/tmp`.

**Preprocessing (all 6 VCFs, run concurrently — trivial load):** het records kept after
split + chrM/chrY/chrL drop:

| callset  | input→split | H1930001 het | H1930002 het | H1930004 het |
|----------|-------------|-------------:|-------------:|-------------:|
| bsgenova | 1166→1256 / 421→465 / 1712→1827 | 792 | 271 | 1211 |
| naive    | 660→663 / 306→312 / 688→697      | 419 | 180 | 457  |

(Split adds records because bsgenova multiallelics — 44 split for H1930002 etc. — become
biallelic; het counts are a touch above the 2026-07-28 estimate because 1/2 sites split
into two hets. All outputs verified: no chrM/chrY/chrL, no multiallelic, all het.)

**HapCUT2 (all 6 run CONCURRENTLY, resource-capped — service undisturbed):** staged the 3
merged BAMs+`.bai` (~11 GB) to `/tmp` (background cp; foreground timed out at 2 min = known
iRODS slowness), then launched all 6 with `setsid nohup` + `nice -n 10` (survives
disconnects, the S8 death cause). extractHAIRS is single-threaded + streams the BAM (flat
low RAM); load average peaked ~2.5, ~154 GiB RAM still free — VSCode-server unaffected.
All 6 `EXIT_0` (`ALL_HAPCUT2_FINISHED overall_rc=0`) in a few minutes.

| callset  | donor    | fragments | blocks | phased SNPs (`\|`) |
|----------|----------|----------:|-------:|-------------------:|
| bsgenova | H1930001 | 6276 | 109 | 359 |
| bsgenova | H1930002 | 2234 |  44 | 129 |
| bsgenova | H1930004 | 6005 | 132 | 382 |
| naive    | H1930001 | 8207 |  61 | 215 |
| naive    | H1930002 | 2456 |  31 | 106 |
| naive    | H1930004 | 7244 |  72 | 250 |

- Block-file phased-SNP counts == `|`-genotype counts in each `.blocks.phased.VCF` (internal
  consistency ✓). Each phased VCF re-parses with bcftools and carries the FULL het site list
  (792/271/1211; 419/180/457) with phaseable sites marked `|`, the rest left `/`.
- As predicted vs S4 (2 blocks/5 SNPs on 1 cell): ~100× more sites → dozens–hundreds of
  blocks. bsgenova phases more SNPs than naive (bigger het input); naive yields MORE
  fragments on H1930001/H1930004 but fewer phased SNPs (fewer, sparser het sites to link).
- Caveat unchanged: extractHAIRS has NO bisulfite mode → C/T & G/A sites bias-prone;
  A/T-backbone blocks solid.

**Outputs** → `phased-from-HapCUT2/from-{bsgenova,naive}/`, per donor: `.fragments`,
`.blocks`, `.blocks.phased.VCF`, `.extractHAIRS.log`, `.HAPCUT2.log` (5×3=15 files each dir,
all size-verified on mount). `/tmp/hapcut2_run` removed.

**Next steps:** optional SHAPEIT5 exploratory arm (autosomes only, panel-present sites, ≤1
HapCUT2 block/chr as scaffold — S4/S10 expect ~0); optional bsgenova-vs-naive phased-block
concordance comparison.

---

## Session 12 (2026-07-29) — Can HapCUT2 output feed SHAPEIT5 for refinement? Assessed + confirmed.

**Question (from PI):** take the HapCUT2 `.blocks.phased.VCF` and run it through SHAPEIT5
to get more refined phasing / phase more of the unphased SNPs.

**Confirmed against the actual binary** (`SHAPEIT5_phase_common --help`, env rebuilt this session):
- Options: `-I/--input`, `-H/--reference` (panel), `-S/--scaffold` (VCF/BCF), `-M/--map`,
  `-R/--region`, `--pedigree`. **No `--use-PS`** (grep count 0). So the ONLY way to inject
  HapCUT2's read-backed phasing is via `--scaffold`. Our phased VCFs DO carry a `PS`
  (phase-set/block-id) tag, but SHAPEIT5 cannot consume `PS` directly.

**Verdict: mechanically possible, but it will NOT refine THIS data — expect ~0 extra sites.**
The read-backed→statistical direction is sound in general (it's the standard hybrid), and the
`--scaffold` interface is real. But three data facts block any gain here:

1. **Scaffold-contract mismatch.** SHAPEIT5 treats a scaffold as haplotypes that are
   *globally consistent across the whole chromosome*. HapCUT2 emits MANY independent blocks
   per chromosome whose relative orientation is arbitrary (e.g. bsgenova/H1930001 has 113
   phased SNPs on chr16 spread over many blocks). Feeding a multi-block VCF as scaffold
   injects hard constraints that are wrong ~50% of the time for every block after the first
   → it can *degrade*, not refine. Only safe if reduced to <=1 block per chromosome.

2. **SHAPEIT5 is a population/panel method; we have neither.** It refines by borrowing
   linkage from a reference panel (`--reference`) or many samples. We have 3 unrelated donors
   phased independently, and our SNPs are rare/de-novo low-coverage single-cell calls that are
   essentially absent from the 1000G panel. No panel overlap => no statistical info to add.
   (S4 already got 0 sites phased for exactly this reason.)

3. **No map / no panel where our phased SNPs actually are.** The `|` genotypes concentrate on
   chr16 pericentromeric/segdup, `chr22_*_random`, `chrUn_*`, `chrM`, alt/`_random` contigs
   (this session: phased SNPs span 31 contigs, top = chr16 113, chr22_KI270733v1_random 62,
   chrUn_KI270744v1 28 ...). There are NO genetic maps and NO panel haplotypes for those
   contigs => SHAPEIT5 cannot even run on them.

**What SHAPEIT5 would need to help:** common variants that exist in a population panel, on
mappable autosomes, with a panel supplied via `-H`. Our callset is the opposite of that.

**Recommendation (unchanged from S10):** keep **HapCUT2 as the reportable phasing result**.
A SHAPEIT5 arm is only worth running as an explicitly-labelled *exploratory / negative-control*
step: autosomes only, sites actually present in the 1000G panel, `-H` = a remote-sliced panel
window, and <=1 HapCUT2 block/chromosome as `--scaffold`. A 0-site result there is a legitimate
scientific finding about statistical phasing of low-coverage bisulfite single-cell callsets —
not a tooling failure. Bottom line for the PI: the idea is right in principle, but this
particular data lacks the panel-catalogued common variants that statistical refinement relies on.

### Session 12 (cont., 2026-07-29) — SNP-count + phasing summary (bsgenova vs naive)

Pseudo-bulk CX45 (10 cells/donor, 3 donors), bisulfite-aware SNP calling then HapCUT2 phasing.

**(A) Raw callset — homozygous vs heterozygous SNPs (per donor):**

| approach | donor    | total | het  | hom |
|----------|----------|------:|-----:|----:|
| bsgenova | H1930001 | 1166  | 869  | 297 |
| bsgenova | H1930002 |  421  | 292  | 129 |
| bsgenova | H1930004 | 1712  | 1378 | 334 |
| bsgenova | **sum**  | **3299** | **2539** | **760** |
| naive    | H1930001 |  660  | 440  | 220 |
| naive    | H1930002 |  306  | 190  | 116 |
| naive    | H1930004 |  688  | 485  | 203 |
| naive    | **sum**  | **1654** | **1115** | **539** |

(naive ~half the calls of bsgenova — it hard-masks conversion-ambiguous bases, so it is the
more conservative caller.)

**(B) HapCUT2 phasing of the heterozygous SNPs.** Input to HapCUT2 = the *preprocessed* het
set (after splitting multiallelics + dropping chrM/chrY/chrL), so it is a little below the raw
het count above. Of that input: phased = genotype written `A|B`; unphased = left `A/B`.

| approach | donor    | het into HapCUT2 | phased | unphased |
|----------|----------|-----------------:|-------:|---------:|
| bsgenova | H1930001 | 792  | 359 | 433 |
| bsgenova | H1930002 | 271  | 129 | 142 |
| bsgenova | H1930004 | 1211 | 382 | 829 |
| bsgenova | **sum**  | **2274** | **870** | **1404** |
| naive    | H1930001 | 419  | 215 | 204 |
| naive    | H1930002 | 180  | 106 |  74 |
| naive    | H1930004 | 457  | 250 | 207 |
| naive    | **sum**  | **1056** | **571** | **485** |

**Why some het SNPs stayed unphased (one sentence):** at ~1x pseudo-bulk depth most het sites
are isolated, so no single read or read-pair covers two het sites at once — there is no
read-backed link to place them in a haplotype block, so they remain unphased.

**(C) Chromosomes of phased vs unphased SNPs (union across the 3 donors):**
- **Phased** — concentrate on a few high-coverage clustered loci: dominated by **chr16**
  (bsgenova 278 / naive 321), then chr17, chr21, chr3, chr10, and a handful of unplaced/random
  contigs that carry local pileups (chr22_KI270733v1_random, chrUn_KI270438v1, chrUn_GL000220v1,
  chrUn_KI270744v1). Phasing happens where sites cluster densely enough to share reads.
- **Unphased** — spread thinly across **nearly all autosomes** (chr1, chr2, chr4, chr5, chr7,
  chr9, chr10, chr12, chr15, chr17, chr19, chr20, chr21 ...) plus a **long tail of chrUn_* /
  *_random / alt contigs**, most with only 1-2 sites each — i.e. isolated singletons with no
  read-neighbour to link to. (chr16 also has an unphased remainder: 133 bsgenova / 88 naive.)

**Does more sequencing depth (merging more snMC-seq cells) help use SHAPEIT5? (one sentence):**
No — deeper coverage would let HapCUT2 phase MORE sites (more reads spanning multiple hets),
but it would NOT make SHAPEIT5 effective, because SHAPEIT5's bottleneck is whether the variants
are common and catalogued in its reference panel + genetic maps, which sequencing depth does
not change.

### Next-session TODO (queued 2026-07-29)
- **Make a git repository of the `workspace/` directory** (init at `workspace/`, commit the
  scripts/data-provenance/logs). NOTE a partial/empty `workspace/.git` skeleton already exists
  (created this session, no HEAD/objects) — inspect and either finish or remove+re-init.
- `workspace/conversation-logs/.git` was **removed this session** (per user) so the logs are NOT
  tracked as their own nested repo; when the `workspace/` repo is created, keep conversation-logs
  out of tracking (e.g. add it to `.gitignore`).
