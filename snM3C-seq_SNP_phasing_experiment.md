# snM3C-seq SNP Phasing Experiment — Working Log

**This is the fresh log for the *snM3C-seq* SNP-phasing experiment. New sessions append here.**

snM3C-seq = single-nucleus **methyl-3C** — the combined assay that captures **both** DNA
methylation (snmC-seq3 chemistry) **and** 3D-genome chromatin contacts (Hi-C / "3C") in the
same nucleus. It is the sibling assay to the plain snMC-seq we worked on previously, from the
**same** Tian W, Zhou J, … Ecker JR, *Science* 2023 human-brain study.

The section below summarizes the completed **snMC-seq** experiment so a future agent has the
full foundation. The full blow-by-blow of that prior work lives in
`snMC-seq_SNP_phasing_experiment.md` (12 sessions, 2026-07-12 → 2026-07-29).

---

## Why move to snM3C-seq — the scientific point (read this first)

The plain snMC-seq experiment reached a clear conclusion: **read-backed phasing (HapCUT2) is
the only usable method** for this data, and **statistical phasing (SHAPEIT5) does not work** —
the SNPs are rare/de-novo low-coverage single-cell calls essentially absent from population
reference panels, so there is no linkage information to borrow.

The bottleneck for HapCUT2 was **read span**: at ~1× pseudo-bulk depth most het sites are
isolated, so no single short read or read-pair covers two het sites at once → no read-backed
link → the site stays unphased (only ~35–45% of het sites phased in snMC-seq).

**This is exactly what the "3C" in snM3C-seq should fix.** Hi-C/3C reads carry **long-range
chromatin contacts** — read pairs whose two ends can be tens of kb to Mb apart. Those
long-range links can bridge het sites that short reads never co-cover, so read-backed phasing
should reach dramatically more sites and build much larger blocks. Concretely:

- `extractHAIRS` has a dedicated **`--hic 1`** mode built for exactly this (Hi-C linked reads),
  which we deliberately did NOT use on snMC-seq (short Illumina reads → standard mode).
- HapCUT2 was originally designed with Hi-C phasing as a first-class use case.

So the hypothesis for this experiment: **the same pipeline, but on snM3C-seq data with
extractHAIRS in Hi-C mode, will phase substantially more SNPs into far larger haplotype
blocks** than the plain snMC-seq result — potentially chromosome-scale.

---

## What the prior snMC-seq experiment established (foundation)

### Dataset & source
- Paper: **Tian W, Zhou J, … Ecker JR, "Single-cell DNA methylation and 3D genome architecture
  in the human brain," *Science* 2023** (DOI 10.1126/science.adf5357, PMID 37824674).
- GEO superSeries **GSE215353**. 517k cells, 3 adult male donors: **H1930001, H1930002,
  H1930004** (note: NOT H1930003). 46 brain dissection regions, 188 (methylation-cluster) cell
  types.
- **GEO/SRA host NO BAMs** — only raw FASTQ (SRA) and processed ALLC methylation tables. The
  aligned **BAMs live only at the NeMO Archive** (`nemo:dat-jx4eu3g` → alignments
  `dat-budp3jg`), one per-cell tarball each:
  - snMC-seq: `<cell>.final.bam.tar`
  - **snM3C-seq: `<cell>.3C.sorted.bam.tar`** — cell names contain **`_3C_`** (this is how you
    filter the manifest for the 3C assay; the plain snMC-seq work filtered `_3C_` OUT).
- Download URL pattern (worked for snMC-seq; expect the analogous path for 3C):
  `https://data.nemoarchive.org/biccn/grant/u01_ecker/ecker/epigenome/sncell/mCseq/human/processed/align/<cell>.<...>.bam.tar`
- Manifest kept on the mount: `nemo_align_manifest.tsv` (524,004 rows, md5
  `954359ea98c3006364cdc7f1233ad012`) + `nemo_align_manifest.README`. The manifest mixes both
  assays plus stray mouse-3C rows and a `Methylome` comment row — filter carefully.

### Alignment & reference (applies to 3C too)
- Aligner = **Bismark `--pbat`** (bisulfite). All BAMs share an identical **456-contig `@SQ`** =
  UCSC hg38 (chr1…chrY + randoms/alts) **+ `chrL` = lambda spike-in** (NC_001416.1, 48,502 bp,
  bisulfite-conversion control).
- Reference built and kept on the mount: **`hg38_chrL.fa`** (3,273,530,352 B, md5
  `612ab3bbf05a3f916093e5cd11b31d5a`) + `.fai` = UCSC hg38 + lambda→`chrL`. Contig
  names+lengths diff-match the BAM `@SQ` exactly.
- **Reference-genome decision (important):** the authors aligned each donor to its OWN per-donor
  **SNP-substituted (personalised) hg38** (`@PG` shows `.../hba-donor/hXXXXXXX`). We do **NOT**
  have those and do **NOT** want them for SNP calling. **Call germline variants against standard
  hg38+chrL**, because a donor-substituted reference has the donor's alleles baked in → germline
  SNPs would be suppressed. Standard hg38 surfaces the germline variants (the goal). Bonus: since
  reads were mapped to a personalised ref, ref-mapping bias at het sites is already reduced.

### The SNP-calling + phasing pipeline (proven on snMC-seq)
1. **Bisulfite-aware SNP calling**, two callers, kept for comparison:
   - **bsgenova** (`workspace/tools/bsgenova`) — probabilistic bisulfite-aware caller. Keeps all
     read evidence via a conversion model → more calls.
   - **naive caller** (`workspace/Python-scripts/bisulfite_aware_naive_SNP_caller.py`, v1.0.0) —
     written by us; emits bsgenova-compatible output. Uses the **hickit.js masking rule**
     (Watson `T` and Crick `A` are conversion-ambiguous → masked); flat error model, no
     methylation transition matrix. Hard-masks ambiguous bases → ~half as many, more
     conservative calls.
   - Both emit `<prefix>.vcf.gz` (`GT:GQ:GQH:DP:DPW:DPC`) + `<prefix>.snv.gz` (fuller per-site
     table). All calls come out **unphased** (`0/1`, no `|`).
2. **4-step preprocessing before phasing** (script:
   `workspace/shell-scripts/HapCUT2_preprocessing.sh`):
   (1) inject `##contig` headers from `hg38_chrL.fa.fai` (callsets have none, breaks tabix &
   HAPCUT2 `--outvcf`); (2) split multiallelics `bcftools norm -m -any` (extractHAIRS defaults
   `--triallelic 0` → silently skips them); (3) keep het only + drop `chrM,chrY,chrL`
   (`bcftools view -g het -t ^chrM,chrY,chrL` — none are diploid; chrM hets = heteroplasmy/NUMT
   artifacts); (4) decompress to a plain `.vcf` for extractHAIRS (must pass the SAME file to
   HAPCUT2).
3. **Read-backed phasing — HapCUT2** (script: `workspace/shell-scripts/HapCUT2_run.sh`):
   `extractHAIRS --bam … --VCF … --region chr<N>` (loop chromosomes to bound memory; fragments
   concatenable) → `HAPCUT2 --fragments … --VCF … --outvcf 1` → `.blocks` + `.blocks.phased.VCF`.
   - **For snM3C-seq, this is where you add `extractHAIRS --hic 1`** (the whole point — see top).
   - Caveat carried from snMC-seq: standard extractHAIRS has **NO bisulfite mode** (verified 3
     ways) → C/T & G/A sites are bias-prone; A/T-backbone blocks are solid.
4. **Statistical phasing — SHAPEIT5: assessed, NOT usable on this data** (0 sites expected). No
   `--use-PS`; only entry is `--scaffold`, and HapCUT2's many independent per-chromosome blocks
   violate the scaffold's global-consistency contract. Our SNPs are absent from the 1000G panel
   and sit on unmappable contigs (chr16 pericentromeric/segdup, chrM, chrUn_*, *_random). Keep
   SHAPEIT5 only as an explicitly-labelled exploratory/negative-control arm. **Deeper coverage
   would NOT rescue SHAPEIT5** (its bottleneck is panel-catalogued common variants, not depth).

### snMC-seq final numbers (pseudo-bulk CX45, 10 cells/donor merged with `samtools merge`)
Raw SNP calls (total / het):
| caller   | H1930001 | H1930002 | H1930004 |
|----------|---------:|---------:|---------:|
| bsgenova | 1166 / 869 | 421 / 292 | 1712 / 1378 |
| naive    | 660 / 440  | 306 / 190 | 688 / 485   |

HapCUT2 phasing (het into HapCUT2 → phased / unphased):
| caller   | H1930001 | H1930002 | H1930004 |
|----------|---------:|---------:|---------:|
| bsgenova | 792 → 359 / 433 | 271 → 129 / 142 | 1211 → 382 / 829 |
| naive    | 419 → 215 / 204 | 180 → 106 / 74  | 457 → 250 / 207  |

Phased SNPs concentrate on a few high-coverage clustered loci (chr16 dominant, then chr17,
chr21, chr3, chr10, some unplaced/random contigs). Unphased SNPs spread thinly as isolated
singletons across nearly all autosomes + a long tail of chrUn_*/*_random/alt contigs. **This
sparsity is the ~1× short-read limitation that snM3C-seq's long-range contacts should relieve.**

---

## Environment & infrastructure — CRITICAL for every session

### The compute node is EPHEMERAL
Conda envs do **NOT** survive between sessions — only the `workspace/` iRODS mount persists.
**Rebuild envs each session.** Recipes live in `workspace/environments/`:
- `00_build_all.sh` (rebuilds all four), `01_bsgenova.sh` (python 3.11 + numpy + pysam +
  samtools — also serves the naive caller), `02_htslib-tools.sh` (bcftools/samtools/htslib 1.14
  + gsl 2.6), `03_hapcut2.sh` (hapcut2 1.3.4 + samtools + pysam), `04_shapeit5.sh` (shapeit5
  5.1.1 + libboost 1.85). `FORCE_RECREATE=1` to wipe+rebuild.
- **KEY channel fix** (in `_common.sh`): the node's base conda lists only `defaults`; appending
  `-c bioconda -c conda-forge` alone mixes ABIs and every env fails to solve. Always use
  `--override-channels --strict-channel-priority -c conda-forge -c bioconda`.
- `conda` is not on PATH in non-login shells → `source /opt/conda/etc/profile.d/conda.sh` first.
- `SHAPEIT5_xcftools` is missing from the bioconda 5.1.1 build (only needed for BCF↔XCF; not in
  our plan). SHAPEIT5 has a live official repo now: `https://github.com/odelaneau/shapeit`.

### The `workspace/` mount is FUSE iRODS (CyVerse) — slow & lossy for large writes
- **Always stage large files in local `/tmp`, run there, then `cp` to the mount and verify.**
- Verify **both size AND md5** after every copy — copies can silently produce a **0-byte file**
  (real failed write, happened to `chr7.b38.gmap.gz`), and immediate post-`cp` checks give false
  mismatches from read-after-write lag (re-read after a short settle).
- Tools that write via `tmpfile + rename()` fail with **`EREMOTEIO`** → write to `/tmp`, then
  `cp`.
- Directory listings lag: a folder that looks empty may have been populated minutes ago. Re-check
  with `find`/`stat`, not the immediate `ls`, before concluding anything about mount state.
- iRODS reads are slow → stage with a **background** `cp` (foreground often times out at ~2 min).
- Run detached with `setsid nohup` + `nice -n 10` so jobs survive disconnects (a disconnect was
  the suspected cause of an early silent job death) and don't disrupt the VSCode-server. Cap
  worker count (e.g. `-P 16` each) to leave cores free.

### Data locations (under `workspace/data/snMC-seq_SNP_phasing_experiment/`)
- `hg38-reference/` — `hg38_chrL.fa` + `.fai` + `README.md`.
- `Science-snMC-seq/` — per-donor BAM folders `H1930001/ H1930002/ H1930004/`, `merged-BAM/`
  (3 per-donor pseudo-bulk merged BAMs + `.bai`), `nemo_align_manifest.tsv` + `.README`.
- `bisulfite-aware-SNPs/from-bsgenova/`, `.../from-naive/` — raw unphased callsets.
- `bisulfite-aware-SNPs/{from-bsgenova,from-naive}-HapCUT2_preprocessed/` — preprocessed VCFs.
- `phased-from-HapCUT2/{from-bsgenova,from-naive}/` — per donor: `.fragments`, `.blocks`,
  `.blocks.phased.VCF`, `.extractHAIRS.log`, `.HAPCUT2.log`.
- Scripts: `workspace/shell-scripts/HapCUT2_preprocessing.sh`, `HapCUT2_run.sh`;
  `workspace/Python-scripts/bisulfite_aware_naive_SNP_caller.py`; envs in
  `workspace/environments/`; phasing aux in `workspace/phasing-resources/` (25 md5-verified b38
  genetic maps + a panel-slicing README; 1000G panels are streamed via remote tabix, never
  downloaded). Tools cloned in `workspace/tools/` (`bsgenova`, `HapCUT2`).

---

## snM3C-seq experiment — status

**NOT STARTED YET.** This log was created 2026-08-12 to open the fresh experiment. Nothing has
been downloaded, called, or phased for snM3C-seq at time of writing.

### Suggested first steps for a future agent (to be confirmed with the user)
1. Filter `nemo_align_manifest.tsv` for **snM3C-seq** rows (cell names containing `_3C_`,
   tarballs `*.3C.sorted.bam.tar`); pick cells for the 3 donors (region choice per user — CX45
   was used for snMC-seq). Download → md5-verify vs manifest → extract → verify.
2. Confirm the 3C BAM `@SQ` still matches `hg38_chrL.fa` (expected: same 456-contig hg38+chrL).
3. Reuse the exact SNP-calling pipeline (bsgenova + naive, standard hg38+chrL) → unphased het
   callsets.
4. Preprocess (same 4 steps), then run **HapCUT2 with `extractHAIRS --hic 1`** — the key
   difference from the snMC-seq run — and compare phased-SNP count / block sizes against the
   snMC-seq baseline in the table above.

---

## Session log (snM3C-seq)

### Session 1 (2026-08-12) — log created
- Read the full snMC-seq log (`snMC-seq_SNP_phasing_experiment.md`, 12 sessions) and the folder
  `README.md`.
- Created this file as the working log for the new snM3C-seq experiment, seeded with a summary
  of the completed snMC-seq foundation (dataset, reference decision, two-caller pipeline,
  HapCUT2 read-backed phasing, SHAPEIT5-not-applicable finding, all env/FUSE gotchas, data
  locations, and final snMC-seq numbers).
- No data work performed yet.

### Session 2 (2026-08-12) — snM3C-seq data availability verified + download plan

**Task:** confirm whether the snM3C-seq (Bismark-aligned BAM) counterpart of our 30 CX45
snMC-seq cells can be downloaded; set up folders; estimate download time. Nothing downloaded.

#### Folders created
- `Science-snM3C-seq/` (already made) + donor subfolders `H1930001/ H1930002/ H1930004/`.

#### Data availability — VERIFIED against the NeMO manifest + live server
- **Same *nucleus*: NOT available.** snMC-seq and snM3C-seq are two **separate assays run on
  different nuclei** — a nucleus is consumed by one assay only. Confirmed directly: **none** of
  our 30 snMC-seq cell barcodes reappear as a `_3C_` cell in the manifest.
- **Same *donor + brain region (CX45)*: available in abundance.** snM3C-seq Bismark BAMs exist
  at NeMO for all 3 donors:

  | donor | snM3C-seq CX45 cells | CX45 sub-regions (3C) | our snMC-seq CX45 sub-regions | overlap |
  |-------|---------------------:|-----------------------|-------------------------------|---------|
  | H1930001 | 11,574 | NAC, CaB, FI, A46 | NAC, CaB | NAC, CaB ✓ |
  | H1930002 | 5,500  | A24, A46          | A25, A38 | none |
  | H1930004 | 2,989  | Pir               | A24, A44, Pir | Pir ✓ |

  Sub-region overlap is **irrelevant for germline SNP phasing** (germline variants are identical
  across all cells of a donor) — just need same-donor CX45 cells, which all 3 have.
- The manifest holds **137,051** `_3C_` rows total (of 524,004). For 3C rows, **column 1 already
  contains the full `.3C.sorted.bam.tar` suffix** (snMC-seq rows store the bare cell name).

#### KEY: the 3C download URL path is DIFFERENT from snMC-seq
- snMC-seq: `.../epigenome/sncell/mCseq/human/processed/align/<cell>.final.bam.tar`
- **snM3C-seq: `.../epigenome/sncell/m3C-seq/human/processed/align/<col1>`** where `<col1>` is
  the manifest column-1 string (already ends in `.3C.sorted.bam.tar`). The naive
  `mCseq/.../align/` path gives **404** for 3C files — must use `m3C-seq/`. Verified 200/OK and
  directory-listable at the `m3C-seq/` path.

#### Chosen download set (10 largest CX45 3C tars per donor, matching the snMC-seq strategy)
Selecting the 10 largest per donor for max pseudo-bulk coverage (as in snMC-seq Session 6):
- H1930001: 8.18 GiB (per-cell ~814–878 MB; sub-regions NAC + A46)
- H1930002: 8.54 GiB (per-cell ~832–928 MB; sub-region A24)
- H1930004: 7.93 GiB (per-cell ~744–909 MB; sub-region Pir)
- **TOTAL 30 tars = 24.65 GiB** (26,471,403,520 bytes). Note 3C tars are ~2× larger than the
  snMC-seq CX45 tars (~800–900 MB vs ~180–570 MB) — the extra Hi-C/3C contact reads.

#### Measured transfer speeds (live probes, this session)
- **NeMO download**: single stream **~25 MB/s** (25.3 MB/s over a 300 MB range fetch); **3
  concurrent streams ~40 MB/s aggregate** (server appears to cap total near ~40 MB/s, so
  parallelism gives ~1.6×, not 3×).
- **iRODS mount write** (the dominant bottleneck): **~8 MB/s** single-stream (800 MB test, 98 s).
  The extracted BAMs (~24 GiB) must be written here for persistence (ephemeral node).

#### Download-time estimate → bottom line ~60–80 min end-to-end
- Network download to `/tmp`: ~11 min (3 parallel streams) to ~18 min (single stream).
- md5-verify + `tar -xf` in `/tmp` (local disk): ~5 min.
- **Copy ~24 GiB of BAMs to the iRODS mount @ ~8 MB/s: ~50 min** ← dominates the wall time.
- Total realistic: **~65–80 min**, mostly the slow mount write (not the download).

#### DOWNLOAD COMPLETE (2026-08-12) — 30/30 BAMs on the mount ✅
- Ran the downloader (setsid nohup): 3 concurrent downloads → size-check →
  md5-verify vs manifest → `tar -xf` → `samtools quickcheck` → **single serial writer** to the
  mount → verify. samtools 1.9 env built for the checks.
- **Recipe persisted:** `workspace/shell-scripts/download_snM3C-seq_CX45.sh` (the exact script,
  parameterised by `DONORS`/`REGION`/`TOPN`, with the mount-verify settle raised per the FUSE
  lesson below). The original ran from `/tmp` which was cleaned, so this saved copy is the
  reusable reference. samtools env recipe: `mamba create -n samtools --override-channels
  --strict-channel-priority -c conda-forge -c bioconda samtools=1.9`.
- Tar internal layout: `<cell>.3C.sorted/<cell>.3C.sorted.bam`.
- **Result: all 30 present and integrity-verified.** Per donor (10 each):
  H1930001 8.18 GiB, H1930002 8.54 GiB, H1930004 7.93 GiB — **TOTAL 24.65 GiB**, exactly the
  predicted size. BAMs live in `Science-snM3C-seq/{H1930001,H1930002,H1930004}/` as
  `<cell>.3C.sorted.bam`. **No `.bai` built yet.**
- **FUSE gotcha (again):** 6 files were flagged `FAIL_MOUNT` by the copier's inline verify —
  the mount md5 hadn't settled within the 3×4 s retry window. All 6 were **re-verified afterward
  and their md5 matches the staged source exactly** — pure read-after-write lag, not real
  corruption. Lesson reaffirmed: give mount md5-verification a longer settle, and note that
  md5summing an ~850 MB file *off* the mount takes ~40–60 s (a 6-file loop can exceed a 2-min
  Bash timeout). `/tmp/m3c_dl` removed after verification.
- Selected cells (largest-10 per donor): H1930001 = NAC + A46 sub-regions; H1930002 = A24;
  H1930004 = Pir. (Sub-region mix is irrelevant for germline SNP phasing.)

#### Next steps
- **[TODO — NEXT SESSION] Create a GitHub repository for the Claude-agent recipe markdown log
  files.** Purpose: version-control / share the recipe + experiment logs (this file, the
  `snMC-seq_*` log, etc.) so future agents and collaborators can pull them. (Ties into the
  still-open snMC-seq TODO of git-tracking `workspace/` while keeping the logs out of any nested
  repo — decide whether the logs live in their own GitHub repo, which this task now favors.)
- Build `.bai` for the 30 BAMs (samtools env already available).
- `samtools merge` per donor → 3 pseudo-bulk snM3C-seq BAMs (mirror the snMC-seq Session-6
  layout, into a `merged-BAM/` folder), then bsgenova + naive SNP calling vs `hg38_chrL.fa`.
- Phase with HapCUT2 — **this time add `extractHAIRS --hic 1`** to exploit the 3C long-range
  contacts (the whole reason for switching to snM3C-seq).

---

## Session 3 (2026-08-18) — snM3C-seq PHASING PIPELINE on the 30 downloaded BAMs

**Task (from user):** run the full phasing pipeline on the 30 already-downloaded snM3C-seq CX45
BAMs (10 per donor: H1930001, H1930002, H1930004), mirroring how the snMC-seq experiment created
folders and named files. Steps, one after another:
1. **Index** — build `.bai` for all 30 BAMs (`samtools index`).
2. **Merge per donor** — `samtools merge` the 10 cells/donor → 3 pseudo-bulk BAMs into a new
   `Science-snM3C-seq/merged-BAM/` folder, named `<donor>_CX45_snM3Cseq_10cells_merged.bam`
   (mirrors snMC's `<donor>_CX45_snMCseq_10cells_merged.bam`).
3. **SNP calling — BOTH callers** vs standard `hg38-reference/hg38_chrL.fa` (the S6 reference
   decision carries over: call germline variants against standard hg38+chrL, NOT the per-donor
   SNP-substituted refs). Reuse the shared reference from the snMC experiment folder.
   - **bsgenova** (`workspace/tools/bsgenova`): `bsextractor.py … | bsgenova.py -P 16`.
   - **naive** (`workspace/Python-scripts/bisulfite_aware_naive_SNP_caller.py`, `-P 16`).
   → unphased het callsets into `bisulfite-aware-SNPs/{from-bsgenova,from-naive}/`.
4. **HapCUT2 preprocessing** (reuse `HapCUT2_preprocessing.sh` unchanged — 4 steps: inject
   `##contig` from `hg38_chrL.fa.fai`; split multiallelics; keep het + drop chrM/chrY/chrL;
   emit plain `.vcf`) → `bisulfite-aware-SNPs/{from-bsgenova,from-naive}-HapCUT2_preprocessed/`.
5. **HapCUT2 phasing WITH Hi-C mode** — the whole point of switching to snM3C-seq:
   **`extractHAIRS --hic 1`** to exploit long-range 3C contacts that bridge het sites short reads
   never co-cover. Phased outputs → `phased-from-HapCUT2/{from-bsgenova,from-naive}/`. Compare
   phased-SNP count / block sizes against the snMC-seq baseline (bsgenova 359/129/382,
   naive 215/106/250 phased per donor).

Runs IN PARALLEL with the snMC-seq 1000-BAM download (see snMC-seq log), in its own
disconnect-safe tmux session `snm3c-phase`.

### Folder structure created (mirrors snMC-seq_SNP_phasing_experiment/)
Under `data/snM3C-seq_SNP_phasing_experiment/`:
- `Science-snM3C-seq/{H1930001,H1930002,H1930004}/` — the 30 BAMs (exist) + `.bai` (this session).
- `Science-snM3C-seq/merged-BAM/` — 3 pseudo-bulk merged BAMs + `.bai` (this session).
- `bisulfite-aware-SNPs/{from-bsgenova,from-naive}/` — raw unphased callsets.
- `bisulfite-aware-SNPs/{from-bsgenova,from-naive}-HapCUT2_preprocessed/` — preprocessed VCFs.
- `phased-from-HapCUT2/{from-bsgenova,from-naive}/` — phased blocks + `.blocks.phased.VCF`.
- Reference reused from `../snMC-seq_SNP_phasing_experiment/hg38-reference/hg38_chrL.fa` (same
  456-contig hg38+chrL; snM3C BAMs share the identical `@SQ` — to be re-confirmed this session).

### Key differences from the snMC-seq run
- BAM basenames end `.3C.sorted.bam` (not `.final.bam`); tar layout was `<cell>.3C.sorted/…`.
- Merged names use `snM3Cseq` not `snMCseq`.
- **`extractHAIRS --hic 1`** at the phasing step (snMC used standard mode). Requires extending
  `HapCUT2_run.sh` with an extractHAIRS-passthrough flag (added this session) OR a Hi-C variant.
- 3C BAMs are ~2× larger (Hi-C contact reads) → merges + bsextractor scans take longer.

### Env (ephemeral node — rebuild each session)
`bsgenova` (01), `htslib-tools` (02), `hapcut2` (03) via `environments/0{1,2,3}_*.sh`
(+ `samtools` 1.9 for index/merge/quickcheck). Channel fix in `_common.sh`.

### CHECKPOINTS
- [ ] envs rebuilt (bsgenova, hapcut2, htslib-tools, samtools)
- [ ] 30 .bai built
- [ ] 3 per-donor merged BAMs + .bai in merged-BAM/
- [ ] bsgenova callsets (3) in from-bsgenova/
- [ ] naive callsets (3) in from-naive/
- [ ] preprocessed VCFs (6) in *-HapCUT2_preprocessed/
- [ ] HapCUT2 --hic 1 phased outputs (6) in phased-from-HapCUT2/

### CHECKPOINT 1 (2026-08-18 ~17:03) — phasing pipeline LAUNCHED & running
- Orchestrator `run_snM3C-seq_phasing.sh` launched in disconnect-safe tmux session
  **`snm3c-phase`** (via `nohup nice -n 10`; survives client disconnect and tmux death).
- Envs verified present (bsgenova py3.11.15/numpy2.4.6/pysam0.24.0/samtools1.24; htslib-tools
  bcftools1.14; hapcut2 extractHAIRS+HAPCUT2 1.3.4). Now in **Stage 0** staging reference +
  30 BAMs to /tmp.
- Runs all 6 stages one-after-another (index → merge → bsgenova → naive → preprocess →
  HapCUT2 `--hic 1`), 3 donors concurrent within a stage. Each stage is resumable via
  `/tmp/snm3c_phase/.stage_<N>_done` markers; verified outputs land in the mirror folders under
  `data/snM3C-seq_SNP_phasing_experiment/`.
- **`HapCUT2_run.sh` extended** with an extractHAIRS-passthrough flag `-x` (backward-compatible;
  snMC behavior unchanged when empty). Stage 6 calls it with `-x "--hic 1"` — the 3C long-range
  contacts that are the whole point of switching to snM3C-seq.
- Monitoring: `/tmp/snm3c_phase/progress.log` (+ per-step `*_$donor.log`) and
  `/tmp/snm3c_phase/runner.log`. Attach: `tmux attach -t snm3c-phase`.
- Results (block/phased-SNP counts per donor per caller) will be recorded here as each stage
  completes, for comparison against the snMC-seq baseline.

### CHECKPOINT 2 (2026-08-18 ~17:39) — v1 pipeline FAILED (queryname sort); fixed in v2

**Important gotcha discovered:** the NeMO snM3C-seq `<cell>.3C.sorted.bam` files are
**`@HD SO:queryname` (name-sorted)** — the Hi-C convention — NOT coordinate-sorted, despite
`.sorted` in the name. (`@SQ` = 456 contigs, matches `hg38_chrL.fa` exactly; only sort order is
the issue.)

**v1 cascade failure (first run, ~17:00–17:28):**
- Stage 1 `samtools index` → `[E::hts_idx_push] Unsorted positions …` → 0 `.bai` built.
- Stage 2 `samtools merge && samtools index` → merge produced a file but the index step failed
  (name-sorted) → `MERGE_FAIL`, nothing persisted to `merged-BAM/`.
- Stages 3/4 (bsgenova/naive) ran on the **missing** merged BAMs but the shell pipe's exit code
  was the last command's, so they falsely logged OK while emitting **0-record** callsets.
- Stage 5 preprocessed the empty callsets; Stage 6 HapCUT2 failed on all 6.
- All bogus 0-record outputs were cleaned from the mount.

**v2 fix (`run_snM3C-seq_phasing.sh` rewritten, relaunched 17:39 in tmux `snm3c-phase`):**
- **New Stage 1: coordinate-sort + index each of the 30 per-cell BAMs** (`samtools sort` →
  `<cell>.3C.coordsort.bam` + `.bai`; P=6 × `-@4`). Coordinate order is also exactly what
  `extractHAIRS --hic 1` requires. This is the real "build .bai for the 30" step.
- Stage 2 merges the **coordinate-sorted** per-cell BAMs per donor → coord-sorted pseudo-bulk
  BAM (+`.bai`), with a read-count check (merged == Σ per-cell) before persisting to `merged-BAM/`.
- **Hard input guards** added to Stages 3/4/6: a caller/phaser is skipped with an explicit
  `*_FAIL` (never a false OK) unless its merged BAM+`.bai` exists and is non-empty; record counts
  are logged after each caller.
- Per-cell coord-sorted BAMs live in `/tmp` as intermediates; persisted analysis artifacts are the
  3 coord-sorted merged BAMs (+`.bai`) plus callsets/preprocessed/phased. (The authentic
  queryname-sorted downloads in the donor folders are left untouched.)
- Monitoring: `/tmp/snm3c_phase/progress.log`, `runner_v2.log`, per-step `*_$donor.log`.
- Results (block/phased-SNP counts per donor per caller, with `--hic 1`) will be recorded here on
  completion and compared to the snMC-seq baseline.

### CHECKPOINT 3 (2026-08-19 15:31) — v2 died mid-Stage-2 on FUSE; v3 launched, resumed at Stage 3

**What actually happened to the v2 run (18:22 Aug 18):** it did NOT complete. It finished
Stages 1–2 correctly and then the *interpreter itself* died:
`run_snM3C-seq_phasing.sh: error reading input file: Remote I/O error`.
Two independent iRODS/FUSE failure modes were at work:

1. **The running script lived on the mount.** `bash` re-reads a script from disk as it
   executes it, so one transient mount read error killed the whole pipeline — after ~43 min of
   sorting and merging. (Same exposure applied to every `bash "$PREP"` / `bash "$HRUN"` call.)
2. **`cp` of the ~6 GB merged BAMs to the mount raised `Remote I/O error`.** Result on the mount:
   H1930001 landed intact (6304569554 B), H1930002 never written, H1930004 **truncated**
   (5078253568 of 5926788102 B). Worse, `mark 2` was gated on *mount* file count, so a mount
   write failure would have forced a full ~30-min re-merge on the next run.

**State found intact in `/tmp` (node up 14 d, `/opt/conda` + `/tmp` both survived):**
- 30/30 per-cell `.3C.coordsort.bam` + `.bai` in `/tmp/snm3c_phase/sorted/`
- 3/3 merged pseudo-bulk BAMs + `.bai` in `/tmp/snm3c_phase/merged/`, all `samtools quickcheck`
  clean, read counts already verified equal to Σ per-cell (98270879 / 103086888 / 95040464)
- all 5 conda envs (`bsgenova`, `hapcut2`, `htslib-tools`, `samtools`, `tmux`) still built
So **Stage 1 and Stage 2 compute never needed to be redone** — only mount publication did.

**v3 changes (`run_snM3C-seq_phasing.sh`, 281 lines, published to `shell-scripts/`):**
- **Code runs from `/tmp/snm3c_local/`, never the mount.** Local copies of the orchestrator,
  `HapCUT2_preprocessing.sh`, `HapCUT2_run.sh`, the naive caller and `bsgenova/{bsextractor,bsgenova}.py`.
  The mount is now only an input source and a publish target — it is off the code path entirely.
- **COMPUTE decoupled from PUBLISH.** Stage markers depend only on `/tmp` artifacts. `cpv()` retries
  each mount copy 3× with backoff and is **non-fatal** (`CP_FAIL` logged); a new **Stage 7**
  re-publishes everything still missing at the end and prints the mount inventory.
- **Stages 5/6 now read local inputs and write local outputs** (`/tmp/snm3c_phase/prep/{bsgenova,naive}`,
  `phased/{bsgenova,naive}`), so a mount hiccup can no longer break phasing. Published afterwards.
- **Idempotent per-item skips** everywhere (`MERGE_SKIP`/`BSG_SKIP`/`NAIVE_SKIP`/prep/phase guards),
  so a re-launch never redoes finished work; Stage 2 quickchecks the local BAMs instead of re-merging.
- Callers bumped to `-P 24` each (3 donors concurrent = 72 of 128 cores; 503 GB RAM).
- Truncated `H1930004` BAM removed from `merged-BAM/` on the mount (intact copy lives in `/tmp`).

**v3 run (launched 15:31:31, `setsid nohup nice -n 10`, pid 205206):**
- envs verified in 3 s; Stage 0 re-staged nothing (ref + 30 raw BAMs already in `/tmp`)
- Stage 1 `already done, skip`; **Stage 2 `MERGE_SKIP ×3` → 3/3 valid → marked** (no re-merge)
- **Stage 3 bsgenova started 15:31:53** on all three real merged BAMs — i.e. past the point where
  v2 died. Monitoring `/tmp/snm3c_phase/runner_v3.log` (+ `progress.log`, per-step `*_$donor.log`).
- Still to come: Stage 4 naive, Stage 5 preprocessing (6), Stage 6 HapCUT2 `--hic 1` (6),
  Stage 7 publish. Block / phased-SNP counts will be recorded here and compared with the snMC-seq
  baseline (bsgenova 359/129/382, naive 215/106/250 phased per donor).

### CHECKPOINTS
- [x] envs rebuilt (bsgenova, hapcut2, htslib-tools, samtools)
- [x] 30 per-cell coordinate-sorted BAMs + .bai
- [x] 3 per-donor merged coord-sorted BAMs + .bai (in /tmp; 1/3 on mount, Stage 7 to re-publish)
- [ ] bsgenova callsets (3)
- [ ] naive callsets (3)
- [ ] preprocessed VCFs (6)
- [ ] HapCUT2 --hic 1 phased outputs (6)

### CHECKPOINT 4 (2026-08-19 15:50) — v3 past the v2 failure point; Stage 3 producing real callsets

- v3 launched **15:31:31** (pid 205206, `setsid nohup nice -n 10`, log
  `/tmp/snm3c_phase/runner_v3.log`).
- **Stages 1–2 correctly resumed, not recomputed.** Envs verified in 3 s; Stage 0 re-staged nothing
  (reference + 30 raw BAMs already in `/tmp`); Stage 1 `already done, skip`; Stage 2 logged
  `MERGE_SKIP ×3` after `samtools quickcheck` on each local merged BAM (6304569554 / 6227682158 /
  5926788102 bytes) → `3/3` → marked. The ~43 min of sorting + merging from the v2 run was reused
  in full; nothing was re-merged.
- **Stage 3 bsgenova started 15:31:53** — i.e. v3 got past the exact point where v2 died — with all
  three `bsextractor.py | bsgenova.py -P 24` pipes confirmed running on the real merged BAMs
  (82 worker processes, load ~8).
- By 15:50 all three donors are emitting **non-empty, growing** `.vcf.gz` + `.snv.gz` into
  `/tmp/snm3c_phase/out/` (H1930001 135 KB, H1930002 410 KB, H1930004 300 KB and climbing) —
  the decisive contrast with the v1 run, which produced 0-record callsets from missing BAMs.
- Still to come: Stage 4 naive, Stage 5 preprocessing (6), Stage 6 HapCUT2 `--hic 1` (6), Stage 7
  publish. Block / phased-SNP counts go here on completion, against the snMC-seq 10-cell baseline
  (bsgenova 359/129/382, naive 215/106/250 phased per donor).
- The 2 merged BAMs that failed to reach the mount in v2 are **not** blocking anything: they are
  valid in `/tmp`, and v3's Stage 7 retries publication at the end.

### NOTE (2026-08-19 16:01) — mount verification caveat affecting this experiment's published copies

While mopping up the snMC-seq download, the iRODS mount was measured directly and found to
**return a wrong size AND wrong bytes for a file read immediately after `cp` returns** — a
205163551 B BAM read back as 209000510 B with a mismatching md5, then correct on the next
attempt. **Performing a full read is itself what settles the object.** Full write-up in the
snMC-seq log, CHECKPOINT 3.

Consequences for this experiment:
- Every `WARN copy-verify size lag` in the v2 run above is explained by this — including the two
  merged-BAM copies. A passing size check is not proof of a good copy on this mount, and a failing
  one is not proof of a bad copy.
- **v3's `cpv()` verifies published copies by SIZE only.** v3 was already running and a running
  bash script must not be edited — that is precisely what killed v2 — so it was left alone.
- Instead, **`shell-scripts/verify_published_artifacts.sh` (new)** md5-compares every artifact on
  the mount against the `/tmp` original the pipeline produced and re-copies mismatches. Run it
  after v3 finishes: `bash verify_published_artifacts.sh snm3c`. That, not v3's own size check, is
  the authoritative content check for the snM3C-seq deliverables.

### CHECKPOINT 5 (2026-08-19 16:25) — Stage 3 progress measured: ~24 % of the genome at 53 min

bsgenova emits variants in reference-`.fai` contig order, so the last emitted `CHROM POS`
cumulated over the contig lengths is a cheap and honest progress meter (no need to guess from
file sizes). At 53 min into Stage 3:

| donor | last position | genome covered | records so far |
|---|---|---|---|
| H1930001 | chr13:111,384,115 | 763 / 3209 Mb = **23.8 %** | 47,905 |
| H1930002 | chr14:20,581,955  | 786 / 3209 Mb = **24.5 %** | 76,612 |
| H1930004 | chr14:52          | 766 / 3209 Mb = **23.9 %** | 61,104 |

All three donors are advancing in step (they run concurrently at `-P 24` each). Linear
extrapolation gives **~3.7 h for Stage 3, finishing ~19:15 UTC**; the `.fai` order is roughly
lexicographic (chr1, chr10, chr11, chr12, chr13, chr14, chr15…), so the remaining share still
contains the large chr2–chr9, and the many small `*_random`/decoy contigs at the tail are cheap.
Stages 4–7 (naive caller, preprocessing, HapCUT2 `--hic 1`, publish) follow automatically.

Note the contrast with the v1 run, which "finished" bsgenova in ~20 min — because it was scanning
missing BAMs and emitting 0 records. A caller that takes hours on a 6 GB pseudo-bulk is the
expected behaviour.

### CHECKPOINT 6 (2026-08-19 16:42) — INCIDENT: my own watcher published a partial callset; fixed

**Self-inflicted, caught within 5 minutes, one file affected.** Recording it in full because the
mistake is a trap anyone re-running these pipelines will hit.

At 16:33 a helper (`watch_snm3c_then_verify.sh` v1) was launched to md5-verify v3's published
artifacts after it finished. It tested completion with:

```bash
grep -q 'ALL_STAGES_DONE' /tmp/snm3c_phase/progress.log
```

**`progress.log` is append-only and SHARED by v1, v2 and v3.** The **v1** run wrote
`=== PIPELINE COMPLETE ===` + `ALL_STAGES_DONE` at **2026-08-18 17:28:35** — immediately after
`Stage 6 phased VCFs=0`, i.e. the bogus run that phased nothing. That line never goes away, so the
grep was a **permanent false positive**. The watcher declared v3 "done" 1 minute after arming, ran
the verifier, and the verifier's job is to re-copy `/tmp` → mount on mismatch — so it began
**publishing v3's in-progress callsets as final results.**

Damage: exactly one file — `from-bsgenova/H1930001_…merged.vcf.gz` at **898840 bytes**, a snapshot
of a scan that was only ~24 % through the genome. It was removed, and the object name was probed
afterwards and confirmed still writable, so v3's real Stage 3 publish is unaffected. v3 itself was
never touched (82 workers throughout, callsets still growing). All other output dirs were untouched
because the verifier spent its first 5 minutes on the merged-BAM checks.

**Root cause in one line: a stage-completion test that reads a log shared across runs, instead of
per-run state.**

**Fixes — completion is now defined by per-run `/tmp` stage markers, never by a log line:**
- `watch_snm3c_then_verify.sh` v2 — completion requires **(a)** the pipeline process to be gone
  **and (b)** `.stage_{3,4,5,6}_done` to all exist. If the process exits without them it logs
  `INCOMPLETE RUN`, prints the last progress lines and **publishes nothing**.
- `verify_published_artifacts.sh` — now carries a `require_complete()` **guard of its own**, so it
  refuses even when invoked by hand. Self-tested just now against the live in-flight run:
  `REFUSING to verify snM3C-seq: incomplete run, missing stage markers: 3 4 5 6`. `FORCE=1`
  overrides deliberately.
- `supervise_snMC_phasing.sh` — had the identical latent flaw (its `progress.log` happens to be
  clean today only because it was created fresh). Rewritten to use the same `complete_p()`
  stage-marker test.

Stage markers are the right authority here: they live in the per-run `/tmp` working tree and each
is written only after its stage verified its own outputs, so they cannot be inherited from a
previous run the way a log line can.

### CHECKPOINT 7 (2026-08-19 18:59) — **Stage 3 COMPLETE: 3/3 bsgenova callsets, ~150× the 10-cell snMC-seq yield**

Stage 3 ran 15:31:53 → 18:59:27 (**3 h 28 min**, against the 3.7 h projected at the 53-min mark —
the `.fai`-position progress meter was accurate).

| donor | bsgenova records | published to mount |
|---|---|---|
| H1930001 | **186,988** | `.vcf.gz` 2.68 MB + `.snv.gz` 3.33 MB |
| H1930002 | **302,301** | `.vcf.gz` 4.35 MB + `.snv.gz` 5.35 MB |
| H1930004 | **230,426** | (publishing) |

`Stage 3 bsgenova callsets (local, non-empty) = 3/3`.

**Against the snMC-seq 10-cell baseline — same donors, same caller, same reference — with read
depth included so the comparison is honest:**

| dataset | reads in pseudo-bulk | bsgenova records | records/dataset ratio | read ratio |
|---|---|---|---|---|
| snMC 10-cell H1930001 | 30,441,510 | 1,166 | — | — |
| snMC 10-cell H1930002 | 19,211,865 | 421 | — | — |
| snMC 10-cell H1930004 | 37,517,754 | 1,712 | — | — |
| snM3C 10-cell H1930001 | 98,270,879 | **186,988** | **160×** | 3.2× |
| snM3C 10-cell H1930002 | 103,086,888 | **302,301** | **718×** | 5.4× |
| snM3C 10-cell H1930004 | 95,040,464 | **230,426** | **135×** | 2.5× |

**The variant gain (135–718×) vastly exceeds the read gain (2.5–5.4×).** A plausible reading: at
~20–37 M reads over a 3.2 Gb reference the genome-wide depth is well under 1×, so almost no site has
the ≥2 reads needed to call a heterozygote; at ~100 M reads (~1–3×) a large fraction of sites cross
that threshold. Variant yield is strongly non-linear in depth in this regime, so a 3× depth increase
producing a >100× callset increase is not implausible — **but the magnitude should not be taken at
face value yet.** Raw callset size is not the deliverable, and part of this could be low-confidence
or spurious calls that the het-filtering in Stage 5 will strip.

**The number that actually matters is the phased-SNP count after Stage 5 het filtering + Stage 6
HapCUT2 `--hic 1`**, against the snMC-seq baseline of **359 / 129 / 382** phased SNPs per donor.
That is the comparison this experiment was set up to make; it will be recorded here on completion.

Stage 4 (naive caller) now running on the same three merged BAMs, then Stage 5 preprocessing,
Stage 6 phasing with `--hic 1`, Stage 7 publish. The marker-gated watcher will md5-verify all
published artifacts once the stage 3–6 markers are all present.

### CHECKPOINT 8 (2026-08-19 20:00) — Stage 6 failed on a doc error in HapCUT2_run.sh; fixed and re-running

**Stage 6 failed all 6 runs in 68 seconds.** Not a data problem — extractHAIRS worked perfectly and
extracted **372,626 / 505,476 / 451,997** (bsgenova) and **372,850 / 471,878 / 423,770** (naive)
Hi-C fragments, confirming the 3C long-range contacts are there. HAPCUT2 then died on every one:

```
ERROR: Invalid fragment file, too many fields at line 0.
If the file is the new HiC-related format (from extractHAIRS --hic option), be sure to
use --hic, --hic_htrans_file, or --nf HapCUT2 options.
```

**Root cause: a false claim in `HapCUT2_run.sh`'s own header**, which said "When Hi-C fragments are
used, HAPCUT2 auto-detects them from the fragment format." **It does not.** HAPCUT2's usage text is
explicit: *"(1) When running extractHAIRS, must use --hic 1 … (2) When running HapCUT2, use --hic 1
if h-trans probabilities are unknown."* **Both** stages need the flag; the script passed it only to
extractHAIRS. `--hic 1` on HAPCUT2 also enables h-trans error modelling from the data, which is
precisely what makes Hi-C phasing work.

**Fixes in `HapCUT2_run.sh`:**
- `--hic 1` is now passed to HAPCUT2 automatically whenever `-x` contains `--hic`, with a log line
  saying so; new `-y` flag passes arbitrary extra args to HAPCUT2 and overrides the auto-derivation.
- The incorrect header claim is replaced by the actual error text and what it cost, so this is not
  re-derived by anyone.

Stages 1–5 markers were all set, so the re-run **resumed directly at Stage 6** — none of the 3 h 28
min of genotyping was repeated. Re-running since 19:59:51.

**Note — v3's Stage 7 published successfully before exiting:**
`mount: merged=3/3 bsg_callsets=3/3 naive_callsets=3/3 prep=6/6 phased=0/6`.
The two merged BAMs that failed to copy in v2 finally landed, so `merged-BAM/` is complete.

### Trap worth remembering: `pgrep -f` matches shells that merely MENTION the pattern

The Stage 6 re-run driver sat blocked for 27 minutes after v3 had already exited. Cause: a **stale
wrapper shell from the 2026-08-18 v2 session (PID 77148, started 17:39:56 Aug 18)** whose command
line literally contained `while pgrep -f run_snM3C-seq_phasing.sh ...`. That text matches the
pattern, so `pgrep -f 'run_snM3C-seq_phasing'` reported a live pipeline for **26 hours**. The
watcher had the same flaw and would also never have fired.

This is the same trap that repeatedly killed the interactive shell today (`pkill -f <pattern>`
matching its own wrapper, exit 144). Fix, applied to `watch_snm3c_then_verify.sh`,
`supervise_snMC_phasing.sh` and `chain_snMC_download_then_phase.sh`:

```bash
pipeline_running(){ pgrep -af '<name>' 2>/dev/null \
  | grep -qE '^[0-9]+ (bash|/bin/bash) /tmp/snm3c_local/<script>\.sh'; }
```

anchoring the match to the real invocation and excluding `-c` wrapper shells. Behaviour-tested
against a decoy shell that mentions the script name: correctly reports NOT running. The stale
PID 77148 was killed.

## RESULTS (2026-08-19 20:41) — snM3C-seq phasing COMPLETE, all 6 runs

`Stage 6 phased VCFs = 6/6`; `mount: merged=3/3 bsg_callsets=3/3 naive_callsets=3/3 prep=6/6
phased=6/6`; `PIPELINE COMPLETE`. All stage markers 1–6 set. Total pipeline time for the productive
run: Stage 3 3 h 28 min + Stage 4 20 min + Stage 5 10 min + Stage 6 41 min.

### Headline: 85–542× more phased SNPs than the snMC-seq 10-cell baseline

Same donors, same callers, same reference, same cell count (10). Only the assay and `--hic 1` differ.

| caller / donor | snMC-seq 10-cell | **snM3C-seq 10-cell (`--hic 1`)** | ratio |
|---|---|---|---|
| bsgenova / H1930001 | 359 | **42,081** | **117×** |
| bsgenova / H1930002 | 129 | **69,952** | **542×** |
| bsgenova / H1930004 | 382 | **54,076** | **142×** |
| naive / H1930001 | 215 | **18,289** | **85×** |
| naive / H1930002 | 106 | **27,884** | **263×** |
| naive / H1930004 | 250 | **21,456** | **86×** |

### Where the gain actually comes from — and where it does NOT

**It comes from het-site discovery, not from better linkage per site.** The snM3C pseudo-bulks carry
2.5–5.4× more reads (~3.0–3.2× genome depth vs ~0.6–1.2× for snMC-seq 10-cell), which pushed
het-site counts from 271–1,211 to 60,099–262,266. Phasing then converts a similar-or-lower
*fraction* of them:

| | phased / het sites | rate |
|---|---|---|
| snMC bsgenova H1930001 | 359 / 792 | 45.3 % |
| snM3C bsgenova H1930001 | 42,081 / 163,894 | **25.7 %** |
| snMC naive H1930001 | 215 / 419 | 51.3 % |
| snM3C naive H1930001 | 18,289 / 60,099 | **30.4 %** |

The baseline's higher *rate* is at least partly a selection effect — with only a few hundred het
sites called, those sit in the deepest, easiest-to-link regions. But it means **`--hic 1` did not
raise the per-site phasing rate above the short-read baseline.**

**Block contiguity did not improve on average, and that is the honest disappointment here:**

| | SNPs per block | largest block (SNPs) | largest block span |
|---|---|---|---|
| snMC bsgenova H1930001 | **3.29** | 18 | — |
| snM3C bsgenova H1930001 | **2.26** | **139** | 44–51 kb |
| snMC naive H1930001 | **3.52** | — | — |
| snM3C naive H1930001 | **2.89** | — | 50.8 kb |

Mean SNPs/block is *slightly lower* than the sparse baseline (2.23–2.28 vs 2.89–3.29 for bsgenova).
The long-range benefit shows up only in the tail: the largest block grew from **18 to 139 SNPs**
(7.7×). So the distribution gained a long tail while remaining dominated by 2-SNP blocks.

### Likely cause, and the concrete follow-up it implies

Largest block spans are only **44–51 kb**. Hi-C contacts routinely span megabases, so the 3C
long-range linkage is barely being exploited. The most probable reason is stated in HapCUT2's own
documentation:

> "Using `--hic_htrans_file` is faster than `--hic` and may yield better results at low read
> coverage (<30x)."

`--hic 1` *estimates h-trans error probabilities from the data*, and it needs matepairs per
insert-size window to do so (`--htrans_read_lowbound` defaults to **500** per window). At the
**~3× depth** measured here — an order of magnitude below the <30× regime the docs call out — those
per-window estimates are poorly constrained, which would cap block extension regardless of how many
contacts exist. extractHAIRS did find **372k–505k fragments**, so the contacts are present; the
limitation appears to be h-trans estimation, not fragment availability.

**Recommended next experiment:** re-run Stage 6 with `--hic_htrans_file` supplying h-trans
probabilities per insert-size bin instead of estimating them, and/or raise depth by using more than
10 cells per donor. `HapCUT2_run.sh` now has a `-y` passthrough, so this is a one-flag change:
`-y "--hic_htrans_file <file>"`. Comparing block N50 between the two would isolate whether h-trans
estimation is the true bottleneck.

**Bottom line:** switching to snM3C-seq delivered a 85–542× increase in phased SNPs, which is a
decisive win on yield and answers the question the experiment was set up to ask. But the specific
promise of 3C — long haplotype blocks from long-range contacts — is **not yet realised**; blocks are
short and slightly more fragmented than the short-read baseline, and the evidence points at
low-coverage h-trans estimation as the reason.

### CHECKPOINT 9 (2026-08-19 20:57) — independent md5 verification of every published artifact: CLEAN

The marker-gated watcher correctly detected genuine completion this time
(`all stage markers present -- run genuinely complete`, 20:42:43) and ran the md5 verification:

```
SUMMARY  verified_ok=72  repaired=0  still_bad=0
```

**All 72 artifacts on the mount are byte-identical to the `/tmp` originals** — the 3 merged BAMs and
`.bai`, 6 raw callsets (`.vcf.gz` + `.snv.gz`), 6×4 preprocessed files, and 6×5 phased files
(`.fragments`, `.blocks`, `.blocks.phased.VCF`, and both logs). Nothing needed repair.

This closes the concern raised in CHECKPOINT 6 / snMC CHECKPOINT 3: v3's `cpv()` verifies published
copies by **size only**, and this mount can report a plausible size for wrong bytes — but in the
event every copy was in fact correct. The size-only check was insufficient *in principle*; the
independent md5 pass is what makes the deliverables trustworthy *in fact*.

**snM3C-seq experiment status: COMPLETE and verified.** Deliverables on the mount under
`data/snM3C-seq_SNP_phasing_experiment/`:
- `Science-snM3C-seq/merged-BAM/` — 3 coordinate-sorted pseudo-bulk BAMs + `.bai`
- `bisulfite-aware-SNPs/from-{bsgenova,naive}/` — 6 raw callsets
- `bisulfite-aware-SNPs/from-{bsgenova,naive}-HapCUT2_preprocessed/` — 6 het-filtered VCFs
- `phased-from-HapCUT2/from-{bsgenova,naive}/` — 6 phased blocks + `.blocks.phased.VCF`

### CHECKPOINT 10 (2026-08-20) — depth-matched snMC-seq vs snM3C-seq at 10 cells, and **the real reason `--hic 1` produced no long blocks**

Prompted by the question "how do the two assays compare at the same cell count?", I recomputed every
number directly from the published artifacts rather than reusing the CHECKPOINT 8 summary. Two of
those numbers turned out to be wrong, and chasing the second one found the actual bottleneck.

**Depth first — the comparison is NOT depth-matched, and the gap is bigger than assumed:**

| donor | snMC-seq reads | snM3C-seq reads | snMC bases | snM3C bases | depth ratio (bases) |
|---|---:|---:|---:|---:|---:|
| H1930001 | 30,441,510 | 98,270,879 | 3.81 Gbp (1.23×) | 7.01 Gbp (2.26×) | 1.84× |
| H1930002 | 19,211,865 | 103,086,888 | 2.39 Gbp (0.77×) | 7.68 Gbp (2.48×) | 3.21× |
| H1930004 | 37,517,754 | 95,040,464 | 4.66 Gbp (1.50×) | 7.01 Gbp (2.26×) | 1.51× |

Mean read length differs and matters: **snMC-seq 124.1–125.1 bp vs snM3C-seq 71.3–74.5 bp**
(sampled per donor from each merged BAM), because the snM3C reads are chimeric-split fragments.
Read *count* ratios (3.2–5.4×) therefore overstate the depth advantage; base ratios are 1.5–3.2×,
i.e. ~0.8–1.5× genome coverage for snMC-seq against ~2.3–2.5× for snM3C-seq.

**Full side-by-side, 10 cells per donor, same reference and callers:**

| caller | donor | assay | het into HapCUT2 | blocks | phased SNPs | rate | SNPs/block | phased per M reads |
|---|---|---|---:|---:|---:|---:|---:|---:|
| bsgenova | H1930001 | snMC | 792 | 109 | 359 | 45.3 % | 3.29 | 11.8 |
| bsgenova | H1930001 | **snM3C** | 163,894 | 18,585 | **42,081** | 25.7 % | 2.26 | **428** |
| bsgenova | H1930002 | snMC | 271 | 44 | 129 | 47.6 % | 2.93 | 6.7 |
| bsgenova | H1930002 | **snM3C** | 262,266 | 31,378 | **69,952** | 26.7 % | 2.23 | **679** |
| bsgenova | H1930004 | snMC | 1,211 | 132 | 382 | 31.5 % | 2.89 | 10.2 |
| bsgenova | H1930004 | **snM3C** | 198,106 | 23,761 | **54,076** | 27.3 % | 2.28 | **569** |
| naive | H1930001 | snMC | 419 | 61 | 215 | 51.3 % | 3.52 | 7.1 |
| naive | H1930001 | **snM3C** | 60,099 | 6,318 | **18,289** | 30.4 % | 2.89 | **186** |
| naive | H1930002 | snMC | 180 | 31 | 106 | 58.9 % | 3.42 | 5.5 |
| naive | H1930002 | **snM3C** | 99,905 | 9,997 | **27,884** | 27.9 % | 2.79 | **270** |
| naive | H1930004 | snMC | 457 | 72 | 250 | 54.7 % | 3.47 | 6.7 |
| naive | H1930004 | **snM3C** | 75,502 | 7,522 | **21,456** | 28.4 % | 2.85 | **226** |

**Is the snM3C win just depth?** No — and this is quantifiable. The snMC-seq depth series measured
on this project (1.23× → 58.5× coverage gave 359 → 976,575 phased SNPs, H1930001/bsgenova) implies
yield scales as **depth^2.05**. Applying that exponent to each donor's base ratio predicts how much
of the cross-assay gap depth alone explains:

| caller | donor | observed ratio | depth-predicted | **residual** |
|---|---|---:|---:|---:|
| bsgenova | H1930001 | 117× | 3.5× | **34×** |
| bsgenova | H1930002 | 542× | 10.9× | **50×** |
| bsgenova | H1930004 | 142× | 2.3× | **61×** |
| naive | H1930001 | 85× | 3.5× | **24×** |
| naive | H1930002 | 263× | 10.9× | **24×** |
| naive | H1930004 | 86× | 2.3× | **37×** |

So **24–61× of the snM3C advantage is not depth** — it is assay quality (library complexity, more
informative fragments per base). Caveat: the exponent was fitted within snMC-seq and extrapolated
across assays, so treat these as order-of-magnitude, not precise.

**CORRECTION to CHECKPOINT 8.** That section reported the largest snM3C block spanning **44–51 kb**.
That is wrong. Measured directly, the largest block (bsgenova/H1930001, 139 SNPs) spans **422 bp** —
the `SPAN:` field and the first/last variant positions agree exactly. Across all six snM3C datasets:

| | mean block span | max block span | blocks > 1 kb |
|---|---:|---:|---:|
| snM3C bsgenova (3 donors) | 29–30 bp | 897–974 bp | **0** |
| snM3C naive (3 donors) | 35–36 bp | 913–1,023 bp | **2 of 27,837** |
| snMC bsgenova (3 donors) | 70–82 bp | 504 bp | 0 |

**Not one snM3C block reaches 1 kb.** Blocks are *shorter* than the snMC-seq baseline's, exactly in
proportion to the shorter reads (70 bp vs 125 bp). Every block is intra-read.

**Root cause — there were never any Hi-C contacts in the fragment file.** extractHAIRS reports
**`PE-fragments 0`** for all six snM3C runs. The merged BAMs contain no paired reads at all: sampled
alignments carry FLAG 0/16 with empty mate fields, and read names look like

```
A00280:553:HMMVVDSX3:3:1675:14696:29543_1:N:0:CATGGTGTAC+TACAGGTCCT-r
A00280:553:HMMVVDSX3:3:1675:14696:29543_2:N:0:CATGGTGTAC+TACAGGTCCT-l
```

i.e. the upstream snm3C (YAP/cemba) pipeline splits chimeric reads and maps every piece
**single-end**, encoding the original pairing in `_1`/`_2` and `-l`/`-r` name suffixes. Mates are
therefore *different read names* to samtools and to extractHAIRS, so `--hic 1` had nothing to pair
and silently degraded to plain short-read phasing.

**The contacts are in the data — just not visible to the tool.** Stripping the suffixes and
re-linking by base name on chr21 alone (donor H1930001) recovers:

```
693,148 intra-chr21 name-linked pairs
165,220 (23.8 %) span > 1 kb
112,256 (16.2 %) span > 100 kb
 53,749 ( 7.8 %) span > 1 Mb
```

That is precisely the megabase-scale linkage HapCUT2 needs, and it was discarded before phasing.

**This supersedes the CHECKPOINT 8 hypothesis.** h-trans estimation was a reasonable guess, but it
cannot be the bottleneck: `--hic_htrans_file` only refines error probabilities for mate pairs that
exist, and there were zero. Supplying an h-trans file would have changed nothing.

**Concrete fix (high value, low cost):** rewrite the merged BAM so mates are real pairs before
Stage 6 — strip the `_[12]:N:0:...-[lr]` suffix to a common QNAME, set FLAG bits 0x1/0x40/0x80,
fill RNEXT/PNEXT via `samtools fixmate` on a queryname-sorted copy, coordinate-sort, then re-run
`extractHAIRS --hic 1`. Success criterion is unambiguous and cheap to read off: **`PE-fragments`
must be non-zero**, and block spans should jump from <1 kb into the 100 kb–Mb range. This is a
few hours on the existing 10-cell BAMs and needs no new downloads.

**Bottom line:** the 85–542× yield win over snMC-seq stands and is mostly real assay quality. But
snM3C-seq's defining advantage — long-range linkage — has still never been tested on this data. It
was not a low-coverage limitation; it was a BAM-format mismatch that made the contacts invisible.

---

## Session 4 (2026-08-20) — TWO PARALLEL TASKS LAUNCHED (read this to resume)

Two independent, disconnection-safe runs were started this session (each in its own
detached tmux session, each also `setsid nohup` so it survives tmux death AND a VSCode /
VSCode-server disconnect). **This file covers Task 1 (snM3C-seq).** Task 2 (snMC-seq
read-downsampling saturation series) is documented in `snMC-seq_SNP_phasing_experiment.md`.

### Environment rebuilt this session (ephemeral node — /tmp and conda envs were wiped)
- Node still up (16 d) but `/tmp/snm3c_phase`, `/tmp/snm3c_local`, and ALL conda envs were
  gone at session start. Rebuilt: envs `bsgenova`, `htslib-tools`, `hapcut2`, `samtools`
  (idempotent recipes in `environments/`). Reference re-staged to `/tmp/snm3c_phase/hg38_chrL.fa`.
- Code re-staged to `/tmp/snm3c_local/` (orchestrators, HapCUT2_preprocessing.sh, callers).
- **`HapCUT2_run.sh` on the mount is a STUCK/CORRUPT iRODS object** — every read returns
  `Remote I/O error` while every other file reads fine (~25 retries failed). It was
  **reconstructed from scratch** this session (same `-b/-v/-o/-n/-x/-y` interface, the
  CHECKPOINT-8 auto-`--hic 1`-to-HAPCUT2 behaviour, per-contig extractHAIRS loop) and saved to
  `/tmp/snm3c_local/shell-scripts/HapCUT2_run.sh`; it is republished to the mount when writable.
  If you find the mount copy still 0-length/unreadable, use the /tmp copy.

### TASK 1a — snM3C read-pairing fix: DEVELOPED AND VALIDATED ✅ (the CHECKPOINT-10 fix)
CHECKPOINT 10 proved `--hic 1` produced no long blocks because the merged BAMs have **zero
real pairs** (`PE-fragments 0`): the upstream YAP/cemba snm3C pipeline splits chimeric reads
at ligation junctions and maps every piece SINGLE-END, encoding pairing only in the read NAME:
```
<base>_<1|2>:N:0:<index>[-<l|m|r>]
   _1/_2  = the two ends of the PE read (the ligation-junction ends)
   -l/-m/-r = chimeric sub-segments of one arm (usually local)
```
Measured on the H1930001 10-cell merged BAM, chr21: 566,789 distinct base molecules,
395,151 with ≥2 alignments (candidate intra-chr21 contacts). Multiplicity 1–6 segments/base.

**The fix (script: `Python-scripts/repair_m3c_pairs.py`, new this session):** group alignments
by the shared `<base>`; PER CONTIG pair the two OUTERMOST segments (min-pos, max-pos → maximum
phasing span) into one proper R1/R2 pair (flags 0x1|0x40 / 0x1|0x80, mate-bits cleared first);
keep remaining/singleton segments as unpaired singles (retain base coverage); never pair across
contigs (can't phase across chromosomes). Pipeline:
```
samtools sort -n  (groups by base QNAME)  ->  repair_m3c_pairs.py  ->
samtools sort -n  ->  samtools fixmate -m  ->  samtools sort (coord)  ->  index  ->
extractHAIRS --hic 1  ->  HAPCUT2 --hic 1
```
`fixmate` fills RNEXT/PNEXT/TLEN/mate-reverse. "properly paired 0%" is EXPECTED and fine — Hi-C
inserts are huge; extractHAIRS `--hic` only needs the mate linkage, not the 0x2 proper-pair bit.

**Validation on chr21 (H1930001 10-cell, 3243 het sites), original vs repaired BAM:**

| | max block span | blocks >1 kb | blocks >100 kb | blocks >1 Mb | phased SNPs |
|---|---:|---:|---:|---:|---:|
| original (single-end) | 281 bp | 0 | 0 | 0 | 1,097 |
| **repaired (paired)** | **27,308,980 bp** | 28 | 14 | **11** | 1,329 |

i.e. from <1 kb (matching CHECKPOINT 10) to **27.3 Mb** blocks — the megabase-scale linkage
snM3C-seq was supposed to give. The fix is confirmed; it is now baked into the orchestrator.

### TASK 1 — full plan (tmux session `m3c-task1`, driver `run_snM3C-seq_task1.sh`)
Runs these stages one after another; resumable via `/tmp/m3c_task1/.stage_*_done` markers;
compute in /tmp, publish to the mount via retrying non-fatal `cpv()`; log
`/tmp/m3c_task1/progress.log`.

- **Step A — re-phase the existing 3× 10-cell snM3C merged BAMs WITH the pairing fix.**
  Inputs: `Science-snM3C-seq/merged-BAM/<donor>_CX45_snM3Cseq_10cells_merged.bam` (all 3 donors,
  on the mount). For each: repair → extractHAIRS `--hic 1` → HAPCUT2 `--hic 1` reusing the
  existing preprocessed het VCFs in `bisulfite-aware-SNPs/from-{bsgenova,naive}-HapCUT2_preprocessed/`.
  Outputs the corrected phasing to `phased-from-HapCUT2-hicfix/from-{bsgenova,naive}/`
  (NOT overwriting the old `phased-from-HapCUT2/` so the before/after stays on record).
  Deliverable: corrected block counts / spans vs the broken CHECKPOINT-9 numbers.

- **Step B — download 90 MORE H1930001 CX45 snM3C-seq BAMs** (cells 11–100 by tar size, same
  donor ⇒ same ancestry, same region CX45; skip the 10 already downloaded). Adapted downloader
  `download_snM3C-seq_more.sh`: manifest filter `_3C_` + `_H1930001_CX45_`, `.3C.sorted.bam.tar`,
  m3C-seq NeMO subtree, md5-verify vs manifest, 3 concurrent net streams + serial mount writer.
  Land in `Science-snM3C-seq-100cell/H1930001/`.

- **Step C — 100-cell pseudo-bulk phasing pipeline for H1930001.** coordinate-sort+index all 100
  per-cell BAMs (the downloads are queryname-sorted `.3C.sorted.bam`, see CHECKPOINT 2) → merge →
  **apply the pairing fix to the merged BAM** → sharded bsgenova + naive genotyping vs
  `hg38_chrL.fa` → HapCUT2 preprocessing → HAPCUT2 `--hic 1`. Outputs under
  `data/snM3C-seq_100cell_experiment/`. **Deliverable: does 10→100 cells improve phasing
  (phased-SNP count, block N50, max span)** vs the 10-cell fixed baseline from Step A.

Baselines to compare against: 10-cell `--hic` (broken) phased SNPs bsgenova 42,081/69,952/54,076,
naive 18,289/27,884/21,456 — but those had <1 kb blocks; the Step-A fixed numbers are the real
baseline for block length.

---

## Session 5 (2026-08-24) — TASK 1 EXECUTION LAUNCHED (read this to resume)

**User request (verbatim intent):** "Take a look at the pairing issue with snM3C-seq data,
resolve it, then run the whole phasing pipeline again on 10 snM3C-seq for each donor and 100
snM3C-seq for only donor-1." Runs as one of two parallel disconnect-safe tasks this session (the
sibling is snMC-seq downsampling, in `snMC-seq_SNP_phasing_experiment.md`).

### Status of the "pairing issue" — already RESOLVED and VALIDATED (Session 4)
The pairing problem (CHECKPOINT 10): the upstream YAP/cemba snm3C pipeline splits chimeric reads
at ligation junctions and maps every segment **single-end**, so `extractHAIRS --hic 1` saw
`PE-fragments 0` and silently degraded to short-read phasing (blocks <1 kb). The fix
(`Python-scripts/repair_m3c_pairs.py`) reconstructs real R1/R2 pairs from the read-name suffixes
per contig, then `samtools fixmate` fills mate fields. **Validated on chr21 (H1930001):
max block span 281 bp → 27.3 Mb, blocks>1Mb 0 → 11.** The fix is baked into the Task-1 driver.
So Session 5 is about *running* the fix end-to-end at scale, not re-developing it.

### What is being run — `run_snM3C-seq_task1.sh` (tmux `m3c-task1`, setsid nohup nice -n 10)
- **Step A — re-phase the 3 existing 10-cell merged BAMs WITH the pairing fix.** Inputs on the
  mount: `Science-snM3C-seq/merged-BAM/<donor>_CX45_snM3Cseq_10cells_merged.bam` (all 3 donors,
  present + valid). For each: repair pairs → extractHAIRS `--hic 1` → HAPCUT2 `--hic 1`, reusing
  the existing preprocessed het VCFs. Outputs → `phased-from-HapCUT2-hicfix/from-{bsgenova,naive}/`
  (the old broken `phased-from-HapCUT2/` is kept for the before/after record). Deliverable: 6
  corrected block sets with real megabase spans vs the broken <1 kb baseline.
- **Step B — download 90 more H1930001 CX45 snM3C BAMs (ranks 11–100 by tar size)** via
  `download_snM3C-seq_more.sh` → `Science-snM3C-seq-100cell/H1930001/` (dir
  `data/snM3C-seq_100cell_experiment/H1930001/`). Same donor ⇒ same ancestry; skips the 10 already
  present.
- **Step C — 100-cell pseudo-bulk pipeline for H1930001** (C1 sort+index+merge 100 cells →
  C2 repair pairs → C3 sharded bsgenova+naive genotyping on the merged BAM → C4 HapCUT2
  preprocessing → C5 HAPCUT2 `--hic 1`). Outputs → `data/snM3C-seq_100cell_experiment/`.
  **Deliverable: does 10→100 cells improve phased-SNP count / block N50 / max span** vs the
  10-cell fixed baseline from Step A.

### Execution environment (node was ephemeral again — rebuilt this session)
- Node up 20 d but **all conda envs gone, /tmp wiped, tmux not installed** at session start.
  Rebuilt via `/tmp/bootstrap.sh` (setsid nohup): envs `bsgenova`, `htslib-tools`, `hapcut2`,
  `samtools`(1.9), `tmux`; reference re-staged to `/tmp/snm3c_phase/hg38_chrL.fa`; code staged to
  `/tmp/snm3c_local/` (drivers, callers, `HapCUT2_run.sh` = the validated reconstruction,
  `repair_m3c_pairs.py`, bsgenova tool tree). SHAPEIT5 deliberately NOT rebuilt (established
  unusable on this callset).
- Resumable via `/tmp/m3c_task1/.stage_{A,B,C1..C5}_done`; compute in /tmp, publish to mount via
  retrying non-fatal `cpv()`; monitor `/tmp/m3c_task1/progress.log`. Resource plan: samtools
  NP=16, merge PWORK=24, sharded genotyping SHARD_CONC=26 × -P2 (~52 procs). Runs concurrently
  with Task 2 (~104 sharded procs combined) on a 256-core / 503 GB node — comfortably within budget.

### Housekeeping done this session (per user request)
- Deleted `logs-sync-20260820.bundle` (untracked git-sync transport artifact).
- Compacted `shell-scripts/` from 20 → 7 top-level: HapCUT2_run triplet collapsed to one canonical
  `HapCUT2_run.sh` (the validated reconstruction); `verify_published_artifacts.sh` +
  `validate_1000_bams.sh` + `tmux_attach.sh` merged into `pipeline_utils.sh` (subcommands
  `verify` / `validate-bams` / `tmux`); nine one-shot/superseded scripts moved to
  `shell-scripts/archive/`. Active drivers untouched.

## STEP A COMPLETE (2026-08-24, 18:35–22:07 UTC) — the read-pairing fix delivers Mb-scale blocks

Re-phased all three 10-cell snM3C merged BAMs with `repair_m3c_pairs.py` applied before
`extractHAIRS --hic 1`. All 6/6 phased sets published to
`phased-from-HapCUT2-hicfix/from-{bsgenova,naive}/`. Pair repair on H1930001:
`pairs=31562445 singles=35145989 dropped_unmapped=0`.

**This overturns the previous session's conclusion.** The earlier writeup blamed h-trans
estimation at low depth for the short blocks and recommended `--hic_htrans_file` as the next
experiment. That was wrong: the bottleneck was the read pairing. With the fix, and h-trans still
estimated from data (`--hic 1`, unchanged), blocks go from sub-kb to tens of megabases.

| donor / caller | blocks | phased SNPs | SNPs/block | largest block (SNPs) | span N50 | largest span | blocks >=1 Mb |
|---|---|---|---|---|---|---|---|
| H1930001 bsgenova | 19,144 | 48,156 | 2.52 | 1,059 | 22.33 Mb | **83.84 Mb** | 274 |
| H1930002 bsgenova | 33,331 | 80,869 | 2.43 | 1,414 | 19.71 Mb | **71.11 Mb** | 554 |
| H1930004 bsgenova | 24,662 | 61,647 | 2.50 | 1,260 | 24.39 Mb | **69.75 Mb** | 488 |
| H1930001 naive | 5,565 | 20,166 | 3.62 | 1,201 | 23.21 Mb | **37.05 Mb** | 82 |
| H1930002 naive | 9,616 | 31,341 | 3.26 | 1,656 | 21.86 Mb | **39.83 Mb** | 203 |
| H1930004 naive | 6,916 | 23,779 | 3.44 | 1,340 | 23.29 Mb | **56.61 Mb** | 119 |
| *pre-fix, same 6 sets* | *6,263–31,046* | *18,289–69,952* | *2.25–2.92* | *133–184* | *63–75 bp* | *897–1,024 bp* | ***0*** |

Largest block overall: **chr16:1,970,966–85,807,326 (83.84 Mb, 549 SNPs)** — essentially the whole
of chr16 (90.3 Mb). Sanity check: **0 blocks span more than one chromosome** in any of the 6 sets.
Phased-SNP counts rise only modestly (+11–16%); the change is in contiguity, not yield.

### Is the long-range signal real? Three checks

**1. Independent reproducibility (fragment split-half).** Donor 1's 357,255 bsgenova fragments were
split odd/even and each half phased separately against the same VCF:

| | blocks >=1 Mb | largest span | span N50 | phased SNPs |
|---|---|---|---|---|
| half 1 (178,628 frags) | 255 | 83.84 Mb | 22.55 Mb | 45,105 |
| half 2 (178,627 frags) | 247 | 63.38 Mb | 21.51 Mb | 45,229 |
| full (357,255 frags) | 274 | 83.84 Mb | 22.33 Mb | 48,156 |

Half the data reproduces ~90% of the Mb-scale blocks and the identical chr16 maximum, so the long
blocks are not an artefact of one particular fragment set.

**2. Long-range phase concordance between the halves.** Adjacent-site junction walks cannot score
this (consecutive shared sites are nearly always <1 kb apart — 0 comparable junctions >=1 Mb), so
`phase_concordance.py --longrange` compares all site *pairs* >=1 Mb apart that are co-phased in both
halves:

| comparison | pairs | discordance |
|---|---|---|
| half1 vs half2, >=1 Mb | 10,182 | **18.48 %** |
| half1 vs half2, >=10 Mb | 8,292 | **18.52 %** |
| null control (alleles randomised) | 10,182 | 50.18 % |

Far from chance, and **flat from 1 Mb to 10 Mb+** — the agreement does not decay with distance,
which is the Hi-C signature. But ~1 in 5 long-range relationships flips between independent halves.
Each half carries half the coverage, so 18.5% is a conservative upper bound on the full-data rate.

**3. Cross-caller check — underpowered, reported for completeness.** bsgenova vs naive on donor 1:
8,073 shared phased sites, 4,927 comparable junctions, 12.89% switch rate overall (12.12% <1 kb).
Only 3 junctions >=1 Mb, so this says nothing about long-range quality.

### Honest reading

The fix delivers genuine chromosome-arm-scale linkage — the specific promise of 3C that the
previous session concluded was "not yet realised". Two caveats belong with it:

- **The long blocks are sparse skeletons.** The 274 blocks >=1 Mb in H1930001/bsgenova hold only
  3,628 SNPs — 7.5% of all phased SNPs, a median of **0.72 SNPs per Mb**. They tile 2.70 Gb (about
  one pass over the genome), but the *bulk* of phased SNPs still sit in 2-SNP local blocks nested
  inside them. Median block is still 2 SNPs / 24 bp.
- **Switch noise is substantial** (>=18.5% at Mb scale, above) and there is still no truth haplotype
  for these donors, so no absolute switch-error rate can be quoted — only this internal bound.

### Tooling added this session (published to `workspace/Python-scripts/`, md5-verified)

- `block_stats.py` — HapCUT2 `.blocks` summariser: block/SNP counts, SNPs-per-block mean/median/N50,
  span median/N50/max, counts of blocks >=100 kb and >=1 Mb. `--tsv` for table output. Spans are
  computed from first/last variant position, not the header `SPAN` field.
- `phase_concordance.py` — cross-checks two `.blocks` files phasing the same BAM: adjacent-junction
  switch rate binned by distance, plus `--longrange` all-pairs orientation discordance at >=1 Mb.

Note: HAPCUT2's own log line `N50 haplotype length is 0.00 kilobases` is not a contradiction — it
is dominated by the thousands of trivial 2-SNP blocks; the span-weighted N50 above is the metric
that reflects the long blocks.

### CHECKPOINTS (Task 1)
- [x] envs rebuilt (bsgenova, htslib-tools, hapcut2, samtools, tmux)
- [x] Step A: 6 fixed `--hic` phased sets in phased-from-HapCUT2-hicfix/ — DONE 22:07 UTC, spans 37–84 Mb, 82–554 blocks >=1 Mb per set (pre-fix: 0)
- [x] Step B: 100 H1930001 CX45 snM3C BAMs present (10 old + 90 new) — DONE 01:25 UTC
- [x] Step C1: 100-cell coord-sorted merged BAM + .bai — DONE 02:03 UTC, 843,150,557 reads (45.5 GB); .bam now publishing to mount (driver only sent the .bai)
- [x] Step C2: paired (repaired) 100-cell BAM — DONE 07:19 UTC, pairs=273,535,308 singles=296,079,941 dropped=0
- [x] Step C3: bsgenova + naive callsets (100-cell) — DONE 10:39 UTC, 4,939,635 / 5,242,025 records
- [x] Step C4: 2 preprocessed het VCFs — DONE 10:47 UTC, het 3,897,860 (bsg) / 4,389,230 (naive)
- [~] Step C5: 2 HAPCUT2 `--hic 1` phased sets — RUNNING since 11:17/11:19 UTC (715,422 / 825,583 components)
- [ ] final outputs copied to workspace/data (per user: always keep a mount copy)

## STEP B + C1–C4 COMPLETE (2026-08-25, 01:25–10:47 UTC) — 100-cell pseudo-bulk built and genotyped

Step A's fix carried straight into the scale-up. Everything up to the phasing call is done; **C5
(HAPCUT2 `--hic 1`) is still running at the time of writing** — results section follows below when
it lands.

### Step B — 90 more H1930001 CX45 snM3C cells (01:25 UTC)
`downloaded(new)=90  existing(old)=10  total=100`. The 90 new per-cell BAMs are on the mount at
`data/snM3C-seq_100cell_experiment/H1930001/`; the original 10 stay in the Step-A tree. One
transient `Remote I/O error` reading `...P1-3-M14-A18.3C.sorted.bam` off the mount during staging —
retried clean, and the final merge read count matched the expected sum exactly, so no cell was lost.

### C1 — merge (02:03 UTC)
100/100 cells coord-sorted + indexed, then merged: **843,150,557 reads** (exactly the expected sum),
45,529,743,420 bytes.

### C2 — pair repair (07:19 UTC, 5 h 16 m)
`repair_m3c_pairs: pairs=273535308 singles=296079941 dropped_unmapped=0` → 45,205,381,144 bytes.
Cross-check: 273,535,308×2 + 296,079,941 = **843,150,557** = the merged read count, so the repair is
read-conservative, as at 10 cells.

### C3/C4 — genotyping + HapCUT2 preprocessing (10:47 UTC)
Sharded on the *original* (unrepaired) merged BAM, as designed: `bsgenova records=4,939,635`,
`naive records=5,242,025`; after preprocessing, `het bsgenova=3,897,860`, `het naive=4,389,230`.
Both callsets + preprocessed VCFs (+`.gz`/`.tbi`/logs) published to the mount and present.

### What 10 → 100 cells actually bought at the *input* stage

| | 10-cell | 100-cell | ratio |
|---|---|---|---|
| reads in merged BAM | 98,270,879 | 843,150,557 | 8.58x |
| het sites, bsgenova | 163,894 | 3,897,860 | **23.8x** |
| het sites, naive | 60,099 | 4,389,230 | **73.0x** |
| extractHAIRS fragments, bsgenova | 357,255 | 12,461,865 | 34.9x |
| extractHAIRS fragments, naive | 341,752 | 19,023,334 | 55.7x |
| variants/read, bsgenova | 5.25 | 5.39 | — |
| variants/read, naive | 7.49 | 4.82 | — |
| graph components entering HAPCUT2, bsgenova | 19,406 | 715,422 | 36.9x |
| graph components entering HAPCUT2, naive | 5,642 | 825,583 | 146x |

Het discovery scaled *faster than linearly in reads* (23.8x / 73x for 8.6x the data) — at 10 cells
most of the genome simply had no coverage to call from, so extra depth converts into new callable
sites rather than deeper pileups at known ones. The naive caller gains most because its 10-cell
callset was the more depth-starved of the two. Note this also means the 100-cell run is **not** the
same phasing problem scaled up: it is a ~24–73x larger variant graph, which is why C5 is taking
hours rather than the ~8 minutes it took per 10-cell set.

Caution for the results section: more het sites is not itself an improvement. Whether the extra
sites are real, and whether they link into longer blocks or just add more 2-SNP islands, is exactly
what the C5 block stats have to answer.

### Deliverable-retention gap found and closed (13:56 UTC)

The driver publishes only the merged BAM's `.bai` to the mount (`run_snM3C-seq_task1.sh` line 166) —
the 45.5 GB `.bam` itself was never attempted, so after the next node wipe the 100-cell pseudo-bulk
would have been unrecoverable without redoing C1 (and C2's 5 h repair). Same class of loss as the
1000-cell BAM in Task 2. Launched `publish_big_file.sh` on it (11 x 4 GiB chunks, fenced
`taskset -c 96-159 nice -n 19 ionice -c 3`), verifying size + quickcheck +
`idxstats reads=843150557` at the end.

**That first attempt died and taught us the poisoning rule is sharper than recorded.** Chunk 0 wrote
fine (4 GiB in 8 min, ~8.5 MB/s); chunk 1's write hit `Remote I/O error`, and attempts 2 and 3 then
`failed to open` the file at all, while `stat` reported *No such file or directory*. So
`merged-BAM/H1930001_CX45_snM3Cseq_100cells_merged.bam` is now **permanently poisoned** — do not
publish there again. Two refinements to what the Task 2 entry recorded:

1. **Poisoning is not specific to `cp`, and not to whole-file writes.** A chunked `dd` pwrite into an
   already-partly-written object poisons the path just as thoroughly.
2. **`publish_big_file.sh`'s 3 same-path retries are worthless.** Once the first write fails the path
   is dead, so attempts 2-3 can only fail to open it. The retry loop just burns ~10 min before
   reporting a failure it could have reported immediately.

The mount itself was fine throughout — a 4 MiB probe write to a *fresh* name in the same directory
succeeded immediately. (A concurrent symptom of the same stall: `posix_spawn '/bin/bash'` returned
`EREMOTEIO` for a command whose cwd was on the mount.)

So `shell-scripts/publish_with_failover.sh` now wraps it: probe the mount with a small fresh-name
write (wait up to 1 h if the mount is genuinely degraded, rather than burning a path), stage the
`.bai` beside the target so the read-count verify can run, then run the chunked publish — and on
*any* failure abandon that path entirely and restart at the next one with fresh chunk state, up to 5
attempts. Destination paths are `merged-BAM/published/`, then `published-2/`, ... which keeps the
canonical filename (the same trick as Task 2's `rebuilt/`). Relaunched 16:12 UTC at
`merged-BAM/published/`; log `/tmp/m3c_task1/publish_merged100_failover.log`.

## 2026-08-25 (16:00-17:00 UTC) — storage audit: Step A validation evidence published

Audited every local deliverable in `/tmp/m3c_task1` against a size index of the published trees
(363 files across the three experiment folders). Results:

| artefact | status |
|---|---|
| Step A phased sets, 6 sets x 5 files (`.blocks`, `.blocks.phased.VCF`, `.fragments`, 2 logs) | 30/30 on the mount |
| Step C callsets (bsgenova + naive, 100-cell) | 4/4 on the mount |
| Step C preprocessed het VCFs (+`.gz`/`.tbi`/logs) | 8/8 on the mount |
| **Step A validation evidence** (split-half reproducibility + block-stats table) | **was missing — 10 files published this session** |
| Step C5 phased sets | not yet — C5 still running; the driver publishes them on completion |
| 100-cell merged BAM (45.5 GB) | publish in flight, chunk 2/11 at 16:35 UTC |
| 100-cell repaired BAM (45.2 GB) | **deliberately not published** — see below |

### Step A validation evidence now on the mount

The "is the long-range signal real?" section of the Step A writeup rested on files that only existed
in /tmp and would have died with the next node wipe — the numbers would have survived in this log
with no way to re-derive them short of re-running the whole check. Published to
`phased-from-HapCUT2-hicfix/validation/`:

- `half1.*` / `half2.*` (`.blocks`, `.blocks.phased.VCF`, `.fragments`, `.HAPCUT2.log`) — the odd/even
  fragment split-half phasings of donor 1 that reproduced ~90% of the Mb-scale blocks and the
  identical chr16 maximum, and that the 18.5% long-range discordance figure was computed from.
- `stepA_block_stats.tsv` — the block-stats table behind the Step A results table.
- `cell_reads.txt` — per-cell read counts for the 100-cell set.

All 10 verified byte-exact against their local copies.

Not published from `compare/`: the `prefix.*.blocks` files, which are copies of the *pre-fix* broken
phased sets already living on the mount under `phased-from-HapCUT2/`.

### USER DECISION: keep the repaired 100-cell BAM on the mount

`repaired/H1930001_CX45_snM3Cseq_100cells_merged.paired.cs.bam` (45.2 GB) is the direct input to
`extractHAIRS`. It is reproducible — `repair_m3c_pairs.py` applied to the merged BAM, one
deterministic 5 h 16 m pass — so it was raised as a keep-or-drop question alongside Task 2's
downsampled BAMs. **The user chose to keep this one and to omit the Task 2 downsampled BAMs.**

Verified locally before publishing: `quickcheck OK`, `idxstats` = **843,150,557 reads**, identical to
the merged BAM it was derived from, which independently re-confirms the repair is read-conservative.
That count is what the publish verifies against at the far end.

Publishing to `snM3C-seq_100cell_experiment/repaired-BAM/published/` via
`publish_with_failover.sh`, **chained to start only after the merged-BAM publish exits**
(`/tmp/m3c_task1/chain_repaired_publish.sh`, log `publish_repaired100.log`). Strictly sequential on
purpose: two concurrent multi-GB writes is what contended during the Task 2 publish, and here a
write failure does not merely slow things down, it poisons the destination path permanently.

### USER DECISION: keep the three 10-cell repaired BAMs too

Also published, on the same request: the Step A repaired BAMs for all three donors, into
`snM3C-seq_SNP_phasing_experiment/repaired-BAM/published/` (mirroring the 100-cell layout). Each
verified locally first:

| donor | size | quickcheck | reads |
|---|---|---|---|
| H1930001 | 6,350,117,015 | OK | 98,270,879 |
| H1930002 | 6,313,601,973 | OK | 103,086,888 |
| H1930004 | 5,931,272,503 | OK | 95,040,464 |

H1930001's 98,270,879 is an exact independent check on the Step A pair-repair figures:
31,562,445 x 2 + 35,145,989 = 98,270,879. Those counts are what each publish verifies against.

Run by `chain_repaired10_publish.sh`, which waits on the 100-cell repaired chain and then does the
three strictly one at a time. Full publish order this session, all sequential:
**merged 100-cell -> repaired 100-cell -> 10-cell x3.** These are ~6 GB each, i.e. small enough that
a plain `cp` would probably survive — but the 100-cell failure poisoned its path only ~4-8 GB into
the object, so size is not a safety argument on this mount. Same chunked failover publisher for all
of them.

Still not published, with the user's agreement: Task 2's 189 GB of downsampled BAMs.

### Verification protocol correction (applies to both tasks)

While publishing 12 small logs on the Task 2 side, 10 of 12 "failed" a size check run immediately
after `cp` returned; re-checked a minute later, all 12 were byte-correct. The settling behaviour
already documented for large chunked writes applies to small files too — **an immediate post-write
`stat` on this mount is not a valid check**, and a verify pass that runs one will report phantom
failures. Retry the verification after a delay instead.

## 2026-08-26 (15:56-16:10 UTC) — publish chains found dead-locked on a zombie; Step C5 bsgenova results

Resumed to find the node **not** wiped for once (tmux sessions, conda envs and `/tmp` all intact
from the 08-25 session), one HAPCUT2 still running, and both publish chains silently stuck.

### The 21-hour stall: `kill -0` succeeds on a zombie

`chain_repaired_publish.sh` and `chain_repaired10_publish.sh` serialised themselves against the
upstream publish with the obvious idiom:

```bash
while kill -0 "$WAITPID" 2>/dev/null; do sleep 60; done
```

The merged-BAM publish **succeeded at 18:42 UTC on 08-25** (all 11 chunks, size + `quickcheck` +
`idxstats reads=843150557` all verified). Both chains were nevertheless still sleeping 21 h later,
and neither repaired BAM had been published.

**Cause: PID 1 in this container does not reap orphans.** The finished publish (pid 958508) was
reparented to init and stayed `<defunct>` forever, and `kill -0` on a zombie *returns success* —
the PID is still in the process table, it just has no process behind it. The wait predicate can
therefore never become false. (Corroborating evidence was in plain sight: `ps` shows dozens of
`[claude.exe] <defunct>` entries up to two days old. Nothing on this node is ever reaped.)

This is a general property of the node, not of these scripts, so **never serialise background work
here by polling a PID.** Either check for the zombie state explicitly (`awk '{print $3}'
/proc/$pid/stat` = `Z`), or — better, and what was done — do not wait across processes at all.

`publish_all_repaired.sh` replaces both chains with a **single process** that runs all four
publishes sequentially in one script. Same strict one-write-at-a-time guarantee (two concurrent
multi-GB writes is what contended during the Task 2 publish, and a failed write poisons the
destination path permanently), with no cross-process PID waiting to get wrong. Relaunched 16:05 UTC
fenced `taskset -c 96-159 nice -n 19 ionice -c 3`; log `/tmp/m3c_task1/publish_all_repaired.log`.
Order: repaired 100-cell (45.2 GB) → 10-cell H1930001 / H1930002 / H1930004 (~6 GB each), ~4 h total.

The two stuck chains were killed before relaunching, deliberately: they were inert, but had that
zombie ever been reaped they would have woken up and written to the *same* `published/` path as the
new driver, and a concurrent write is exactly what poisons a path permanently.

### Step C5 (bsgenova, 100 cells) — the contiguity question, answered

bsgenova finished 08-25 20:33 (9 h 15 m of Max-Likelihood-Cut). HAPCUT2's own summary:
**N50 haplotype length 19,854 kb**, 715,422 non-trivial components, max-degree 182,458,
2,275,683 connected variants, coverage-per-variant 32.7. Full block stats vs the Step A 10-cell
phasing of the same donor and caller:

| | 10-cell (Step A) | 100-cell (Step C5) | ratio |
|---|---|---|---|
| blocks | 19,144 | 707,196 | 36.9x |
| phased SNPs | 48,156 | 2,236,344 | **46.4x** |
| het sites in | 163,894 | 3,897,860 | 23.8x |
| **phasing rate** | **29.4%** | **57.4%** | **1.95x** |
| mean SNPs/block | 2.52 | 3.16 | 1.25x |
| SNP-N50 / max SNPs in a block | 2 / 1,059 | 3 / 7,043 | — / 6.6x |
| span N50 | 22,331,364 bp | 22,429,295 bp | **1.004x** |
| max block span | 83,836,361 bp | 127,161,910 bp | 1.52x |
| blocks >= 1 Mb | 274 | 84,750 | **309x** |
| blocks 100 kb - 1 Mb | 210 | 53,699 | 256x |

(Note `blocks_ge_100kb` in `block_stats.py` is the half-open bucket [100 kb, 1 Mb), not a cumulative
count — that is why it reads lower than the >= 1 Mb column in both rows.)

**Answer to the question C5 was set up to settle: the extra het sites link, they do not just pile up
as 2-SNP islands.** The 08-25 caution was that 24x more het sites is not itself an improvement. It
turns out to be one — but the improvement is in *how many* long blocks there are, not in how long
they are:

- **Contiguity per block did not improve at all.** Span N50 moved 22.33 → 22.43 Mb, a 0.4% change.
  Ten cells was already enough to produce Mb-scale blocks; that was the read-pairing fix's win, and
  it was already saturated. The 3C contact range sets the block length, and more cells do not
  lengthen a contact.
- **The number of long blocks exploded**: 274 → 84,750 blocks of >= 1 Mb, 309x for 8.6x the reads.
  What extra depth buys is *more regions that have enough het sites to be linkable at all*, not
  longer links in regions that already were.
- **Phasing rate nearly doubled** (29.4% → 57.4%), and this is the cleanest contrast with Task 2.
  In snMC-seq the rate *fell* with depth (46.6% at 20 cells → 42.2% at 1000) and was called out as
  a low-depth selection artefact. Here it rises steeply with depth, because a 3C fragment can reach
  a distant het site the moment that site becomes callable — so a newly discovered site usually
  arrives already linkable, rather than as an isolated island. That difference is attributable to
  the Hi-C contacts and not to depth per se, since Task 2 had far more depth and did worse.
- The median block is still tiny (SNP-N50 = 3). The distribution is strongly bimodal: a large mass
  of 2-3 SNP fragments plus a heavy tail of Mb-spanning blocks, and 127 Mb is now the largest single
  block — roughly half of chr1, i.e. chromosome-arm scale.

### Still open

- **C5 naive is still running** — 28 h 39 m of HAPCUT2 wall-clock as of 16:00 UTC, vs 9 h 15 m for
  bsgenova. It carries a larger graph (4,389,230 het sites / 825,583 components / 19,023,334
  fragments) and a *lower* variants-per-read (4.82 vs 5.39), i.e. more nodes joined by weaker
  evidence, which is the expensive direction for Max-Likelihood-Cut. No `.blocks` file is written
  until it finishes. Not yet a cause for concern, but it is the long pole; the driver
  (`run_snM3C-seq_task1.sh`, pid 13717) is alive and will publish the C5 phased sets on completion.
- Publishing of the four repaired BAMs is in flight (see above).

## 2026-08-27 (18:41-  UTC) — Task 1 complete; C5 naive results; cross-caller long-range check; teardown

Resumed to a node that was again **not** wiped (tmux sessions, `/opt/conda` envs and `/tmp` all
intact). Nothing was running: every driver launched on 08-25/08-26 had reached a terminal marker
while the session was down. This entry records what those unattended runs produced, one new
analysis, the publication of everything still local, and the teardown of `/tmp`.

### What finished unattended

| run | marker | finished (UTC) |
|---|---|---|
| Step C5 naive HAPCUT2 | `C5_OK naive` | 08-26 17:19 |
| Task 1 driver overall | `TASK1_ALL_DONE` | 08-26 17:24 |
| 4x repaired-BAM publish | `REPAIRED_ALL_OK` (failures=0) | 08-26 23:42 |
| C5 phased-set publish | `C5_BSG_SOME_FAILED` (failures=1) | 08-26 23:49 |

The single reported publish failure was **not real**. `publish_c5_bsgenova.sh` waits 300 s after
issuing its copies and then size-checks each one; it reported
`VERIFY FAIL ... blocks.phased.VCF (src=405484240 dst=missing)`. Re-checked 19 h later the file is
present at exactly 405,484,240 B, and all ten C5 files (both callers) match their local sizes byte
for byte. This is the third instance of the same thing, and it sharpens the rule already in this
log: **300 s is not a long enough settle window on this mount.** A verify pass has to be a separate
step run much later — minutes are not reliably enough — or the driver will keep manufacturing
phantom failures and, worse, tempt a future session into re-publishing over a path that is
actually fine (which is how paths get poisoned).

The single-process `publish_all_repaired.sh` rewrite worked exactly as intended: four sequential
multi-GB publishes, no cross-process PID waiting, no repeat of the 21-hour zombie deadlock. The
100-cell repaired BAM again lost its first path and landed at `repaired-BAM/published-2/`; the
three 10-cell BAMs went to `published/` on attempt 1.

> **Hazard left on the mount.** `repaired-BAM/published/H1930001_CX45_snM3Cseq_100cells_merged.paired.cs.bam`
> is the **truncated 30,261,313,536 B remnant** of the poisoned first attempt, sitting next to a
> correctly-named `.bai`. The good 45,205,381,144 B copy is in `published-2/`. With no icommands on
> this node the remnant cannot be renamed or removed, so a `README-DO-NOT-USE.txt` has been written
> beside it. Anyone globbing `repaired-BAM/*/…paired.cs.bam` will otherwise hit the bad one first.

### Step C5 naive — the second caller at 100 cells

HAPCUT2 took **~30.0 h** of Max-Likelihood-Cut (11:20:38 08-25 → 17:18:46 08-26) against **9.2 h**
for bsgenova, as anticipated from its larger, weaker graph. Its own summary: N50 haplotype length
**36,246 kb**, 825,583 non-trivial components, max-degree 182,451, 2,981,517 connected variants,
coverage-per-variant 36.79.

Against the Step A 10-cell naive phasing of the same donor and caller:

| | 10-cell (Step A) | 100-cell (Step C5) | ratio |
|---|---|---|---|
| blocks | 5,565 | 818,118 | 147x |
| phased SNPs | 20,166 | 2,937,592 | **146x** |
| het sites in | 60,099 | 4,389,230 | 73.0x |
| **phasing rate** | **33.6%** | **66.9%** | **1.99x** |
| mean SNPs/block | 3.62 | 3.59 | 0.99x |
| SNP-N50 / max SNPs in a block | 3 / 1,201 | 4 / 14,317 | — / 11.9x |
| span N50 | 23,209,324 bp | 23,622,262 bp | **1.018x** |
| max block span | 37,047,463 bp | 155,444,709 bp | 4.20x |
| blocks >= 1 Mb | 82 | 121,977 | **1,487x** |
| blocks 100 kb - 1 Mb | 61 | 70,813 | 1,161x |

Same shape as bsgenova, and it reproduces the C5 conclusion independently: **span N50 is flat
(+1.8%) while the count of Mb-scale blocks explodes.** Two callers with different variant sets,
different fragment sets and a 3x difference in runtime agree that extra depth buys *more linkable
regions*, not *longer links*. That was the single most important claim from 08-26 and it is no
longer resting on one callset.

### bsgenova vs naive at 100 cells — the ordering inverts

| | bsgenova | naive | winner |
|---|---|---|---|
| het sites in | 3,897,860 | 4,389,230 | naive +12.6% |
| fragments | 12,461,865 | 19,023,334 | naive +52.7% |
| mean variants per read | 5.39 | 4.82 | bsgenova |
| coverage per variant | 32.7 | 36.8 | naive |
| blocks | 707,196 | 818,118 | — |
| **phased SNPs** | 2,236,344 | **2,937,592** | **naive +31.4%** |
| **phasing rate** | 57.4% | **66.9%** | **naive** |
| span N50 | 22.43 Mb | **23.62 Mb** | naive +5.3% |
| max block span | 127.16 Mb | **155.44 Mb** | naive +22.2% |
| blocks >= 1 Mb | 84,750 | **121,977** | naive +43.9% |
| HAPCUT2 wall-clock | **9.2 h** | 30.0 h | bsgenova 3.3x faster |

**At 10 cells bsgenova phased 2.4x more SNPs than naive (48,156 vs 20,166). At 100 cells naive
phases 31% more.** The crossover is the headline result of this session.

The mechanism is visible in the het counts. Going 10 -> 100 cells, bsgenova's het callset grew
23.8x but naive's grew 73.0x. bsgenova's posterior model needs depth before it will call a site
heterozygous, so at 10 cells it converts sparse coverage into calls that naive's fixed thresholds
throw away — a real advantage, and the reason it was made the primary caller. But that advantage is
a *low-depth* advantage. By 100 cells (8.6x the reads, ~33-37x per variant) the depth that bsgenova
needed is simply present, naive's thresholds clear on their own, and naive ends up with the larger
callset. Crucially the extra sites are not junk: naive's phasing *rate* is also higher (66.9% vs
57.4%), so a larger share of what it calls goes on to be linked, and its long-block counts are
higher too.

This does not retire bsgenova. It is 3.3x cheaper here, it is the better caller in exactly the
low-coverage regime this project cares about, and — see below — the two callers disagree enough at
long range that neither should be treated as ground truth. What it does mean is that **the choice
of caller is depth-dependent, and the 10-cell result must not be extrapolated to pseudo-bulk.**

### Are the Mb-scale blocks still sparse skeletons?

Step A found the 10-cell long blocks were skeletons: 7.5% of phased SNPs at a median 0.72 SNPs/Mb.
Rerun at 100 cells (`stepC_bigblock_content.txt`):

| | bsgenova | naive |
|---|---|---|
| SNPs inside blocks >= 1 Mb | 486,440 (**21.8%** of phased) | 868,941 (**29.6%**) |
| median SNPs per Mb | 0.73 | 0.79 |
| total bp spanned by those blocks | 914.6 Gb | 1,409.2 Gb |

The long blocks now hold 3-4x the share of the phasing they did at 10 cells, so they are less
marginal. But **median density is unchanged** (0.72 -> 0.73/0.79 SNPs per Mb): they are the same
sparse skeletons, there are just vastly more of them. The third row is the caution — 915 Gb and
1,409 Gb of span against a 3.1 Gb genome means these blocks **overlap roughly 300-450 deep**. The
Step A figure of "2.70 Gb tiled, about one pass over the genome" does not survive to this depth, and
the count of blocks >= 1 Mb must not be read as genome coverage. It is a heavily nested pile.

### Cross-caller long-range concordance — the underpowered check, rerun

Step A's third validity check (bsgenova vs naive) had 4,927 comparable junctions but only **3** at
>= 1 Mb, so it was reported for completeness and nothing was concluded. At 100 cells there is real
power. `phase_concordance.py bsgenova naive --longrange`:

```
bsgenova: 2,229,600 phased sites | naive: 2,930,846 | shared: 995,392
comparable junctions: 466,992   switch rate 6.89%
      <1 kb   459,946    6.73%
   1-100 kb     7,028   17.19%
 100-1000 kb        15   26.67%
long-range pairs >=1 Mb : 19,383   discordance 53.20%
long-range pairs >=10 Mb:  7,487   discordance 42.29%
```

**Short-range phasing is in good shape: 6.73% switch rate below 1 kb**, less than half the 12.89%
the two callers managed at 10 cells, and it degrades smoothly with distance (17.2% at 1-100 kb).

The long-range numbers as printed say 53.2%, i.e. chance. **They are a sampling artefact and should
not be quoted.** `long_range()` strides its outer loop to cap pair count, which at this scale leaves
only ~10 anchor sites per chromosome; all ~19,383 pairs come from a few hundred anchors, so one
anchor sitting in a block the two callers oriented oppositely contributes hundreds of "switches".
The effective sample size is the anchor count, not the pair count.

Re-measured with `longrange_structure.py`, which pairs **each** shared site with the nearest
co-phased site >= 1 Mb downstream (one observation per anchor, ~53k independent-ish pairs):

| estimator | pairs | discordance |
|---|---|---|
| strided anchors (as printed above) | 19,383 | 53.20% |
| **one pair per anchor** | **53,539** | **28.01%** |
| same-caller split-half, 08-25 | 10,182 | 18.48% |
| null control (alleles randomised), 08-25 | 10,182 | 50.18% |

**28.0% is the number to use.** It is far from the 50.2% null, so the Mb-scale cross-caller
agreement is real signal — but it is meaningfully worse than the 18.5% the same caller achieved
against itself on split fragments, which is what should be expected: that comparison shared a
variant set and this one does not (only 995,392 of 2.2-2.9 M sites are even shared).

A companion check asked whether the residual disagreement is whole blocks flipped end-to-end or
site-level noise, by histogramming each anchor's own discordance fraction. Over 52 anchors with
>= 5 pairs the distribution is **unimodal and broad (mode 5-35%, nothing at all above 95%)**, not
piled at 0% and 100%. So this is not a few misoriented blocks; the long-range links are individually
noisy. Caveat: 52 anchors is thin, and this sub-check is suggestive rather than conclusive.

**Note for future sessions: `phase_concordance.py --longrange` is not safe to quote at pseudo-bulk
scale.** Its striding was tuned for the ~48k-site 10-cell callsets. Either raise
`max_pairs_per_chrom` until the anchor count is a meaningful fraction of the sites, or use the
one-pair-per-anchor estimator.

### Honest reading, updated

- The 08-26 conclusion (more depth = more long blocks, not longer ones) is **confirmed by a second
  independent caller**. This is now the best-supported claim in Task 1.
- The long blocks carry 3-4x more of the phasing than at 10 cells, but are the same sparse
  skeletons and overlap ~300-450 deep. "84,750 blocks >= 1 Mb" is not 84,750 Mb-scale regions.
- Cross-caller long-range agreement is **real but noisy: ~28% discordance vs a 50% null**. There is
  still no truth haplotype for these donors, so this remains an internal bound, not a switch-error
  rate.
- Caller choice is depth-dependent and inverts between 10 cells and 100.

### Published this session

Everything still local at the start of the session was pushed to the mount before teardown.

Task 1, into `snM3C-seq_100cell_experiment/phased-from-HapCUT2/`:
`stepC_block_stats_naive.tsv`, `stepC_block_stats_bsgenova.tsv` (already there),
`stepC_bigblock_content.txt`, `stepC_concordance_bsgenova_vs_naive.txt`,
`stepC_longrange_structure.txt`; driver + publish logs into `run-logs/`.

**Task 2 — a false alarm, recorded because the mistake is repeatable.** An inventory of
`/tmp/mc_downsample` against `data/snMC-seq_SNP_phasing_experiment/` found no per-fraction callsets,
no phased sets and no `saturation_summary.csv`, which looked like the entire scientific output of
Task 2 living only in `/tmp`, one node-wipe from being lost. It is not: the downsample series was
published on 08-25 to its **own top-level tree**, `data/snMC-seq_downsample_series/`, not into the
`snMC-seq_SNP_phasing_experiment/` folder the rest of Task 2 uses. A file-by-file comparison
confirms **139/139 files present and byte-identical** (24 callsets, 48 preprocessed VCFs, 60 phased
files, `saturation_summary.csv`, plus the 6 `.bai`). The naive callsets are simply renamed on
publication (`...frac0.02.naive.snv.gz` locally becomes `from-naive/...frac0.02.snv.gz`), which is
what made a naive name-based diff look like a miss.

Lesson: **check `workspace/data/` for a sibling tree before concluding something is unpublished.**
The four experiment trees are `snM3C-seq_100cell_experiment`, `snM3C-seq_SNP_phasing_experiment`,
`snMC-seq_downsample_series` and `snMC-seq_SNP_phasing_experiment`, and Task 2's outputs are split
across the last two. (Full re-audit in the snMC-seq log, 2026-08-27 entry.)

The only genuinely unpublished Task 2 artefacts were the **6 `sharded_callers.log` files** (~21 KB
each), one per fraction; published this session into each fraction's callset directory. Still not
published, deliberately: the 189 GB of downsampled BAMs (user's standing decision, `.bai` kept) and
the 16,446 per-contig shard intermediates under `gt/*/shards/` (regenerable, and 16k tiny files is
the worst possible workload for this mount).

New tooling, published to `workspace/Python-scripts/`:
- `bigblock_content.py` — what fraction of phased SNPs live in blocks >= 1 Mb, their density and
  total span. Written to test the "sparse skeleton" finding at a new depth.
- `longrange_structure.py` — unbiased one-pair-per-anchor long-range discordance, plus the
  per-anchor spread histogram that distinguishes flipped blocks from noisy links.

### Teardown

`/tmp` cleared of this experiment (487 GB: `m3c_task1` 177 GB, `mc_downsample` 310 GB, plus the
staged reference and code) after verifying every file above against the mount. Deleted content is
either published, or a regenerable intermediate (per-cell coord-sorted BAMs, per-contig shards,
staged reference, staged code — all reproducible from published inputs).

### CHECKPOINTS (Task 1) — final
- [x] envs rebuilt (bsgenova, htslib-tools, hapcut2, samtools, tmux)
- [x] Step A: 6 fixed `--hic` phased sets — DONE 08-24 22:07 UTC
- [x] Step B: 100 H1930001 CX45 snM3C BAMs — DONE 08-25 01:25 UTC
- [x] Step C1: 100-cell merged BAM, 843,150,557 reads — DONE 08-25 02:03 UTC
- [x] Step C2: repaired 100-cell BAM — DONE 08-25 07:19 UTC
- [x] Step C3: bsgenova + naive callsets — DONE 08-25 10:39 UTC
- [x] Step C4: 2 preprocessed het VCFs — DONE 08-25 10:47 UTC
- [x] Step C5: 2 HAPCUT2 `--hic 1` phased sets — DONE 08-25 20:33 (bsg) / 08-26 17:19 (naive)
- [x] final outputs copied to workspace/data, both tasks
- [x] `/tmp` cleared

## NEXT UP — planned 2026-08-27, work paused for a few days

Both experiments are complete, published and torn down (see the 2026-08-27 entry). Nothing is
running; `/tmp` is empty; the node may be wiped before the next session, which is now harmless.

Two pieces of work are queued, in this order:

**1. Scale-up to donors 2 and 4, both assays.** So far only donor **H1930001** has a large-bulk
phasing (1000 cells snMC-seq, 100 cells snM3C-seq); H1930002 and H1930004 exist only as 10-cell
sets. Plan is to repeat the pseudo-bulk pipeline at larger cell counts for **H1930002 and H1930004
in both snMC-seq and snM3C-seq**, giving three donors at comparable depth.

Points to carry in from the completed work, so they are not rediscovered:
- **Caller choice is depth-dependent.** bsgenova wins at 10 cells, naive wins at 100 (naive phased
  31% more SNPs and had the higher phasing rate). At large bulk, run **both** and compare; do not
  assume the 10-cell ordering holds.
- **Run the snM3C read-pairing fix** (`repair_m3c_pairs.py`) before phasing, and pass
  `extractHAIRS --hic 1` *and* `HAPCUT2 --hic 1`. This is what produced Mb-scale blocks at all.
- **Naive HAPCUT2 is the long pole**: ~30 h at 100 cells vs ~9 h for bsgenova, and it scales with
  graph size, not read count. Budget for it and launch it first.
- Expect the same shape of result: more depth → *more* Mb-scale blocks, not *longer* ones.
- Fence every launch (`taskset` disjoint blocks, `nice -n 19`, `ionice -c 3`) — the node hosts the
  user's VSCode-server and saturation disconnects them.

**2. Coverage and homozygosity of the results already obtained.** Not yet started. This is analysis
over the *existing* published callsets and phasings — no new phasing needed:
- **Coverage**: per-site and genome-wide depth behind the calls, and how coverage relates to which
  sites got phased. `coverage-per-variant` is already in each `HAPCUT2.log` (32.7 bsgenova / 36.8
  naive at 100 cells) as a starting point, but per-site distributions need computing from the BAMs.
- **Homozygosity**: the het/hom breakdown of the callsets, runs of homozygosity, and whether
  apparent homozygous stretches explain the gaps between phase blocks. Note the raw-vs-het counts
  already recorded per callset (e.g. 100-cell: bsgenova 4,939,635 raw → 3,897,860 het; naive
  5,242,025 → 4,389,230) — the difference is the homozygous+filtered fraction.
- Both analyses can run off the mount copies; the source BAMs are all published.

Also still open from earlier sessions, unchanged: there is **no truth haplotype** for these donors,
so no absolute switch-error rate can be quoted — only the internal bounds recorded in the log
(~18.5% same-caller split-half, ~28% cross-caller, both at >= 1 Mb, against a 50% null).

### Correction (2026-08-27, late): `rm` DOES work through the FUSE mount

Earlier entries state that because there are no icommands on this node (`ils`/`imv`/`irm`/`icp`
absent, no `~/.irods`) there is "no server-side rename or delete", and a poisoned path "can never be
reclaimed". The first half is right; the conclusion drawn from it was too strong.

Removing a duplicate session-entry file from `NCEMS-projects-conversation-logs/session-entries/`
with a plain `rm` **succeeded** — the FUSE layer implements unlink even though the icommands are
missing. So ordinary POSIX deletion is available on this mount; it is only the iRODS-native
operations (server-side copy/move, checksums, replica management) that are not.

Two things this does *not* establish, and neither has been tested:
- whether `rm` succeeds on a **poisoned** path (one whose write failed part-way) — the failure mode
  there is at the object level, not the directory-entry level, and it may behave differently;
- whether `rm` + re-write **reclaims** a poisoned path for reuse. Even if the unlink succeeds, the
  safe assumption remains "publish under a new name", because a re-poisoned path costs hours.

**Consequence for the truncated remnant.** `repaired-BAM/published/H1930001_CX45_snM3Cseq_100cells_merged.paired.cs.bam`
(30,261,313,536 B, incomplete) may now be removable after all, freeing ~30 GB. It has **not** been
removed — deleting 30 GB from the data store is irreversible and was left for the user to approve.
The `README-DO-NOT-USE.txt` beside it stands in the meantime, and the complete 45,205,381,144 B copy
in `published-2/` is verified.

### Resolved (2026-08-27): the truncated remnant was deleted — `rm` works on a poisoned path too

With the user's approval, `repaired-BAM/published/H1930001_CX45_snM3Cseq_100cells_merged.paired.cs.bam`
(the incomplete 30,261,313,536 B remnant of the 08-26 failed publish) was removed. **The `rm`
returned in 5.5 s and succeeded**, which answers the question left open in the correction above:
POSIX unlink works on this mount even for a path whose write previously failed. A poisoned path is
therefore *not* permanently unreclaimable at the directory level, and ~30 GB of quota was recovered.

Safety steps taken before deleting, worth repeating for any future remnant:
1. Confirmed the surviving copy's size (45,205,381,144 B) **and** ran `samtools quickcheck` on it —
   a size match alone would not prove the BAM is readable.
2. Re-ran both checks *after* the delete, to be sure nothing else had been disturbed. Both passed.

What is still **not** established, and should not be assumed: that a poisoned path can be safely
**re-used**. Unlinking the directory entry is not the same as clearing whatever object-level state
made the original write fail. **Keep publishing under a fresh name** (`published-2/`, `rebuilt/`,
…) — the cost of a re-poisoned path is hours, and the cost of a new name is one directory level.

`published/` now holds only the valid `.bai` and a `README.txt` recording the above; the hazard
notice it replaces is gone, since there is no longer a bad BAM to warn about.

## 2026-08-30 — Task B (donor-2/donor-4, 100 cells) still has not started: two causes, one fixed

Zero snM3C cells have been downloaded for H1930002 or H1930004. Both donors still sit at their
original **10 cells**; `grep -c 'OK mount' /tmp/m3c_d2d4/H1930002_download.log` returns 0, and
`data/snM3C-seq_100cell_experiment/H1930002/` is empty. Two independent things went wrong on 08-29.

**1. The NeMO Archive went down 15 minutes into the download** (2026-08-29 16:12 UTC) and is *still*
down at the time of writing — `nginx` answers with HTTP 502 after a 60 s upstream timeout. Full
diagnosis in the snMC-seq log's 2026-08-30 entry, including the re-audit that confirms both
downloaders build correct URLs (m3C manifest rows already carry `....3C.sorted.bam.tar`, so
`download_snM3C-seq_more.sh` is right to use col-1 verbatim rather than appending a suffix).

**2. When the gate relaunched both tasks at 20:22, Task B silently failed to launch (rc=126).**
The gate held the tmux binary path in a variable named **`TMUX`**, which tmux interprets as its own
server socket path — the first `new-session` bound a socket on top of the binary, and the second
call tried to exec it. Task A survived (it launched first) and went on to finish H1930002's snMC
200-cell run; Task B was simply gone, with nothing in any log but `rc=126`. **Never name a shell
variable `TMUX`.** The env has been rebuilt and `nemo_gate_v2.sh` uses `TMUXBIN` + `unset TMUX`.

So the 100-cell snM3C scale-up is untouched, not partially done — when it runs it starts from DL
(ranks 11-100 per donor) with nothing to clean up. `/tmp/m3c_d2d4` holds only the aborted runner log
and empty stage dirs; there are no `.<donor>.<stage>_done` markers, so no stage will be skipped.
Manifest coverage is confirmed for both donors at CX45: 5,500 3C rows for H1930002 and 2,989 for
H1930004, comfortably more than the 90 each that Task B needs.

Task B is armed to launch automatically on NeMO's recovery, in tmux session `m3c-d2d4`, fenced to
cores 32-95 with `nice -n 19 ionice -c 3` (gate: `nemo-gate`, log `/tmp/nemo_gate2.log`).

## 2026-08-31 23:20 UTC — donor-2 snM3C 100-cell COMPLETE; naive-over-bsgenova replicates

Task B finished H1930002 end to end and moved straight on to H1930004's download. Timings, for
budgeting the donor-4 half: DL 90 cells 3h13m (21:33->00:46), coord-sort+index 100 cells 16m, merge
22m (868,532,873 reads / 44,470,678,729 B, count exact), `repair_m3c_pairs.py` 4h11m (44.76 GB),
sharded genotyping 2h05m, preprocessing 7m, and `--hic 1` phasing 15h31m for both callers in
parallel. **Total ~26 h per donor**, well under the ~2 days feared from donor-1's serial run.

The merged BAM published clean on attempt 1 (size match, `quickcheck OK`, `idxstats reads=868532873
OK`).

### Results, against donor-1's 100-cell run

| donor | caller | raw | het | blocks | phased | max span | blocks >= 1 Mb |
|---|---|---|---|---|---|---|---|
| H1930001 | bsgenova | 4,939,635 | 3,897,860 | 715,422 | 2,236,344 | — | — |
| H1930001 | naive | 5,242,025 | 4,389,230 | 825,583 | 2,937,592 | — | — |
| H1930002 | bsgenova | 5,278,393 | 4,203,619 | 784,565 | 2,455,452 | 141,135,646 | 93,133 |
| H1930002 | naive | 5,771,712 | 4,885,323 | 930,226 | 3,336,467 | 163,890,201 | 139,511 |

**naive beats bsgenova at 100-cell snM3C in donor-2, replicating donor-1** — 36 % more phased SNPs
and 50 % more Mb-scale blocks. Set against the 200-cell snMC result the same week, where bsgenova
won in *both* donors by 6.5x and 2.6x, the picture is now firm and worth stating plainly:

**caller ranking tracks the assay/depth regime, not the donor.** Two donors agree within each
regime and the regimes disagree with each other. So the rule stands — run both callers on any new
experiment, and never carry a ranking across assays.

Donor-2 is also uniformly slightly richer than donor-1 (7 % more raw bsgenova calls, 10 % more
naive, 13 % more phased), consistent across every column, which is what you would expect from a
genuine donor difference rather than a run artefact.

`max_span` of 141 Mb / 164 Mb are chromosome-scale, and 93k-140k blocks exceed 1 Mb — the same
sparse-skeleton shape found for donor-1: depth buys *more* Mb-scale blocks, not longer ones.

### Verification

Shard scratch went to `/tmp/m3c_d2d4/out/H1930002/shards` (per-donor), confirming the 08-30
shared-scratch fix is active in Task B before donor-4 reaches genotyping. Donor-2's callset md5
(`095442a4`) differs from donor-1's (`64dec31f`).

H1930004 download began 23:20:32; on donor-2's timings, expect its summary around 2026-09-02 01:20.

## 2026-09-01 — Task 2 (coverage & homozygosity of existing results): both questions answered

Started the second queued item while donor-4's phasing runs. This is analysis over the
*already published* callsets and phasings — no new phasing, no BAM reads. It turns out almost
everything needed was already in the callsets: both callers emit `DP`/`DPW`/`DPC` and `GT` per
site, and the `.blocks` file carries the full `GT:GQ:GQH:DP:DPW:DPC` string for every phased
variant. So per-site depth never required going back to the 40 GB BAMs.

Three new scripts in `Python-scripts/`, stdlib-only like the rest:
- `callset_depth.py` — genotype-class breakdown and DP histogram of a callset VCF, plus
  per-window tallies. Accumulates depth as a histogram, so memory is flat.
- `phasing_vs_depth.py` — joins the preprocessed het VCF (what HapCUT2 was offered) against the
  `.blocks` file (what came out phased), on `(chrom, pos)`. Emits phasing rate by DP bin, an
  inter-block gap table, and per-window het/phased counts.
- `roh_analysis.py` — classifies fixed-width windows as `normal` / `low_het` / `roh` / `blind`
  and reports the phasing rate in each.

Run over all three donors × both callers at snM3C 100-cell (donor-4's callsets are published
already even though its phasing is still running). ~4 min total, fenced to cores 0-31 with
`nice -n 19 ionice -c 3` so it could not touch the running phasing job.

### Finding 1 — phasing success is a clean monotonic function of the site's own depth

Phasing rate by DP bin, all four completed donor×caller runs:

| DP bin | d1 bsgenova | d1 naive | d2 bsgenova | d2 naive |
|---|---|---|---|---|
| 10–14 | 0.460 | 0.466 | 0.465 | 0.478 |
| 15–19 | 0.522 | 0.561 | 0.527 | 0.570 |
| 20–24 | 0.569 | 0.626 | 0.572 | 0.637 |
| 25–29 | 0.611 | 0.674 | 0.612 | 0.683 |
| 30–39 | 0.665 | 0.718 | 0.665 | 0.727 |
| 40–49 | 0.730 | 0.768 | 0.727 | 0.775 |
| 50–74 | 0.812 | 0.854 | 0.800 | 0.840 |
| 75–99 | 0.940 | 0.978 | 0.918 | 0.960 |
| 100–149 | 0.978 | 0.989 | 0.977 | 0.987 |
| 150–199 | 0.996 | 0.995 | 0.996 | 0.994 |
| ≥200 | 0.9995 | 0.998 | 0.9998 | 0.998 |

The curve is essentially the same in every column: ~46–48 % at DP 10–14, rising monotonically to
effectively 100 % by DP 150. **Phasing is depth-limited at the level of the individual site.**
The `coverage-per-variant` figure in `HAPCUT2.log` (32.7 / 36.8) averages over all of this and
hides it.

Note the minimum DP in every callset is 10 — both callers apply a min-depth-10 filter, so the
`<10` bin is empty by construction, not by chance.

### Finding 2 — naive beats bsgenova because of *where* its calls sit on that curve, not because its calls phase better

At equal depth the two callers phase almost identically (naive is ahead by only 2–6 points per
bin). What differs is the depth distribution of the sites each caller emits:

| | d1 bsgenova | d1 naive | d2 bsgenova | d2 naive |
|---|---|---|---|---|
| het sites offered | 3,897,860 | 4,389,230 | 4,203,619 | 4,885,323 |
| share at DP < 20 | 43.4 % | 20.0 % | 40.2 % | 18.3 % |
| median DP (callset) | 22 | 27 | 23 | 28 |

bsgenova puts ~40 % of its het calls in the DP 10–19 range where phasing is a coin flip; naive
puts under 20 % there. Standardising bsgenova's callset to naive's per-bin rates — i.e. asking
"what if bsgenova's sites phased as efficiently as naive's, depth for depth?" — recovers only:

- donor-1: 2,397,728 vs 2,239,973 actual = **+157,755**, which is 22.6 % of naive's 699,086 lead
- donor-2: 2,660,255 vs 2,459,141 actual = **+201,114**, which is 22.9 % of naive's 878,919 lead

So roughly **23 % of naive's advantage is per-site phasing efficiency and 77 % is having more
het calls, at greater depth**. The two donors agree to within 0.3 points, which is a tighter
replication than anything else in this experiment.

This also gives the mechanism behind the "caller ranking tracks the regime, not the donor" rule
recorded on 08-31: whichever caller's depth distribution best overlaps the phasable range wins,
and that overlap changes with assay and depth. It is not a property of the callers in the
abstract.

### Finding 3 — there are no runs of homozygosity, so they cannot explain anything

The hypothesis queued on 08-27 was that homozygous stretches might explain the gaps between
phase blocks. They do not, because there are no such stretches.

Window classification at 100 kb over chr1–22,X,Y (30,894 windows), donor-1 bsgenova, with the
other five runs within a factor of two on every line:

| class | windows | % genome | phasing rate |
|---|---|---|---|
| normal | 28,844 | 93.4 % | 0.572 |
| low_het (het frac < 25 %) | 192 | 0.62 % | 0.356 |
| roh (called, zero het) | 13 | 0.04 % | n/a |
| blind (< 5 calls) | 1,845 | 5.97 % | 0.298 |

The 13 `roh` windows are 13 *isolated* single windows — the longest run of homozygosity in the
whole genome is 100 kb, and each such window holds only 5–24 calls. Real ROH in a human genome
is multi-Mb. There is nothing here.

Two checks that this is not an artefact of the threshold or of het false positives:
- Relaxing the definition from "zero het" to "het fraction below 25 %" still captures only
  0.6 % of the genome (~18 Mb), non-contiguous. The conclusion is threshold-robust, which
  matters because a modest per-site het false-positive rate would erase short true ROH under
  the strict definition.
- The 5.97 % `blind` fraction is 184.5 Mb, which is about what hg38's assembly N-gaps and
  centromere models come to. Blind windows are the reference's gaps, not lost data.

### Finding 4 — the unphased sites are not in gaps; they are interleaved inside the blocks

Because the `--hic 1` blocks are Mb- to chromosome-scale, they overlap heavily, and the phased
skeleton leaves almost no genomic gaps at all:

| | d1 bsgenova | d1 naive | d2 bsgenova | d2 naive |
|---|---|---|---|---|
| unphased het sites | 1,657,887 | 1,450,171 | 1,744,478 | 1,547,263 |
| …inside the span of a block | 99.77 % | 99.82 % | 99.79 % | 99.81 % |
| …in an inter-block gap | 0.16 % | 0.12 % | 0.15 % | 0.13 % |
| …beyond the outermost block | 0.07 % | 0.06 % | 0.06 % | 0.06 % |
| inter-block gaps | 1,369 | 1,180 | 1,426 | 1,232 |
| total gap span | 37.8 Mb | 35.5 Mb | 36.5 Mb | 32.3 Mb |

Over 99.7 % of unphased het sites sit *inside* the span of a block that failed to incorporate
them. Only ~2,600 sites genome-wide fall in a true gap between blocks. Taken with Finding 1,
the picture is that the blocks are a sparse skeleton threaded across the whole chromosome, and
the sites they skip are skipped for lack of a read-backed link at their own low depth — not
because they lie in some unreachable region.

This closes the 08-27 question in the negative and reframes it: **there is no "gap" problem to
solve. There is a per-site depth problem.** The practical consequence for any future scale-up is
that adding cells helps by lifting individual sites over the DP threshold, which is consistent
with the repeated observation that more depth buys *more* Mb-scale blocks rather than longer
ones.

### One caveat found along the way — extreme-depth pileups in satellite regions

The naive callsets have a mean DP of ~113 against bsgenova's ~26, while the *medians* are close
(27 vs 22). The mean is dragged by a small number of windows with astronomical apparent depth,
all in centromeric/satellite sequence — chr3:93.4 Mb has a window mean DP of 41,924, and
chr16:46.3 Mb, chr2:32.9 Mb, chr16:34.5 Mb, chr5:49.6 Mb and chrY:56.7 Mb follow. These are
collapsed repeats: reads from many repeat copies stacking on one locus, which then look
heterozygous.

Sites at DP ≥ 200 are phased at ~100 % essentially by construction, so it is worth knowing how
much of the result they account for. They do not account for much: 35,808 of naive's 2,939,059
phased sites in donor-1 (1.2 %) and 7,859 of bsgenova's 2,239,973 (0.35 %). The excess is 27,949
sites, about **4 % of naive's 699,086 lead**. So the caller comparison is not driven by repeat
artefacts and Finding 2 stands — but a depth ceiling (say, drop DP > 200) would be a cheap and
defensible filter on any future run, and would remove a class of sites that are phased with
high confidence and probably wrong.

### Definitional note on the phased counts

`phasing_vs_depth.py` joins on `(chrom, pos)`, so at a multi-allelic position split across
several VCF records every record at a phased position counts as phased. HapCUT2's own totals are
per record. This makes the numbers here ~0.16 % higher than the runner's `SUMMARY` lines
(donor-1 bsgenova: 2,239,973 here vs 2,236,344 logged). It is a definitional difference, not a
disagreement about the data; the het-site totals match the log exactly (3,897,860 for donor-1
bsgenova, and so on for all four).

Separately, HapCUT2 preprocessing drops 0.55 % (bsgenova) to 0.96 % (naive) of callset het sites
before phasing, which is why the callset het counts here are slightly above the
`het_into_HapCUT2` figures in the runner log.

### Donor-4 pre-registration

Donor-4's merged BAM has 780,730,691 reads against donor-2's 868,532,873 — 10 % fewer — and its
callset medians are correspondingly lower (bsgenova 21 vs 23, naive 26 vs 28). Since Finding 1
says phasing rate is a stable function of DP, donor-4 should phase a somewhat *smaller* fraction
than donor-2 at both callers. Recording that here before the run lands, as a check on whether
the DP curve actually predicts out of sample.

Donor-4 also breaks the raw-count ordering for the first time: bsgenova 4,995,524 raw vs naive
4,884,362, where donors 1 and 2 both had naive ahead by 6–9 %. On het counts naive is still
narrowly ahead (4,134,601 vs 4,029,953).

## 2026-09-02 05:37 UTC — TASK B COMPLETE: donor-4 done, three donors in both assays

H1930004's snM3C 100-cell run finished, closing the scale-up. Timings tracked donor-2 closely and
came in slightly faster throughout (781 M reads vs 869 M): DL 3h05m, coord-sort 10m, merge 11m
(780,730,691 reads / 40,865,821,525 B, count exact), repair 3h38m, genotyping 1h57m, preprocessing
8m, `--hic 1` phasing 21h08m. Merged BAM published clean on attempt 1.

### snM3C-seq 100-cell, `--hic 1`, all three donors

| donor | caller | raw | het | blocks | phased | max span | blocks >= 1 Mb |
|---|---|---|---|---|---|---|---|
| H1930001 | bsgenova | 4,939,635 | 3,897,860 | 715,422 | 2,236,344 | — | — |
| H1930001 | naive | 5,242,025 | 4,389,230 | 825,583 | 2,937,592 | — | — |
| H1930002 | bsgenova | 5,278,393 | 4,203,619 | 784,565 | 2,455,452 | 141,135,646 | 93,133 |
| H1930002 | naive | 5,771,712 | 4,885,323 | 930,226 | 3,336,467 | 163,890,201 | 139,511 |
| H1930004 | bsgenova | 4,995,524 | 4,013,332 | 741,599 | 2,347,656 | 150,102,174 | 104,017 |
| H1930004 | naive | 4,884,362 | 4,093,725 | 771,136 | 2,707,575 | 174,816,698 | 133,633 |

**naive wins in all three donors on snM3C** (+31 %, +36 %, +15 % phased over bsgenova), while
**bsgenova won in both donors on 200-cell snMC** (6.5x, 2.6x). Three donors, two assays, no
exceptions: the caller ranking is a property of the assay/depth regime, not the sample.

One detail worth remembering because it looks alarming and is not: **donor-4's RAW counts invert**
(bsgenova 4,995,524 > naive 4,884,362, the only donor where that happens) but the inversion
disappears at the het stage (naive 4,093,725 > bsgenova 4,013,332), because naive emits a higher
heterozygous fraction (83.8 % vs 80.3 %). Raw counts are a poor predictor of phasing yield —
donor-2's naive callset was 9 % larger than its bsgenova one yet phased 36 % more SNPs.

Donor-4 also has the **largest max spans of any donor** (150 Mb / 175 Mb) despite the smallest read
count, which is further evidence for the sparse-skeleton reading: span is set by long-range Hi-C
links, not depth.

### New mount lesson: `readdir` can omit a file that is present and byte-correct

A file-count audit reported donor-4's bsgenova `.fragments` as unpublished — it was absent from
`ls` and from `find`, in a directory whose other four files were listed normally. It was not
missing. `stat` on the full path returned the correct size (671,406,491 B, matching local exactly),
the file read fine, and `md5sum` matched the local copy
(`d8ec6086288cfea15a229393ca3826bd`). **After the full read, the entry appeared in `ls`.**

So this mount can serve a stale directory listing that omits a real entry, and a *full read of the
path* refreshes it. This is the mirror image of the 2026-08-25 lesson ("an immediate post-write
`stat` is not a valid check"): there, a listing showed wrong metadata for a file that was fine; here,
a listing omits a file that is fine.

**Audit rule for this mount: never conclude a file is missing from `ls`/`find` alone.** Confirm with
`stat` on the full path, and prefer a size or md5 check against the local source. A name-based diff
is the weakest possible evidence here — this is now the third distinct way it has misled
(sibling-tree publication, caller-name renaming on publication, and now readdir omission).

### Final state

25 published files for each of H1930002 and H1930004 in each assay (H1930001 has 28 in snM3C,
carrying the extra `published-2/` and repaired-BAM artefacts from the 08-26 recovery). All twelve
snM3C callsets/phasings are distinct by md5 across donors.

**The scale-up planned on 2026-08-27 is complete: three donors, both assays, both callers.**
Item 2 of that plan, the coverage and homozygosity analysis, was started on 09-01 (entry above) and
is finished in the 09-02 16:05 entry below. It carried an additional motivating question: why donor-4
yields 81 % more bsgenova snMC calls than donor-2 on 46 % more reads, and why donor-2's naive snMC
callset (410,865) is such an outlier. Both are answered in the snMC-seq log's 09-02 entry.

## 2026-09-02 16:05 UTC — donor-4 settles the pre-registration: the DP curve transfers across donors, not across regimes

Donor-4's phasing published at 05:33-05:37, so the two runs the 09-01 analysis had to skip were
re-run (`run_covhom_jobs.sh`, ~2 min, fenced to cores 0-63). Nothing else changed; donors 1 and 2
are unaltered. The completed analysis spans 30 runs across the three experiments, all error-free.

| | d4 bsgenova | d4 naive |
|---|---|---|
| het sites offered | 4,013,332 | 4,093,725 |
| phased | 2,351,363 | 2,708,901 |
| phasing rate | 0.5859 | 0.6617 |
| share of het calls at DP < 20 | 46.4 % | 22.2 % |
| unphased sites inside a block's span | 99.80 % | 99.80 % |
| inter-block gaps / total gap span | 1,394 / 36.4 Mb | 1,193 / 31.3 Mb |

**Findings 2, 3 and 4 all replicate in the third donor.** Standardising donor-4's bsgenova callset
to naive's per-bin rates recovers 67,341 sites, **18.8 %** of naive's 357,538-site lead, against
22.6 % (donor-1) and 22.9 % (donor-2) — so ~80 % of naive's advantage is again where its calls sit
on the depth curve, not how well they phase. 99.80 % of unphased het sites are inside a block span.
The longest run of homozygosity is 200 kb (`roh` = 8 windows for bsgenova, 21 for naive).

### The pre-registered prediction was half right, and the half that failed is the informative half

09-01 predicted, before these numbers existed, that donor-4 would phase a **smaller fraction than
donor-2 at both callers**, because its BAM is 10 % smaller and its callset medians are 2 lower.

- naive: **0.6617 vs donor-2's 0.6833.** Correct.
- bsgenova: **0.5859 vs 0.5850.** Wrong — flat, in fact a hair higher, on a shallower callset.

Applying donor-2's per-bin rates to donor-4's depth distribution explains why: it under-predicts
bsgenova by 2.30 % and over-predicts naive by 1.71 %. Donor-4's bsgenova sites phase slightly
*better* at equal depth than donor-2's, which is just enough to cancel its depth deficit. Using
donor-1 as the reference instead gives -2.72 % and +0.49 %, and predicting donor-2 from donor-1
gives -0.36 % and -1.14 %.

**So the depth curve is donor-portable to within ~1-3 %** — good enough to budget a run, not good
enough to call a 0.2-point difference between two donors. The residual is real and caller-specific,
and it is the reason a "smaller BAM must phase worse" argument fails at this margin.

### The same curve does not survive a change of assay

Taking donor-2's snM3C per-bin rates and applying them to donor-2's own 200-cell snMC callset —
same donor, same caller, same genome, different assay and no `--hic 1` — predicts 412,028 phased
sites against 137,155 actual, an error of **+200 %**.

That is the boundary condition on Finding 1, and it should be stated with the finding from now on:
**per-site DP predicts phasing well within a fixed (assay, mode, depth) regime and not at all
across regimes.** The snMC-seq log's 09-02 entry works this out on the snMC side, where the picture
is different in every respect that matters — including Finding 4, which inverts outright.

### Published

`analyses/coverage-homozygosity-2026-09-02/` — `summaries/` holds the nine summary tables
uncompressed; `covhom_outputs_2026-09-02.tar.gz` (43,786,820 B, md5 `3778a152dc7bb37483f7cb46f2895540`,
verified by full read after the write) holds every per-run table plus the runners and job lists.
`predict.py` and `mktable.py` joined `workspace/Python-scripts/`, and the generalised runner is
`workspace/shell-scripts/run_covhom_jobs.sh`.

**Both items of the 08-27 plan are now complete.**

### Mount lesson: `>>` DESTROYS the file — never append in place on this mount

Writing this entry with `cat >> snM3C-seq_SNP_phasing_experiment.md` **silently destroyed the whole
log.** The append landed at the correct offset (125,033) and the new text was written intact, but
**every byte before it had been replaced by NULs.** The file was still the right size, `ls` looked
normal, and only `wc -l` reporting 60 lines gave it away.

Recovered in full from a `/tmp` copy taken before the edit, plus `tr -d '\0'` on the wreckage to
pull the new entry back out; the restored file differs from the backup only by the three lines this
session intended to reword. Both logs now match their local copies by md5.

**Rule: build the complete file locally in `/tmp`, then `cp` it over the mount copy, then verify by
md5.** Full-file writes work — the `cp` and the in-place Python rewrites in this session all
verified clean, and so did the 43 MB tarball. It is specifically `O_APPEND` that is broken here.
Take a `/tmp` backup before touching a file on the mount; without one this log would have been
unrecoverable beyond the 08-18 commit.

This is the fourth distinct way this mount has misled: post-write `stat` (08-25), publication under
a sibling tree and renamed callers (08-26/27), `readdir` omitting a byte-correct file (09-02), and
now silent NUL-fill on append.

## 2026-09-02 17:00 UTC — complete snM3C-seq results table, and housekeeping

The 10-cell snM3C callsets had never been through the coverage analysis, so they were run today
(6 runs, `-hicfix` blocks, i.e. the `--hic 1` phasing). That completes the matrix: **every
snM3C-seq callset in this project now has measured raw / het / offered / phased counts from one
consistent tool**, rather than counts assembled from different runner logs over five weeks.

### snM3C-seq, all donors, both callers, `--hic 1`

| donor | cells | caller | raw calls | het calls | het offered to HapCUT2 | phased | phasing rate |
|---|---|---|---|---|---|---|---|
| H1930001 | 10 | bsgenova | 186,988 | 165,413 | 163,894 | 48,239 | 0.2943 |
| H1930001 | 10 | naive | 73,089 | 61,426 | 60,099 | 20,192 | 0.3360 |
| H1930002 | 10 | bsgenova | 302,301 | 263,936 | 262,266 | 81,033 | 0.3090 |
| H1930002 | 10 | naive | 123,590 | 101,713 | 99,905 | 31,374 | 0.3140 |
| H1930004 | 10 | bsgenova | 230,426 | 199,739 | 198,106 | 61,774 | 0.3118 |
| H1930004 | 10 | naive | 93,657 | 77,093 | 75,502 | 23,809 | 0.3153 |
| H1930001 | 100 | bsgenova | 4,939,635 | 3,919,539 | 3,897,860 | 2,239,973 | 0.5747 |
| H1930001 | 100 | naive | 5,242,025 | 4,431,872 | 4,389,230 | 2,939,059 | 0.6696 |
| H1930002 | 100 | bsgenova | 5,278,393 | 4,227,826 | 4,203,619 | 2,459,141 | 0.5850 |
| H1930002 | 100 | naive | 5,771,712 | 4,930,636 | 4,885,323 | 3,338,060 | 0.6833 |
| H1930004 | 100 | bsgenova | 4,995,524 | 4,029,953 | 4,013,332 | 2,351,363 | 0.5859 |
| H1930004 | 100 | naive | 4,884,362 | 4,134,601 | 4,093,725 | 2,708,901 | 0.6617 |

Two things the 10-cell rows add to the record:

- **bsgenova's low-depth advantage is a callset-size advantage, not a phasing advantage.** At 10
  cells it calls 2.4-2.6x more het sites than naive in every donor, but its phasing *rate* is lower
  in all three (0.294/0.309/0.312 against 0.336/0.314/0.315). The 08-27 rule stands — bsgenova is
  the caller to use at low depth — but the reason is that it finds sites naive discards, not that
  those sites phase better.
- **Going 10 -> 100 cells roughly doubles the phasing rate in every donor and caller** (0.29-0.34 to
  0.57-0.68), and it is at 100 cells that the ordering flips to naive. Both tiers are `--hic 1`, so
  this is depth alone, with linkage held constant.

### Why "het calls" and "het offered" differ, in both directions

They are not the same quantity and the gap goes both ways, which had not been pinned down before:
HapCUT2 preprocessing **drops** a fraction of het sites (0.2-1.0 %), and it **splits multi-allelic
sites across several records**. Measured on H1930002's 200-cell bsgenova callset: 847,372
preprocessed records at 842,738 unique positions (4,634 positions duplicated), every one of them
het in the callset — none from hom-ALT, none from outside the callset — against 844,504 het
positions called. So the record count can exceed the het-call count while sites are still being
lost. Phased counts in these tables join on `(chrom, pos)` and therefore run ~0.16 % above the
per-record `SUMMARY` totals in the runner logs; the earlier numbers are not wrong, they count
differently.

### Housekeeping

- The analysis bundle moved from `analyses/coverage-homozygosity-2026-09-02/` to
  **`workspace/coverage-homozygosity-computation/`**, and now includes the 10-cell snM3C set:
  36 runs in four sets, twelve summary tables under `summaries/`, everything else in
  `covhom_outputs_2026-09-02.tar.gz` (47,571,925 B, md5 `5e0d479abde2907846758f00a58add0b`,
  verified by full read after the write, and re-verified after the old copy was deleted).
- **`workspace/data/README.md`** now documents the whole `data/` hierarchy: what each experiment
  folder holds, what every file extension means, which BAMs were published and which deliberately
  were not, and the two mount cautions.
- `/tmp` is clear. Nothing is running on the node.
