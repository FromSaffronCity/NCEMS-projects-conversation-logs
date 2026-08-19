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
