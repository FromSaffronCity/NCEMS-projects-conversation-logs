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
