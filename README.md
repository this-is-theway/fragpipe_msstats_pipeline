# FragPipe + MSstats Proteomics Pipeline

A [Nextflow](https://www.nextflow.io/) (DSL2) pipeline for label-free quantitative proteomics:
raw MS files are searched with **FragPipe**, and the resulting quantification table is analyzed
statistically with **MSstats** and turned into a self-contained **HTML report** (RMarkdown).

```
 .raw/.mzML/.d files ──► (RAW_TO_MZML) ──► (GENERATE_MANIFEST) ──► FRAGPIPE ──► MSSTATS_REPORT ──► msstats_report.html
```

---

## Contents

- [Requirements](#requirements)
- [Repository layout](#repository-layout)
- [Input files](#input-files)
- [Usage](#usage)
- [How it works](#how-it-works)
- [Parameters](#parameters)
- [Output](#output)
- [Containers](#containers)
- [Troubleshooting / notes](#troubleshooting--notes)

---

## Requirements

| Requirement | Notes |
|---|---|
| [Nextflow](https://www.nextflow.io/) ≥ 23.04 | pinned in `manifest.nextflowVersion` in `nextflow.config` |
| Java 11+ | required by Nextflow |
| [Apptainer/Singularity](https://apptainer.org/) (or Docker) | processes are containerized via the `container` directive; build the images from `containers/*/*.def` (Apptainer) or `containers/*/Dockerfile` (Docker) before running |
| A local **FragPipe** installation | the `FRAGPIPE` process calls the FragPipe headless binary directly (path set by `params.fragpipe_bin`), it is not bundled inside a container by default |
| Python 3 | used by `GENERATE_MANIFEST` and by FragPipe itself (`params.python_bin`) |
| (optional) [Slurm](https://slurm.schedmd.com/) | only needed if you run with `-profile slurm` |

R/Bioconductor dependencies for the MSstats/report step (MSstats-family package, `tidyverse`,
`ggplot2`, `ggrepel`, `DT`, `openxlsx`, `rmarkdown`, `PCAtools`, `Biostrings`, …) are installed
inside the `containers/msstats` image — see [Containers](#containers) — so no local R setup is
required as long as containers are enabled.

---

## Repository layout

```
main.nf                      Pipeline entry point / workflow logic
nextflow.config               Default params, process resources, profiles, container map
modules/
  fragpipe.nf                 RAW_TO_MZML, GENERATE_MANIFEST, FRAGPIPE processes
  msstats.nf                  MSSTATS_REPORT process
report/
  run_msstats.Rmd              Main RMarkdown report template
  scripts/                      R scripts sourced by the report (load, preprocess, QC, differential expression, annotation)
  img/logo.png                  Default report logo
containers/
  fragpipe/, msstats/, thermorawfileparser/, diann/   Dockerfile / Apptainer .def build recipes
```

---

## Input files

| File | Required | Format |
|---|---|---|
| MS raw data (`--input`) | yes | Glob pattern, e.g. `'/data/raw/*.raw'`. Accepts `.raw` (auto-converted to `.mzML`), `.mzML`, `.mzXML`, or `.d` (Bruker timsTOF — enables ion-mobility/DIA-NN options in the FragPipe workflow automatically) |
| FASTA database (`--fasta`) | yes | Protein FASTA used for the search and (optionally) for annotation in the report |
| FragPipe workflow file (`--fragpipe_workflow`) | yes | A `.workflow` file exported from the FragPipe GUI. The pipeline appends the FASTA path and Python binary path to it automatically |
| Experimental design (`--experimental_design`) | yes (unless `--fragpipe_manifest` is given) | CSV with columns: `file, experiment, bioreplicate, data_type`. `file` must match (or match after stripping `.raw`→`.mzML`) a name under `--input`. Used to build the FragPipe `.fp-manifest` |
| Comparisons file (`--comparisons_file`) | yes | CSV with columns: `condition_A, condition_B` — one row per pairwise comparison to run in MSstats |
| Highlight genes (`--highlight_genes_file`) | no | Plain text, one gene symbol per line; these are highlighted in volcano plots |
| RMarkdown template (`--rmd_template`) | yes | Normally `report/run_msstats.Rmd` |
| R scripts dir (`--scripts`) | yes | Normally `report/scripts` |
| Logo (`--logo_path`) | yes | Image shown in the report header, e.g. `report/img/logo.png` |
| Pre-made manifest (`--fragpipe_manifest`) | no | Skip `GENERATE_MANIFEST` and supply your own `.fp-manifest` |
| Pre-computed `msstats.csv` (`--msstats_csv`) | no | Skip raw-file conversion + FragPipe entirely and go straight to the MSstats/report step (see [Usage](#usage)) |

---

## Usage

### 1. Build/obtain the containers

```bash
# Apptainer / Singularity (recommended on HPC/Slurm)
apptainer build containers/msstats/msstats.sif               containers/msstats/msstats.def
apptainer build containers/thermorawfileparser/thermorawfileparser_1.4.5.sif  containers/thermorawfileparser/thermorawfileparser.def
apptainer build containers/fragpipe/fragpipe_24.0.sif         containers/fragpipe/<def-or-recipe>

# or Docker, if you adapted nextflow.config to use Docker images instead
```

> Place the built `.sif` images where `params.container_dir` (default: `${projectDir}/containers`)
> expects them — see the [Troubleshooting](#troubleshooting--notes) note about container paths
> before running on a new machine.

### 2. Full run (raw files → report), local executor

```bash
nextflow run main.nf -profile local \
  --input                '/data/raw/*.raw' \
  --fasta                 /data/db/proteome.fasta \
  --fragpipe_workflow     /data/config/dda_lfq.workflow \
  --experimental_design   /data/design/experimental_design.csv \
  --comparisons_file      /data/design/comparisons.csv \
  --highlight_genes_file  /data/design/genes_to_highlight.txt \
  --rmd_template          report/run_msstats.Rmd \
  --scripts               report/scripts \
  --logo_path             report/img/logo.png \
  --outdir                results \
  --report_title          "My Proteomics Project" \
  --report_author         "Your Name"
```

### 3. Full run on a Slurm cluster

```bash
nextflow run main.nf -profile slurm \
  --slurm_queue        <queue-name> \
  --slurm_queue_size   50 \
  --slurm_cluster_options '--account=my_account' \
  --input ... --fasta ... --fragpipe_workflow ... --experimental_design ... \
  --comparisons_file ... --rmd_template report/run_msstats.Rmd --scripts report/scripts \
  --logo_path report/img/logo.png
```

> `slurm_queue`, `slurm_queue_size`, and `slurm_cluster_options` are referenced by the `slurm`
> profile but have no defaults in `nextflow.config` — you must supply them on the command line
> (as above) or in a custom config/params file.

### 4. MSstats-only mode (reuse an existing `msstats.csv`)

Skip raw-file conversion, manifest generation, and the FragPipe search entirely by passing
`--msstats_csv`:

```bash
nextflow run main.nf -profile local \
  --msstats_csv           results/fragpipe/fragpipe_output/msstats.csv \
  --fasta                 /data/db/proteome.fasta \
  --comparisons_file      /data/design/comparisons.csv \
  --highlight_genes_file  /data/design/genes_to_highlight.txt \
  --rmd_template          report/run_msstats.Rmd \
  --scripts               report/scripts \
  --logo_path             report/img/logo.png \
  --experimental_design   /data/design/experimental_design.csv
```

(`--experimental_design` and `--fragpipe_workflow` are still checked for existence by the
pipeline's parameter validation even in this mode unless provided as dummy/placeholder paths —
pass the files you actually used to produce the `msstats.csv`.)

### Useful generic options

```bash
nextflow run main.nf -profile local ... -resume      # resume from cache after a failed/partial run
nextflow run main.nf -profile local ... -with-dag dag.html
```

---

## How it works

1. **`RAW_TO_MZML`** *(per file, skipped for non-`.raw` input)*
   Converts vendor `.raw` files to `.mzML` using ThermoRawFileParser (`-f=2`, indexed mzML).
   Already-converted/other formats (`.mzML`, `.mzXML`, `.d`) pass through unchanged.

2. **`GENERATE_MANIFEST`** *(skipped if `--fragpipe_manifest` is supplied)*
   A small embedded Python script matches each `file` entry in `--experimental_design` against
   the staged MS files (handling the `.raw` → `.mzML` rename), and writes a tab-separated
   FragPipe `.fp-manifest` (`path  experiment  bioreplicate  data_type`).

3. **`FRAGPIPE`**
   Runs FragPipe headless (`${fragpipe_bin}/fragpipe --headless`) with:
   - the supplied `.workflow` file, augmented with the resolved FASTA path and
     `fragpipe-config.bin-python`
   - if any input is Bruker `.d` data, the workflow is further augmented to enable
     ion-mobility / DIA-NN spectral-library prediction
   - the generated/supplied manifest, `--ram` and `--threads` taken from
     `params.fragpipe_ram` / `params.fragpipe_threads`
   The quantitative output consumed downstream is `fragpipe_output/**/msstats.csv`
   (MSstats-format table that FragPipe emits).

4. **`MSSTATS_REPORT`**
   Renders `report/run_msstats.Rmd` with `rmarkdown::render()`, passing all analysis
   parameters (see below). The Rmd, in turn, sources the scripts in `report/scripts/` in order:
   - `01_load_data.R` — `MSstats::FragPipetoMSstatsFormat()` on the `msstats.csv`
   - `02_preprocessing.R` — `dataProcess()`: log2 transform, median normalization,
     top-N feature selection, Tukey Median Polish (TMP) summarization, model-based
     imputation for censored values
   - `03_quality_control.R` — per-sample intensity boxplots and a PCA (via `PCAtools`) on the
     protein-level matrix
   - differential expression, for each row of `--comparisons_file`:
     - `imputation_method = "AFT"` → `04_differential.R` + `comparison_template.Rmd`
     - `imputation_method = "ludovic"` (default) → `05_differential_ludovic.R` +
       `comparison_template_ludovic.Rmd`
     Each comparison gets fold-change/p-value testing with FDR correction, a volcano plot
     (with optional gene highlighting from `--highlight_genes_file`), and FASTA-based
     annotation.
   The final self-contained `msstats_report.html`, plus per-comparison `.xlsx`/`.csv` result
   tables and logs, are published to `--outdir`.

Nextflow's built-in `report`, `timeline`, `trace`, and `dag` outputs are also written under
`<outdir>/pipeline_info/` for every run.

---

## Parameters

### Input / Output

| Parameter | Default | Description |
|---|---|---|
| `--input` | *(none — required)* | Glob pattern to MS files, e.g. `'/data/raw/*.raw'` |
| `--fasta` | *(none — required)* | Protein FASTA database |
| `--experimental_design` | *(none — required unless `--fragpipe_manifest` given)* | CSV: `file, experiment, bioreplicate, data_type` |
| `--outdir` | `results` | Output directory |

### FragPipe

| Parameter | Default | Description |
|---|---|---|
| `--fragpipe_workflow` | *(none — required)* | Path to a FragPipe `.workflow` file |
| `--fragpipe_manifest` | `null` | Optional pre-made `.fp-manifest`; if set, skips `GENERATE_MANIFEST` |
| `--fragpipe_threads` | `8` | Threads passed to FragPipe (`--threads`), also used for `process_high` CPU allocation |
| `--fragpipe_ram` | `20` | RAM (GB) passed to FragPipe (`--ram`), also used for `process_high` memory allocation |
| `--fragpipe_bin` | `/fragpipe_bin/fragpipe-24.0/fragpipe-24.0/bin` | Directory containing the local FragPipe executable |
| `--python_bin` | `python3` | Python binary written into the workflow's `fragpipe-config.bin-python` |

### MSstats / R analysis

| Parameter | Default | Description |
|---|---|---|
| `--msstats_csv` | `''` | If set, skips raw-file conversion and FragPipe entirely and uses this pre-computed MSstats-format CSV |
| `--rmd_template` | *(none — required)* | Path to `report/run_msstats.Rmd` |
| `--scripts` | *(none — required)* | Path to the `report/scripts` directory |
| `--comparisons_file` | *(none — required)* | CSV: `condition_A, condition_B`, one row per pairwise comparison |
| `--highlight_genes_file` | `''` | Optional one-gene-symbol-per-line text file to highlight in volcano plots |
| `--imputation_method` | `ludovic` | Differential-expression path to use: `ludovic` (uses `05_differential_ludovic.R` / `comparison_template_ludovic.Rmd`) or `AFT` (uses `04_differential.R` / `comparison_template.Rmd`); also used as the imputation model label shown in the report |
| `--fdr_threshold` | `0.05` | Adjusted p-value (FDR) cutoff for significance |
| `--log2fc_threshold` | `2` | Log2 fold-change cutoff for significance |

### Report

| Parameter | Default | Description |
|---|---|---|
| `--report_title` | `Proteomics Analysis Report` | HTML report title |
| `--report_author` | `Author Name` | HTML report author / byline |
| `--logo_path` | *(none — required)* | Path to an image shown in the report header |

### Pipeline / resources

| Parameter | Default | Description |
|---|---|---|
| `--help` | `false` | (reserved; no help text is currently printed) |
| `--max_cpus` | `12` | Soft cap referenced for resource-capping customization |
| `--max_memory` | `48.GB` | Soft cap referenced for resource-capping customization |
| `--max_time` | `72.h` | Soft cap referenced for resource-capping customization |
| `--container_dir` | `${projectDir}/containers` | Base directory used to locate container images (see [note](#troubleshooting--notes)) |

### Slurm profile (`-profile slurm` only, no built-in defaults)

| Parameter | Description |
|---|---|
| `--slurm_queue` | Slurm partition/queue name |
| `--slurm_queue_size` | Executor `queueSize` (max concurrently submitted jobs) |
| `--slurm_cluster_options` | Extra `sbatch` options, e.g. `--account=...` |

### Process resource labels (`nextflow.config`)

| Label | cpus | memory | time |
|---|---|---|---|
| `process_low` | 2 | 4 GB | 1 h |
| `process_medium` | 12 | 48 GB | 8 h |
| `process_high` | `params.fragpipe_threads` | `params.fragpipe_ram` GB | 48 h |

---

## Output

```
<outdir>/
├── fragpipe/                       Raw FragPipe output (publishDir copy of fragpipe_output/**)
│   └── fragpipe_output/.../msstats.csv
├── msstats/
│   ├── msstats_report.html         Final self-contained HTML report
│   ├── *.xlsx                      Per-comparison / protein / peptide-level result tables
│   ├── *_files/                    Report assets (figures, etc.)
│   └── logs/                       MSstats / rendering logs (*.log)
└── pipeline_info/
    ├── execution_report_<timestamp>.html
    ├── execution_timeline_<timestamp>.html
    ├── execution_trace_<timestamp>.txt
    └── pipeline_dag_<timestamp>.html
```

---

## Containers

| Directory | Image | Used by |
|---|---|---|
| `containers/thermorawfileparser/` | ThermoRawFileParser 1.4.5 (via micromamba) | `RAW_TO_MZML` |
| `containers/msstats/` | `rocker/r-ver:4.6.0` + tidyverse/Bioconductor/`MSstatsAlt` | `MSSTATS_REPORT` (note: the process also hardcodes `container 'bioconductor/bioconductor_docker:devel'` directly in `modules/msstats.nf` — see note below) |
| `containers/fragpipe/` | *(recipe not included — add your own FragPipe build)* | `GENERATE_MANIFEST`, `FRAGPIPE` |
| `containers/diann/` | DIA-NN 2.x (Academia) | bundled for DIA-NN–based FragPipe workflows (DIA-NN is invoked by FragPipe itself for DIA data; not a separate Nextflow process in this pipeline) |

---

## Troubleshooting / notes

- **Container path mismatch**: in `nextflow.config`, the `withName` container paths are written
  with single-quoted Groovy strings (e.g.
  `container = '${params.container_dir}/containers/fragpipe/fragpipe_24.0.sif'`). Single-quoted
  strings in Groovy are **not** interpolated, so `${params.container_dir}` is taken literally
  rather than substituted, and the path also duplicates `containers/` on top of
  `params.container_dir` (which already points at the `containers/` folder). If your run fails
  trying to pull/mount a container at a path that looks like `${params.container_dir}/...`,
  either hardcode the real `.sif` paths there or switch those lines to double-quoted strings.
- **`MSSTATS_REPORT` container directive**: `modules/msstats.nf` sets
  `container 'bioconductor/bioconductor_docker:devel'` directly on the process, which overrides
  the `withName: 'MSSTATS_REPORT'` entry in `nextflow.config`. Build/point at the image in
  `containers/msstats/` if you want the pinned `MSstatsAlt`-based image to actually be used.
- **FragPipe is not fully containerized by default**: `FRAGPIPE` and `GENERATE_MANIFEST` call
  `${params.fragpipe_bin}/fragpipe` directly, so even when a container is attached, FragPipe
  itself must be installed/licensed and reachable at `params.fragpipe_bin` inside that
  container (or on the host, if running without containers).
- **DIA-NN / Bruker `.d` data**: if any staged input ends in `.d`, the pipeline automatically
  appends `diann.predict-spectral-library=true` and switches the workflow to ion-mobility mode —
  make sure your `.workflow` file is a DIA-compatible FragPipe workflow if you intend to use this.
