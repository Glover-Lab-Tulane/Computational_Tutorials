# Cypress HPC Single-Cell Analysis Pipeline

This repository contains comprehensive documentation and configuration files for processing single-cell RNA-seq data on the **Cypress High Performance Computing (HPC)** cluster. It covers the end-to-end workflow from environment setup and raw data processing to interactive visualization.

## 🧬 Project Overview

The pipeline supports two primary data streams:

1. **10x Genomics:** Using Cell Ranger for alignment and quantification.
2. **Honeycomb Bio (HIVE):** Using the BeeNet Plus suite for TCM generation.

All downstream analysis is unified within **R/Seurat**, hosted on a containerized RStudio Server environment for high-memory interactive analysis.

---

## 📂 Documentation Modules

### 1. Environment & Infrastructure

Before running the pipeline, the base environment must be established on the Lustre parallel filesystem.

* **[Installation Guide (install_R.md)](install_R.md):** Steps for installing a custom Anaconda3 distribution in project space and initializing the shell.
* **[Environment Config (R_env.yml)](R_env.yml):** The YAML definition for the R environment, ensuring reproducible package versions (Seurat, Tidyverse, etc.).

### 2. Data Processing Pipelines

Documentation for converting raw FASTQ files into gene-expression matrices.

* **[10x Genomics Workflow (10X_Seurat.md)](10X_Seurat.md):** SLURM job submission for `cellranger count` and initial Seurat QC.
* **[BeeNet Plus Workflow (fastq_Seurat.md)](fastq_Seurat.md):** Downloading BeeNet references, generating TCMs, and merging multi-sample HIVE datasets.

### 3. Interactive Analysis

* **[RStudio Server on HPC (R_studio_server.md)](R_studio_server.md):** Advanced guide on deploying a Singularity container, mapping R kernels, and establishing SSH tunnels for browser-based analysis.

---

### Resource Allocation Summary

| Task | Suggested Partition | Recommended Memory |
| --- | --- | --- |
| **Cell Ranger** | `centos7` | 240GB - 250GB |
| **BeeNet Analyze** | `centos7` | 128GB+ |
| **Seurat (Merging)** | `idev` (Interactive) | 250GB+ |

---

## 🛠 Prerequisites

* Access to the **Cypress HPC** cluster.
* A valid project account under `/lustre/project/sglover3/`.
* Singularity 3.9.0 (available via `module load`).
