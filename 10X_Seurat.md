# Single-Cell RNA-Seq Pipeline: Cell Ranger to Seurat

This workflow describes the transition from raw FASTQ files to a filtered Seurat object on the Cypress HPC.

## 1. Cell Ranger Installation

For high-throughput genomics, Cell Ranger is installed in the project directory to handle large reference genomes and high I/O workloads.

* **Version:** 10.0.0 (Download from https://www.10xgenomics.com/support/software/cell-ranger/downloads#download-links)
* **Path:** `/lustre/project/sglover3/cellranger-10.0.0`
* **Reference Genome:** GRCh38-2024-A (Download from https://www.10xgenomics.com/support/software/cell-ranger/downloads#reference-downloads)

---

## 2. SLURM Job Submission (`cellranger count`)

The following script automates the alignment, filtering, and barcode counting. Note that we disable BAM file generation (`--create-bam false`) to save significant disk space and reduce runtime. Run this script using the command `sbatch count_job.sh`.

### Batch Script (`count_job.sh`)

```bash
#!/bin/bash
#SBATCH --partition=centos7
#SBATCH --qos=normal
#SBATCH --job-name=cellranger_count
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G

# Define Path
CELLRANGER=/lustre/project/sglover3/cellranger-10.0.0/cellranger

# Example: Processing Sample 232c2
$CELLRANGER count \
  --id=232c2 \
  --create-bam false \
  --fastqs=/lustre/project/sglover3/progress_report_analysis/LC_COVID_10X_DATA/232c2 \
  --sample=232CTVLT4_2_0455058144_TulaneC2 \
  --transcriptome=/lustre/project/sglover3/refdata-gex-GRCh38-2024-A \
  --localcores=16 \
  --localmem=240

```

> [!TIP]
> **Scaling Up:** For projects with many samples, use a **SLURM Job Array** or a `for` loop within the script to iterate through directories, rather than hardcoding individual blocks.

---

## 3. Downstream Analysis with Seurat

Once Cell Ranger completes, it produces a `raw_feature_bc_matrix`. We use the **Seurat** R package to aggregate these samples and perform initial QC.

### Data Integration & QC Script

The following R code loads multiple sample paths, merges them into a single object, and calculates mitochondrial metrics.

```r
library('Seurat')
library('tidyverse')

# 1. Define named paths for all samples
sample_paths <- c(
  lc1 = "/lustre/project/sglover3/232lc1/outs/raw_feature_bc_matrix",
  lc2 = "/lustre/project/sglover3/232lc2/outs/raw_feature_bc_matrix",
  ... # add additional samples
)

# 2. Load and create Seurat Object
data <- Read10X(data.dir = sample_paths)
S <- CreateSeuratObject(counts = data, project = "232l", min.cells = 3, min.features = 200)

# 3. Quality Control (Mitochondrial content)
# High MT% often indicates dying cells or poor library quality
S[["percent.mt"]] <- PercentageFeatureSet(S, pattern = "^MT-")

# 4. Visualization
# Generate Violin plots to inspect distribution of features and MT reads
VlnPlot(S, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Correlation between total counts and unique features
FeatureScatter(S, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  geom_smooth(method = 'lm')

```

### QC Metric Definitions

| Metric | Description |
| --- | --- |
| **nFeature_RNA** | Number of unique genes detected in each cell. |
| **nCount_RNA** | Total number of molecules detected within a cell. |
| **percent.mt** | Percentage of reads mapping to mitochondrial genes. |

---

## 4. Key Considerations

* **Memory Management:** Cell Ranger is memory-intensive. Always match `--localmem` in the command to the `#SBATCH --mem` request.
* **Parallelization:** Ensure `localcores` matches your `cpus-per-task` to prevent CPU oversubscription on the node.
* **Path Management:** Always use absolute paths in SLURM scripts to avoid "File not found" errors during batch execution.
