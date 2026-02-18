This documentation covers the integration of **Honeycomb Bio BeeNet** for processing scRNA-seq data and the subsequent workflow for merging these datasets into a unified Seurat object.

---

# Single-Cell RNA-Seq: BeeNet Plus & Seurat Integration

This guide details the transition from BeeNet-processed FASTQ files to a merged, multi-sample Seurat object on the Cypress HPC.

## 1. BeeNet Processing (`analyze`)

BeeNet is utilized for processing Honeycomb Bio's unique single-cell data formats. The `analyze` command performs alignment and generates a **Transcript Count Matrix (TCM)**.

### Command Structure

We utilize the project storage path for the `beenet` binary and reference genomes. Note the use of `--num-barcodes` to refine the expected cell count.

```bash
# Basic Analysis
/lustre/project/sglover3/beenet analyze \
  --sample-name=290C1 \
  --out=out_290C1_1 \
  --ref=20210603_GRCh38.104 \
  --num-barcodes=25000 \
  transfer_WALKUP-19750_..._290C1_*fastq.gz

```

### Key Parameters

* **`--sample-name`**: Unique identifier for the output.
* **`--ref`**: Reference genome directory (GRCh38.104).
* **`--num-barcodes`**: Limits processing to the top N barcodes (useful for filtering background noise).
* **Output**: Generates a `.tsv.gz` file (Transcript Count Matrix) in the specified output directory.

### Downloading the Reference Genome

BeeNet provides a built-in utility to download and prepare reference genomes. On Cypress, we store these in the shared project directory to ensure accessibility and save user home directory space.

### Command

```bash
# Download the GRCh38 human reference (version 104)
/lustre/project/sglover3/beenet download-ref 20210603_GRCh38.104

```

> [!NOTE]
> The reference download and indexing process is a one-time setup step. Ensure you have sufficient disk space in `/lustre/project/sglover3/` before starting, as reference genomes typically exceed 20GB.

---

## 2. Downstream R Analysis: Loading BeeNet TCMs

BeeNet outputs data in a Tab-Separated Value (TSV) format, which requires a different loading approach than standard Cell Ranger 10X outputs.

### Single Sample Initialization

```r
library(Seurat)
library(data.table)
library(Matrix)
library(tidyverse)

# Define paths to BeeNet .tsv.gz outputs
sample_paths <- c(
  b1 = "/lustre/project/sglover3/.../006_20251203123233_TCM.tsv.gz",
  ...
)

# Load TSV and convert to matrix
counts_data <- read.delim(sample_paths[1], sep = "\t", header = TRUE, row.names = 1)
counts_matrix <- as.matrix(counts_data)

# Create Seurat Object
S1 <- CreateSeuratObject(counts = counts_matrix, project = "BeenetB1", min.cells = 3, min.features = 200)
save(S1, file='beenetb1.RData')

```

---

## 3. Multi-Sample Merging

To analyze across biological replicates or conditions, multiple Seurat objects are combined into a single object.

### The Merge Workflow

The `merge` function allows for the aggregation of multiple objects while appending unique cell IDs to prevent name collisions across samples.

```r
# Assuming S1, S2, S3, S4 have been created similarly to S1 above
MS <- merge(S1, 
            y = c(S2, S3, S4), 
            add.cell.ids = c("B1", "B2", "B3", "B4"), 
            project = "CombinedBeenet")

# Save merged object for future analysis
save(MS, file='beenetmerged.RData')

```

### Data Management Table

| Object | Source Path | Cell ID Prefix |
| --- | --- | --- |
| `S1` | `.../006_..._TCM.tsv.gz` | B1 |
| `S2` | `.../0061_..._TCM.tsv.gz` | B2 |
| `S3` | `.../0062_..._TCM.tsv.gz` | B3 |
| `S4` | `.../1921_..._TCM.tsv.gz` | B4 |

---

## 4. Technical Notes

* **Memory Constraints:** Loading large TSVs into a standard R matrix via `as.matrix()` can be RAM-intensive. For very large datasets, consider using `fread()` from `data.table` and converting directly to a **Sparse Matrix** using `as(as.matrix(data), "dgCMatrix")` to save memory.
* **Storage:** `.RData` files are compressed by default, but for large merged objects, `saveRDS()` is often preferred for better compatibility with single-object loading (`readRDS()`).
