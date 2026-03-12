# Batch Processing: 10x HDF5 (.h5) to Merged Seurat Object

This workflow automates the ingestion of multiple `.h5` files from a research directory (e.g., nasal mucosa data), initializing individual Seurat objects and merging them into a single analytical atlas.

## 1. The HDF5 Advantage

Unlike the standard three-file output (matrix, features, barcodes), the `.h5` format is:

* **Self-contained:** All data required for a Seurat object is in one file.
* **Faster I/O:** Optimized for reading large datasets on HPC file systems.
* **Standardized:** Commonly used for data sharing in NCBI GEO/SRA (e.g., PMID: 34354210).

---

## 2. Automated Pipeline Script

The following script utilizes `lapply` to iterate through a directory, ensuring consistent preprocessing (filtering) across all samples before merging.

### Batch Script (`h5_batch_merge.R`)

```r
library(Seurat)
library(tidyverse)

# 1. Directory Setup
data_dir <- "/lustre/project/sglover3/nose_data/PMID34354210_10X/"
h5_files <- list.files(path = data_dir, pattern = "\\.h5$", full.names = TRUE)

# 2. Iterative Object Creation
seurat_list <- lapply(h5_files, function(file_path) {
 
  # Extract sample ID from filename for project labeling
  sample_name <- tools::file_path_sans_ext(basename(file_path))
 
  # Read the H5 file (native Seurat function)
  counts <- Read10X_h5(file_path)
 
  # Create Seurat object with standard QC filters
  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
 
  return(obj)
})

# 3. Aggregation (Merging)
# Combining individual objects into one "Atlas" object
if (length(seurat_list) > 1) {
  merged_object <- merge(
    x = seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = sapply(seurat_list, function(x) Project(x))
  )
} else {
  merged_object <- seurat_list[[1]]
}

# 4. Persistence
save(merged_object, file = 'merged_nasal_atlas_PMID34354210.RData')

```

---

## 3. Workflow Details

### Data Loading

The `Read10X_h5()` function automatically parses the HDF5 hierarchy. If the file contains multiple genomes or modalities (e.g., Gene Expression + Antibody Capture), this function will return a list of matrices.

### Merging & Cell Identification

When merging, it is critical to use `add.cell.ids`. This appends the sample name to each barcode (e.g., `GSM4964341_AAACCCA...`), preventing overlap if identical barcodes exist in different samples.

| Parameter | Function |
| --- | --- |
| `min.cells = 3` | Excludes genes that are not expressed in at least 3 cells. |
| `min.features = 200` | Excludes "empty" droplets or low-quality cells with <200 genes. |
| `Project(x)` | Sets the `orig.ident` metadata column to the sample name. |

---

## 4. Resource Allocation for HPC

Batch merging multiple `.h5` files is memory-intensive. For the Nasal Atlas dataset:

* **Node Type:** Use an `idev` session or a high-memory SLURM partition.
* **Memory Requirement:** Approximately **128GB - 250GB RAM** depending on the total number of cells across all files.
