# Single-Cell Data Interoperability: H5AD to Seurat

This workflow addresses the common challenge of moving single-cell data between Python (Scanpy) and R (Seurat) ecosystems while preserving critical metadata.

## 1. Prerequisites

To read `.h5ad` files directly in R, you must have a Python environment with the `anndata` library installed, which R accesses via the `reticulate` bridge.

* **Python Library:** `anndata`
* **R Packages:** `Seurat`, `reticulate`, `anndata` (R wrapper)
* **Cypress Environment:** Use the previously configured `R_env`.

---

## 2. Conversion Workflow

The conversion involves importing the Python `anndata` module into R, reading the file as an AnnData object, and then manually reconstructing a Seurat object.

### Conversion Script (`h5ad_to_seurat.R`)

```r
library(Seurat)
library(reticulate)
library(anndata)

# 1. Point to the specific Conda environment containing Python/anndata
use_condaenv("R_env", required = TRUE)

# 2. Import the Python module and read the dataset
ad <- import("anndata")
adata <- ad$read_h5ad("/lustre/project/sglover3/merge_files/gingiva_healthy_data/Mucosal_Atlas_Bryd_10X/eb4d0514-78b7-4328-905c-d379918eeeae.h5ad")

# 3. Create Seurat Object
# Note: adata$X is transposed (Cells x Genes) in Python; 
# Seurat requires (Genes x Cells), hence the t() function.
S1 <- CreateSeuratObject(
  counts = t(adata$X),
  meta.data = as.data.frame(adata$obs)
)

# 4. Persistence
save(S1, file='Mucosal_Atlas_Bryd_10X.RData')

```

---

## 3. Data Mapping Details

When converting manually, it is important to understand how the components translate:

| AnnData Component | Seurat Component | Description |
| --- | --- | --- |
| `adata$X` | `counts` | The raw or normalized expression matrix. |
| `adata$obs` | `meta.data` | Cell-level annotations (e.g., cell type, batch). |
| `adata$var` | `features` | Gene-level metadata. |
| `adata$obsm` | `reductions` | Embeddings like PCA or UMAP (not included in the basic script). |

---

## 4. Technical Considerations

* **Memory Usage:** Loading large `.h5ad` files via `reticulate` can be very RAM-intensive as the data is essentially duplicated in memory during the transition from Python to R. Ensure you are on a high-memory compute node (e.g., via `idev`).
* **Sparse vs. Dense:** If `adata$X` is stored as a sparse matrix in Python, `reticulate` usually handles this correctly, but you may occasionally need to cast it using `Matrix::Matrix(..., sparse = TRUE)` to maintain efficiency.
* **Transposition:** Python's AnnData stores observations as rows (), while Seurat expects features as rows (). The `t()` function is mandatory.
