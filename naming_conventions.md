# Data Management & Naming Conventions

To manage large-scale single-cell datasets efficiently, we utilize a **Project-Centric** organizational structure. This system prevents data fragmentation and simplifies path management in SLURM and R scripts.

## 1. Directory Hierarchy

All project data resides within `/project_name/mucosal_atlas/`. This structure enforces a logical separation between raw inputs, processed outputs, and software environments.

```text
/project_name/mucosal_atlas/
├── bin/                        # Custom analysis scripts (R, Bash, Python)
├── data/
│   ├── raw/                    # Immutable source files (Read-Only)
│   │   ├── 10x/                # Raw FASTQs from 10x platforms
│   │   ├── beenet/             # Raw FASTQs from Honeycomb HIVE
│   │   └── external/           # Downloaded .h5ad or .h5 files
│   └── reference/              # Genome indices (GRCh38, etc.)
├── results/
│   ├── cellranger/             # Output from cellranger count
│   ├── beenet/                 # TCM files from BeeNet Plus
│   └── seurat_objects/         # Final processed .RData and .rds files
├── docs/                       # Project documentation and GitHub Wiki
└── logs/                       # Standard output/error files from SLURM

```

---

## 2. Naming Convention (Metadata-Rich)

Files should be named using a consistent, "slug-based" format. This allows for rapid filtering and pattern matching in R (e.g., using `grep` or `list.files`).

### The Standard Template:

`[YYYYMMDD]_[Tissue]_[Condition]_[Technology]_[SampleID].[ext]`

| Segment | Allowed Values (Examples) | Description |
| --- | --- | --- |
| **Date** | `20260312` | Year, Month, Day of processing. |
| **Tissue** | `gingiva`, `nasal`, `ileum` | The specific mucosal source. |
| **Condition** | `healthy` | Experimental state of the donor. |
| **Technology** | `10x`, `beenet`, `h5ad` | The platform used for sequencing/storage. |
| **SampleID** | `donor05`, `232c2` | Unique identifier for the specific sample. |

**Example:** `20260312_nasal_healthy_h5_GSM4964341.h5`

---

## 3. Implementation in Workflow

Standardization allows for "Clean Code" practices. Instead of hardcoding 20 individual file paths, we use globbing patterns:

```r
# Example: Batch loading all healthy nasal data from 10x
data_path <- "/project_name/mucosal_atlas/results/cellranger/"
nasal_files <- list.files(path = data_path,
                          pattern = "nasal_healthy_10x",
                          full.names = TRUE)

# The pipeline can now iterate through nasal_files automatically

```

---

## 4. Best Practices for Cypress HPC

* **Read-Only Raw Data:** Never modify files in `data/raw/`. Copy or symlink them if they need to be manipulated.
* **Avoid Spaces:** Never use spaces or special characters (like `#`, `&`, or `(`) in filenames. Use underscores `_` or hyphens `-`.
* **Symlinking Large References:** To save space, symlink the `reference/` folder to the central lab reference directory instead of duplicating the ~30GB files.
* **Log Rotation:** Keep your `logs/` folder clean. Use descriptive names for SLURM logs, e.g., `#SBATCH --output=logs/cr_count_%j.log`.
