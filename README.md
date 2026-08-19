# Aging Mouse Heart -- Multi-Omic Single-Cell Analysis

Co-expression, cell-cell communication, and gene regulatory network analysis of aging in the
mouse heart, using single-cell RNA-seq data from *Tabula Muris Senis* (FACS + droplet).

The pipeline combines four complementary analyses run on the **same** annotated single-cell
object (same cells, same `final_celltype` labels, same `sample_id`/age structure), so results
can be cross-validated gene-by-gene and cell-type-by-cell-type across methods:

| # | Analysis | Language / tool | What it answers |
|---|----------|------------------|------------------|
| 1 | Gene co-expression networks | R -- WGCNA / hdWGCNA | Which genes co-vary, in which cell type, and how do those modules change with age? |
| 2 | Cell-cell communication | Python -- LIANA+ / Tensor-cell2cell | Who signals to whom, and how does that change with age? |
| 3 | Differential communication | R -- MultiNicheNet / NicheNet | Which specific ligand-receptor pairs are prioritized for old vs. young, and what upstream ligands could explain a given aging gene module? |
| 4 | Gene regulatory networks | Python -- pySCENIC | Which transcription factors and their regulons drive the observed expression changes? |

Analyses 2, 3, and 4 all reuse the outputs of analysis 1 (age-associated modules and hub genes)
as their cross-validation reference.

## Repository contents

```
.
├── WGCNA_24_hdwgcna_HVG_d_17_08.Rmd         # Step 1: preprocessing, integration, pseudobulk,
│                                             #   WGCNA/hdWGCNA network, module-age association,
│                                             #   hub genes, functional enrichment
├── multinichenet_wgcna_integration_17_08.Rmd # Step 3: MultiNicheNet differential communication
│                                             #   + cross-validation against WGCNA modules
├── WGCNA_to_LIANA_workflow_17_08.ipynb       # Step 2: LIANA+ cell-cell communication
│                                             #   (steady-state, young-vs-old, tensor-cell2cell)
├── WGCNA_to_SCENIC_workflow.ipynb            # Step 4: pySCENIC gene regulatory networks
├── report/                                   # LaTeX report (Overleaf-ready) documenting the
│                                             #   full methodology, with figure/table slots
├── input_data/                               # (not tracked) raw/intermediate data -- see below
├── outputs/                                  # (not tracked) result tables, figures, cached objects
└── README.md
```

> Notebooks/scripts are meant to be run **in order** (1 → 2 → 3 → 4): each step reads at least
> one output file produced by an earlier step (see *Execution order* below).

## Data availability

Raw data are not included in this repository. Heart single-cell data (FACS + droplet) come from
the public **Tabula Muris Senis** atlas. Prior-knowledge resources used by MultiNicheNet/NicheNet
(mouse ligand-receptor network, ligand-target matrix, weighted networks, ligand-TF matrix) and by
SCENIC (mouse cisTarget ranking databases + motif annotation table) are downloaded automatically
or manually as described in each script/notebook (see *Requirements* below for exact sources).

## Execution order and hand-off files

1. **`WGCNA_24_hdwgcna_HVG_d_17_08.Rmd`** (R)
   Reads the raw `.h5ad` files (FACS + droplet), integrates and QCs them, builds the pseudobulk
   matrix, constructs the co-expression network, and computes module-age correlation, trajectory
   shape, hub genes, and functional enrichment.
   **Produces** (consumed by later steps): `hub_genes_by_module.csv`,
   `all_genes_by_module.csv`, `module_age_correlation.csv`, `module_trajectory_shape.csv`,
   `heart_for_liana.h5ad` (raw-count export with `final_celltype`/`sample_id`/`age` metadata for
   the Python notebooks), plus cached GO/KEGG/Reactome enrichment tables.

2. **`WGCNA_to_LIANA_workflow_17_08.ipynb`** (Python)
   Reads `heart_for_liana.h5ad`. Runs LIANA+ (steady-state, young-vs-old, per-mouse tensor
   construction, Tensor-cell2cell factorization) and cross-validates the age-associated
   communication factor against the WGCNA hub-gene set from step 1.
   **Produces**: `liana_res_age_module_relevant.csv` (consumed by step 4), plus tensor/factor
   outputs and figures.

3. **`multinichenet_wgcna_integration_17_08.Rmd`** (R, same session/environment as step 1
   recommended -- it reuses several objects such as `module_trait_table`, `hub_genes_table`,
   `gene_celltype_module`, and `ordered_age_levels` directly from step 1's R session)
   Builds a `SingleCellExperiment`, defines a single oldest-vs-youngest contrast, and runs
   `multi_nichenet_analysis()` (pseudobulk DE, ligand activity, LR prioritization). Cross-validates
   prioritized LR pairs against age-associated WGCNA modules, and runs a mechanistic ligand-activity
   analysis using WGCNA hub genes as the target gene set.
   **Produces**: `multinichenet_output_age.rds`, `multinichenet_age_LR_wgcna_overlap.csv`, circos
   plots, LR dotplots, and a multi-sheet Excel "golden table" summarizing all prioritized
   interactions.

4. **`WGCNA_to_SCENIC_workflow.ipynb`** (Python)
   Reads `heart_for_liana.h5ad` (and, optionally, `liana_res_age_module_relevant.csv` from step 2
   for the final cross-method check). Runs GRNBoost2 → cisTarget → AUCell to build regulons,
   computes regulon-age correlation and cell-type specificity (RSS), and cross-validates regulons
   against WGCNA (TF/hub-gene overlap, target-gene-set enrichment) and against LIANA+ (regulon
   targets vs. age-relevant ligands/receptors).
   **Produces**: `scenic_results/` (adjacencies, regulons, AUCell matrix, RSS table, regulon-age
   correlation, regulon-vs-WGCNA enrichment table).

## Requirements

### R (steps 1 and 3)

R ≥ 4.2 is recommended. Install via CRAN/Bioconductor/GitHub as needed:

```r
# CRAN
install.packages(c(
  "tidyverse", "Seurat", "Matrix", "reticulate", "patchwork", "ggraph", "igraph",
  "uwot", "ggrepel", "ggforce", "ggpubr", "ggtext", "ggalluvial", "pheatmap",
  "viridis", "scales", "circlize", "RColorBrewer", "gridExtra", "magick",
  "visNetwork", "openxlsx", "jsonlite", "knitr", "conflicted", "devtools", "remotes"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "WGCNA", "clusterProfiler", "ReactomePA", "org.Mm.eg.db", "GOSemSim",
  "enrichplot", "SingleCellExperiment", "S4Vectors", "ComplexHeatmap"
))

# GitHub
devtools::install_github("smorabit/hdWGCNA")
remotes::install_github("saeyslab/harmony")        # or CRAN "harmony"
devtools::install_github("saeyslab/nichenetr")
remotes::install_github("saeyslab/multinichenetr")
```

`multinichenetr` additionally requires the NicheNet v2 mouse prior-model `.rds` files
(`lr_network_mouse_*.rds`, `ligand_target_matrix_nsga2r_final_mouse.rds`,
`weighted_networks_nsga2r_final_mouse.rds`, `ligand_tf_matrix_nsga2r_final_mouse.rds`), available
from the NicheNet Zenodo record -- download once and point the `.Rmd` to their local path.

### Python (steps 2 and 4)

Python ≥ 3.9, separate environments recommended for LIANA+/cell2cell and pySCENIC (dependency
conflicts are common between the two stacks).

```bash
# Step 2 -- LIANA+ / Tensor-cell2cell
pip install liana cell2cell decoupler scanpy anndata gseapy \
            numpy pandas scipy statsmodels matplotlib seaborn plotnine
# Optional, only if use_gpu = True for tensor factorization:
pip install torch tensorly

# Step 4 -- pySCENIC
pip install pyscenic arboreto ctxcore "dask==2024.2.1" "distributed==2024.2.1" \
            scanpy numpy pandas scipy statsmodels matplotlib seaborn
```

pySCENIC additionally requires, downloaded once into `resources/` (URLs are fetched automatically
by the first cells of the notebook):
- `mm_mgi_tfs.txt` (mouse transcription factor list)
- Two mouse cisTarget ranking databases (`mm10_500bp_up_100bp_down_full_tx_v10_clust...feather`,
  `mm10_10kbp_up_10kbp_down_full_tx_v10_clust...feather`)
- `motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl` (motif-to-TF annotation table)

### Hardware notes
- Tensor-cell2cell's elbow/rank analysis (step 2) and GRNBoost2 (step 4) are the two most
  compute-intensive stages; both can run on CPU but benefit substantially from more cores/RAM.
  GPU is optional for tensor factorization (`use_gpu = True`, requires `torch`/`tensorly`).
- `multi_nichenet_analysis()` (step 3) and the WGCNA network construction (step 1) are
  comparatively light but can still take from minutes to a few hours depending on dataset size.

## Report

The `report/` folder contains a LaTeX write-up of the full methodology,, a
Discussion/Limitations outline, and a Data and Code Availability section pointing back to this
repository.

## Citation


If you use this pipeline, please cite the underlying methods: WGCNA/hdWGCNA, LIANA+ /
Tensor-cell2cell, NicheNet/MultiNicheNet, and SCENIC/pySCENIC (see the References section of the
report in `report/` for full citations), in addition to the *Tabula Muris Senis* data source.
