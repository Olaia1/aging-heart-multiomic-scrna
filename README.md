# Aging Mouse Heart -- Single-Cell Analysis

Co-expression, cell-cell communication, and gene regulatory network analysis of aging in the
mouse heart, using single-cell RNA-seq data from *Tabula Muris Senis* (FACS + droplet).

The pipeline combines four complementary analyses run on the **same** annotated single-cell
object (same cells, same `final_celltype` labels, same `sample_id`/age structure), so results
can be cross-validated gene-by-gene and cell-type-by-cell-type across methods:

| # | Folder | Language / tool | What it answers |
|---|--------|------------------|------------------|
| 1 | `hdWGCNA/` | R -- WGCNA / hdWGCNA | Which genes co-vary, in which cell type, and how do those modules change with age? |
| 2 | `LIANA+&tensor2cell/` | Python -- LIANA+ / Tensor-cell2cell | Who signals to whom, and how does that change with age? |
| 3 | `MultiNicheNet/` | R -- MultiNicheNet / NicheNet | Which specific ligand-receptor pairs are prioritized for old vs. young, and what upstream ligands could explain a given aging gene module? |
| 4 | `SCENIC/` | Python -- pySCENIC | Which transcription factors and their regulons drive the observed expression changes? |

Folders 2, 3, and 4 all reuse the outputs of folder 1 (age-associated modules and hub genes) as
their cross-validation reference.

## Repository contents

```
.
├── hdWGCNA/
│   ├── WGCNA_24_hdwgcna_HVG_d_17_08.Rmd        # preprocessing, integration, pseudobulk
│   │                                            #   construction, co-expression network,
│   │                                            #   module-age association, hub genes, enrichment
│   ├── export_to_LIANA.R                        # exports the shared raw-count .h5ad
│   │                                            #   (final_celltype/sample_id/age) consumed by
│   │                                            #   the LIANA+ and SCENIC notebooks
│   └── outputs/
│       ├── hub_genes_by_module.csv
│       ├── all_genes_by_module.csv
│       ├── module_age_correlation.csv
│       ├── module_trajectory_shape.csv
│       ├── GO_BP_by_module_and_celltype...csv
│       ├── KEGG_by_module_and_celltype...csv
│       ├── Reactome_by_module_and_celltype...csv
│       ├── module_jaccard_overlap_heatmap...
│       ├── figure_module_age_regression...
│       ├── celular_composition_per_module...
│       ├── alluvial_module_celltype_composition...
│       ├── all_module_trajectories_combined...
│       ├── GO_BP_master_heatmap_all...
│       ├── Reactome_master_heatmap...
│       └── trajectory&pathway_complete.../         # per-module multi-panel publication figures
│
├── LIANA+&tensor2cell/
│   ├── WGCNA_to_LIANA_workflow.ipynb             # LR inference, tensor decomposition,
│   │                                             #   WGCNA cross-validation
│   └── output/
│       ├── liana_res_annotated_with_wgcna...csv   # steady-state/young-vs-old LR results,
│       │                                          #   annotated with WGCNA module membership
│       ├── liana_res_age_module_relevant...csv    # LR pairs overlapping age-associated modules
│       ├── factor_loadings_Contexts.csv           # tensor factor loadings (per mouse)
│       ├── factor_loadings_Ligand-Receptor...csv
│       ├── factor_loadings_Sender_Cells...csv
│       ├── factor_loadings_Receiver_Cells...csv
│       ├── ccc_network_edges_Factor_1.csv
│       ├── top_lr_pairs_Factor_1.csv
│       ├── top_lr_Factor_1_wgcna_annotated...csv
│       ├── sender_receiver_loadings_Factor_1...csv
│       ├── progeny_ulm_pvals_all_factors...csv     # PROGENy pathway scoring of factor loadings
│       ├── GSEA/                                   # pathway enrichment of factor loadings
│       │   ├── GSEA-Adj-Pvals.csv
│       │   ├── GSEA-Dotplot.pdf
│       │   ├── enriched_pathways_significant...csv
│       │   ├── depleted_pathways_significant...csv
│       │   └── top5_pathways_per_factor...csv
│       └── plots/
│           ├── 01_tensor_factors_plot.pdf
│           ├── 02_context_boxplot_clean.pdf
│           ├── 03_ccc_network_Factor_1.pdf
│           ├── 04_top_lr_pairs_Factor_1_heatmap...
│           ├── 05_sender_receiver_loadings...
│           ├── 06_progeny_barplot_Factor_1...
│           └── 07_progeny_heatmap_all_factors...
│
├── MultiNicheNet/
│   ├── multinichenet_wgcna_integration_17_08.Rmd  # differential communication, WGCNA
│   │                                              #   cross-validation, hub-gene ligand activity
│   ├── multinichenet_wgcna_integration...html      # rendered notebook
│   └── output/
│       ├── golden_table_LR_pairs.csv               # master table of prioritized LR interactions
│       ├── all_candidates_summary.csv
│       ├── module_celltype_candidates.csv
│       ├── multinichenet_age_LR_wgcna_overlap...csv
│       ├── GO_LR_wgcna_overlap_from_prioritized...csv
│       ├── ligand_activity_robustness_top...csv
│       ├── ligand_activity_hubgenes_module.../      # per-module ligand-activity results
│       │                                           #   (WGCNA hub genes as target set)
│       └── plots/
│
├── SCENIC/
│   ├── WGCNA_to_SCENIC_workflow.ipynb             # GRN inference, regulon-age correlation,
│   │                                              #   WGCNA/LIANA+ cross-validation
│   └── output/                                     # adjacencies (GRNBoost2), pruned regulons
│                                                   #   (cisTarget), AUCell activity matrix, RSS
│                                                   #   table, regulon-age correlation, and
│                                                   #   regulon vs. WGCNA/LIANA+ cross-validation
│                                                   #   tables
│
├── report/                                         # LaTeX methodology report (Overleaf-ready)
└── README.md
```

> Some file names above are truncated as shown in the repository file browser. Notebooks/scripts
> are meant to be run **in order** (`hdWGCNA/` → `LIANA+&tensor2cell/` → `MultiNicheNet/` →
> `SCENIC/`): each step reads at least one output file produced by an earlier step (see
> *Execution order* below).

## Data availability

Raw data are not included in this repository. Heart single-cell data (FACS + droplet) come from
the public **Tabula Muris Senis** atlas. Prior-knowledge resources used by MultiNicheNet/NicheNet
(mouse ligand-receptor network, ligand-target matrix, weighted networks, ligand-TF matrix) and by
SCENIC (mouse cisTarget ranking databases + motif annotation table) are downloaded automatically
or manually as described in each script/notebook (see *Requirements* below for exact sources).

## Execution order and hand-off files

1. **`hdWGCNA/WGCNA_24_hdwgcna_HVG_d_17_08.Rmd`** (R)
   Reads the raw `.h5ad` files (FACS + droplet), integrates and QCs them, builds the pseudobulk
   matrix, constructs the co-expression network, and computes module-age correlation, trajectory
   shape, hub genes, and functional enrichment.
   **Produces**: everything under `hdWGCNA/outputs/` (`hub_genes_by_module.csv`,
   `all_genes_by_module.csv`, `module_age_correlation.csv`, `module_trajectory_shape.csv`, cached
   GO/KEGG/Reactome enrichment tables, and the per-module trajectory/pathway figures).

2. **`hdWGCNA/export_to_LIANA.R`** (R)
   Exports the shared raw-count object (`final_celltype`/`sample_id`/`age` metadata) as a single
   `.h5ad` file, consumed by both Python notebooks below.

3. **`LIANA+&tensor2cell/WGCNA_to_LIANA_workflow.ipynb`** (Python)
   Runs LIANA+ (steady-state, young-vs-old, per-mouse tensor construction, Tensor-cell2cell
   factorization) and cross-validates the age-associated communication factor against the WGCNA
   hub-gene set from step 1.
   **Produces**: `LIANA+&tensor2cell/output/` -- annotated/filtered LR result tables, tensor
   factor-loading tables, the `GSEA/` pathway-enrichment outputs, and the `plots/` figures
   (`liana_res_age_module_relevant...csv` is consumed by step 5).

4. **`MultiNicheNet/multinichenet_wgcna_integration_17_08.Rmd`** (R, same session/environment as
   step 1 recommended -- it reuses several objects such as `module_trait_table`,
   `hub_genes_table`, `gene_celltype_module`, and `ordered_age_levels` directly from step 1's R
   session)
   Builds a `SingleCellExperiment`, defines a single oldest-vs-youngest contrast, and runs
   `multi_nichenet_analysis()` (pseudobulk DE, ligand activity, LR prioritization). Cross-validates
   prioritized LR pairs against age-associated WGCNA modules, and runs a mechanistic ligand-activity
   analysis using WGCNA hub genes as the target gene set.
   **Produces**: `MultiNicheNet/output/` -- `golden_table_LR_pairs.csv` (master prioritization
   table), the WGCNA-overlap tables, per-module ligand-activity results, and circos/dotplot
   figures.

5. **`SCENIC/WGCNA_to_SCENIC_workflow.ipynb`** (Python)
   Reads the `.h5ad` from step 2 (and, optionally, `liana_res_age_module_relevant...csv` from
   step 3 for the final cross-method check). Runs GRNBoost2 → cisTarget → AUCell to build
   regulons, computes regulon-age correlation and cell-type specificity (RSS), and cross-validates
   regulons against WGCNA (TF/hub-gene overlap, target-gene-set enrichment) and against LIANA+
   (regulon targets vs. age-relevant ligands/receptors).
   **Produces**: `SCENIC/output/` -- adjacencies, regulons, AUCell matrix, RSS table, regulon-age
   correlation, and regulon-vs-WGCNA/LIANA+ cross-validation tables.

## Requirements

### R (steps 1 and 4)

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

### Python (steps 3 and 5)

Python ≥ 3.9, separate environments recommended for LIANA+/cell2cell and pySCENIC (dependency
conflicts are common between the two stacks).

```bash
# Step 3 -- LIANA+ / Tensor-cell2cell
pip install liana cell2cell decoupler scanpy anndata gseapy \
            numpy pandas scipy statsmodels matplotlib seaborn plotnine
# Optional, only if use_gpu = True for tensor factorization:
pip install torch tensorly

# Step 5 -- pySCENIC
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
- Tensor-cell2cell's elbow/rank analysis (step 3) and GRNBoost2 (step 5) are the two most
  compute-intensive stages; both can run on CPU but benefit substantially from more cores/RAM.
  GPU is optional for tensor factorization (`use_gpu = True`, requires `torch`/`tensorly`).
- `multi_nichenet_analysis()` (step 4) and the WGCNA network construction (step 1) are
  comparatively light but can still take from minutes to a few hours depending on dataset size.

## Report

The `report/` folder contains a LaTeX write-up of the full methodology. 
