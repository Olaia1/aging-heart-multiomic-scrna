# --- Export for LIANA+ (Python) 
seurat_combined_filt$final_celltype <- meta$final_celltype

liana_cells <- rownames(seurat_combined_filt@meta.data)[
  !is.na(seurat_combined_filt@meta.data$final_celltype)
]
seurat_liana <- subset(seurat_combined_filt, cells = liana_cells)

counts_liana  <- Matrix::t(GetAssayData(seurat_liana, layer = "counts"))
lognorm_liana <- Matrix::t(GetAssayData(seurat_liana, layer = "data"))

obs_liana <- seurat_liana@meta.data[, c("final_celltype", "sample_id", "age", "tech")]
obs_liana$final_celltype <- as.character(obs_liana$final_celltype)
obs_liana$sample_id      <- as.character(obs_liana$sample_id)
obs_liana$age            <- as.character(obs_liana$age)

var_liana <- data.frame(row.names = colnames(counts_liana))

adata_liana <- ad$AnnData(
  X   = reticulate::r_to_py(as(lognorm_liana, "CsparseMatrix")),  # log-normalized, ready for LIANA
  obs = obs_liana,
  var = var_liana
)
adata_liana$layers["counts"] <- reticulate::r_to_py(as(counts_liana, "CsparseMatrix"))

# UMAP embedding, for exploratory plots on the Python side
if ("umap" %in% Reductions(seurat_liana)) {
  umap_mat <- Embeddings(seurat_liana, "umap")[colnames(seurat_liana), , drop = FALSE]
  adata_liana$obsm[["X_umap"]] <- reticulate::r_to_py(umap_mat)
}

ad$AnnData$write_h5ad(adata_liana, "heart_for_liana.h5ad")


write.csv(gene_celltype_module, "all_genes_by_module.csv", row.names = FALSE)
