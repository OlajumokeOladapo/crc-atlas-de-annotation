#!/usr/bin/env Rscript

# ==============================================================================
# 03. Reclustering of original unresolved Adenoma/Polyp CCA populations
# ==============================================================================
#
# This script starts from the ORIGINAL annotated full Adenoma CCA object and
# extracts only the original unresolved populations. It then reclusters those
# cells using the same CCA integration framework.
#
# IMPORTANT:
#   - This is a downstream refinement step.
#   - It does NOT alter the original cca_snn_res.0.5 assignments.
#   - It does NOT rerun FindVariableFeatures() because the integrated atlas
#     construction uses retained normalized data/HVGs and the count layers in
#     the integration object are not the raw-count source for DE.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(data.table)
})

atlas_root <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"

full_file <- file.path(
  atlas_root,
  "Final_Annotated_Seurat_Objects",
  "Adenoma_polyp_CCA_FINAL_ANNOTATED.rds"
)

out_dir <- file.path(
  atlas_root,
  "Unresolved_Reclustering_By_Condition",
  "Adenoma_polyp"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

out_file <- file.path(
  out_dir,
  "Adenoma_polyp_UNRESOLVED_CCA_RECLUSTERED.rds"
)

cluster_col <- "cca_snn_res.0.5"

# These are the original CCA clusters that remained unresolved after the
# first-pass DESeq2 annotation.
old_unresolved_clusters <- c(
  "0", "1", "2", "4", "5", "7",
  "15", "17", "18", "23", "24"
)

# Reclustering settings used in the final workflow.
dims_use <- 1:20
resolution_use <- 0.5

# ==============================================================================
# 1. LOAD FULL ANNOTATED ADENOMA OBJECT
# ==============================================================================

stopifnot(
  file.exists(full_file)
)

adenoma_full <- readRDS(
  full_file
)

stopifnot(
  cluster_col %in%
    colnames(
      adenoma_full[[]]
    )
)

cat(
  "Full Adenoma cells:",
  ncol(adenoma_full),
  "\n"
)

# ==============================================================================
# 2. EXTRACT ORIGINAL UNRESOLVED CELLS
# ==============================================================================

old_cluster_ids <- as.character(
  adenoma_full[
    [
      cluster_col
    ]
  ][
    ,
    1
  ]
)

unresolved_cells <- colnames(
  adenoma_full
)[
  old_cluster_ids %in%
    old_unresolved_clusters
]

cat(
  "Original unresolved cells selected:",
  length(unresolved_cells),
  "\n"
)

unresolved <- subset(
  adenoma_full,
  cells = unresolved_cells
)

stopifnot(
  ncol(unresolved) ==
    length(unresolved_cells)
)

# Keep the original CCA cluster for provenance.
unresolved$original_CCA_cluster <- as.character(
  unresolved[
    [
      cluster_col
    ]
  ][
    ,
    1
  ]
)

# ==============================================================================
# 3. USE EXISTING NORMALIZED EXPRESSION AND RETAINED HVGs
# ==============================================================================

DefaultAssay(unresolved) <- "RNA"

recluster_features <- VariableFeatures(
  adenoma_full
)

if (
  length(recluster_features) < 100
) {
  stop(
    "Too few retained variable features in the full Adenoma object."
  )
}

recluster_features <- intersect(
  recluster_features,
  rownames(unresolved)
)

cat(
  "Features used for unresolved-cell reclustering:",
  length(recluster_features),
  "\n"
)

# CCA integration expects split normalized-expression layers. The original
# Adenoma CCA object generally already contains data.* layers. Only split if
# normalized layers are not already split.
rna_layers <- Layers(
  unresolved[["RNA"]]
)

cat(
  "RNA layers before reclustering:\n"
)

print(
  rna_layers
)

has_split_data <- any(
  grepl(
    "^data\\.",
    rna_layers
  )
)

if (!has_split_data) {

  if (
    !"study_final" %in%
      colnames(
        unresolved[[]]
      )
  ) {
    stop(
      "RNA data layers are not split and study_final metadata is unavailable."
    )
  }

  unresolved[["RNA"]] <- split(
    unresolved[["RNA"]],
    f = unresolved$study_final
  )
}

# ==============================================================================
# 4. SCALE + PCA
# ==============================================================================

unresolved <- ScaleData(
  unresolved,
  features = recluster_features,
  verbose = FALSE
)

unresolved <- RunPCA(
  unresolved,
  features = recluster_features,
  npcs = max(dims_use),
  reduction.name = "pca",
  verbose = FALSE
)

# ==============================================================================
# 5. CCA INTEGRATION
# ==============================================================================

unresolved <- IntegrateLayers(
  object = unresolved,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  features = recluster_features,
  normalization.method = "LogNormalize",
  dims = dims_use,
  verbose = FALSE
)

# ==============================================================================
# 6. GRAPH CLUSTERING
# ==============================================================================

unresolved <- FindNeighbors(
  unresolved,
  reduction = "integrated.cca",
  dims = dims_use,
  verbose = FALSE
)

unresolved <- FindClusters(
  unresolved,
  resolution = resolution_use,
  verbose = FALSE
)

unresolved$unresolved_cluster_res0.5 <- as.character(
  Idents(unresolved)
)

cat(
  "\nNew unresolved-cell clusters:\n"
)

print(
  sort(
    table(
      unresolved$unresolved_cluster_res0.5
    )
  )
)

# ==============================================================================
# 7. UMAP
# ==============================================================================

unresolved <- RunUMAP(
  unresolved,
  reduction = "integrated.cca",
  dims = dims_use,
  reduction.name = "umap.cca.unresolved",
  reduction.key = "CCAUNRES_",
  verbose = FALSE
)

# ==============================================================================
# 8. SAVE
# ==============================================================================

saveRDS(
  unresolved,
  out_file,
  compress = FALSE
)

fwrite(
  data.table(
    cell = colnames(unresolved),
    original_CCA_cluster =
      unresolved$original_CCA_cluster,
    unresolved_cluster_res0.5 =
      unresolved$unresolved_cluster_res0.5
  ),
  file.path(
    out_dir,
    "Adenoma_unresolved_original_to_recluster_mapping.csv"
  )
)

cat(
  "\nSaved unresolved reclustered object to:\n",
  out_file,
  "\n"
)
