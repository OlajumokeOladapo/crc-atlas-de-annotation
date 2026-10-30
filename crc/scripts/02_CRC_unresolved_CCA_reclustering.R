############################################################
# CRC UNRESOLVED-CELL CCA RECLUSTERING
#
# Starting object:
#   Final_Annotated_Seurat_Objects/
#   CRC_CCA_FINAL_ANNOTATED.rds
#
# Original unresolved CCA clusters:
#   3, 5, 8, 9, 10, 15, 16, 18, 19, 20, 25, 26, 27, 28
#
# Important:
# - These are ORIGINAL GLOBAL CCA cluster IDs.
# - New unresolved_cluster_res0.5 IDs are LOCAL recluster IDs.
# - Existing full-object HVGs are reused.
# - FindVariableFeatures/VST is NOT rerun on the subset because
#   split Seurat v5 RNA layers previously triggered loess failures.
# - PCA uses 20 PCs.
# - CCA integration is performed with IntegrateLayers().
############################################################


suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
})


# ============================================================
# 1. PATHS
# ============================================================

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset = "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

input_file <- file.path(
    atlas_root,
    "Final_Annotated_Seurat_Objects",
    "CRC_CCA_FINAL_ANNOTATED.rds"
)

output_dir <- file.path(
    atlas_root,
    "Unresolved_Reclustering_By_Condition",
    "CRC"
)

dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

output_rds <- file.path(
    output_dir,
    "CRC_UNRESOLVED_CCA_RECLUSTERED_res0.5.rds"
)


# ============================================================
# 2. SETTINGS
# ============================================================

original_cluster_col <- "cca_snn_res.0.5"

unresolved_original_clusters <- c(
    "3",
    "5",
    "8",
    "9",
    "10",
    "15",
    "16",
    "18",
    "19",
    "20",
    "25",
    "26",
    "27",
    "28"
)

npcs <- 20
dims_use <- 1:20
resolution_use <- 0.5


# ============================================================
# 3. LOAD FULL CRC OBJECT
# ============================================================

stopifnot(
    file.exists(input_file)
)

crc <- readRDS(
    input_file
)

DefaultAssay(crc) <- "RNA"

if (!original_cluster_col %in% colnames(crc[[]])) {
    stop(
        "Missing original CRC cluster column: ",
        original_cluster_col
    )
}


# ============================================================
# 4. CHECK ORIGINAL UNRESOLVED CLUSTER COUNTS
# ============================================================

original_ids <- as.character(
    crc[[original_cluster_col]][, 1]
)

original_unresolved_counts <- table(
    original_ids[
        original_ids %in%
            unresolved_original_clusters
    ]
)

cat(
    "\nOriginal unresolved CRC CCA clusters:\n"
)

print(
    original_unresolved_counts
)

cat(
    "\nTotal unresolved cells: ",
    sum(original_unresolved_counts),
    "\n",
    sep = ""
)


# ============================================================
# 5. REUSE EXISTING FULL-ATLAS HVGs
# ============================================================

recluster_features <- VariableFeatures(
    crc[["RNA"]]
)

cat(
    "\nExisting RNA HVGs: ",
    length(recluster_features),
    "\n",
    sep = ""
)

# Fallback: if the stored RNA HVG list is unexpectedly small,
# recover genes represented in the existing PCA loadings.
if (length(recluster_features) < 500) {

    pca_candidates <- c(
        "pca",
        "pca.cca"
    )

    pca_hit <- pca_candidates[
        pca_candidates %in%
            Reductions(crc)
    ]

    if (!length(pca_hit)) {
        stop(
            "Fewer than 500 RNA HVGs and no usable PCA reduction found."
        )
    }

    loadings <- Loadings(
        crc[[pca_hit[1]]]
    )

    recluster_features <- rownames(
        loadings
    )

    cat(
        "Recovered features from ",
        pca_hit[1],
        " loadings: ",
        length(recluster_features),
        "\n",
        sep = ""
    )
}

recluster_features <- intersect(
    recluster_features,
    rownames(crc[["RNA"]])
)

if (length(recluster_features) < 500) {
    stop(
        "Too few reclustering features after matching to RNA assay."
    )
}


# ============================================================
# 6. SUBSET ORIGINAL UNRESOLVED CELLS
# ============================================================

keep_cells <- colnames(crc)[
    original_ids %in%
        unresolved_original_clusters
]

crc_unresolved <- subset(
    crc,
    cells = keep_cells
)

rm(crc)
invisible(gc())

DefaultAssay(
    crc_unresolved
) <- "RNA"

cat(
    "\nSubset cells: ",
    ncol(crc_unresolved),
    "\n",
    sep = ""
)

if (ncol(crc_unresolved) != 270276) {
    warning(
        "Expected 270,276 unresolved CRC cells; found ",
        ncol(crc_unresolved),
        ". Verify that the input object/version is the same one used in the analysis."
    )
}


# ============================================================
# 7. SCALE EXISTING HVGs
# ============================================================

# Do NOT rerun FindVariableFeatures() here.
crc_unresolved <- ScaleData(
    crc_unresolved,
    features = recluster_features,
    verbose = TRUE
)


# ============================================================
# 8. PCA
# ============================================================

crc_unresolved <- RunPCA(
    crc_unresolved,
    assay = "RNA",
    features = recluster_features,
    npcs = npcs,
    reduction.name = "pca.unresolved",
    verbose = TRUE
)


# ============================================================
# 9. CCA INTEGRATION OF UNRESOLVED SUBSET
# ============================================================

crc_unresolved <- IntegrateLayers(
    object = crc_unresolved,
    method = CCAIntegration,
    orig.reduction = "pca.unresolved",
    new.reduction = "integrated.cca.unresolved",
    assay = "RNA",
    features = recluster_features,
    normalization.method = "LogNormalize",
    dims = dims_use,
    verbose = TRUE
)


# ============================================================
# 10. UMAP / NEIGHBORS / CLUSTERING
# ============================================================

crc_unresolved <- RunUMAP(
    crc_unresolved,
    reduction = "integrated.cca.unresolved",
    dims = dims_use,
    reduction.name = "umap.cca.unresolved",
    verbose = TRUE
)

crc_unresolved <- FindNeighbors(
    crc_unresolved,
    reduction = "integrated.cca.unresolved",
    dims = dims_use,
    graph.name = c(
        "CRC_unresolved_nn",
        "CRC_unresolved_snn"
    ),
    verbose = TRUE
)

crc_unresolved <- FindClusters(
    crc_unresolved,
    graph.name = "CRC_unresolved_snn",
    resolution = resolution_use,
    cluster.name = "unresolved_cluster_res0.5",
    verbose = TRUE
)


# ============================================================
# 11. VERIFY NEW LOCAL RECLUSTER COUNTS
# ============================================================

recluster_counts <- sort(
    table(
        as.character(
            crc_unresolved$unresolved_cluster_res0.5
        )
    ),
    decreasing = FALSE
)

cat(
    "\nNew local unresolved recluster counts:\n"
)

print(
    recluster_counts
)

cat(
    "\nTotal cells after reclustering: ",
    sum(recluster_counts),
    "\n",
    sep = ""
)


# ============================================================
# 12. SAVE COUNTS AND OBJECT
# ============================================================

write.csv(
    data.frame(
        unresolved_cluster_res0.5 =
            names(recluster_counts),
        n_cells =
            as.integer(recluster_counts)
    ),
    file.path(
        output_dir,
        "CRC_unresolved_recluster_counts_res0.5.csv"
    ),
    row.names = FALSE
)

saveRDS(
    crc_unresolved,
    output_rds,
    compress = FALSE
)

capture.output(
    sessionInfo(),
    file = file.path(
        output_dir,
        "sessionInfo_reclustering.txt"
    )
)

cat(
    "\nSaved:\n",
    output_rds,
    "\n",
    sep = ""
)
