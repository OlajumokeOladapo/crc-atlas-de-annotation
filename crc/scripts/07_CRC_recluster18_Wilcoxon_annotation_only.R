############################################################
# CRC LOCAL RECLUSTER 18
# CELL-LEVEL WILCOXON MARKER ANALYSIS — ANNOTATION ONLY
#
# Comparison:
#   local unresolved recluster 18 vs all other local
#   unresolved reclusters.
#
# This analysis is descriptive/annotation-supporting only.
# Recluser 18 had only one eligible paired biological replicate,
# so no formal replicate-aware pseudobulk DE was performed.
############################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratObject)
    library(data.table)
})

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset = "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

input_file <- file.path(
    atlas_root,
    "Unresolved_Reclustering_By_Condition",
    "CRC",
    "CRC_UNRESOLVED_CCA_RECLUSTERED_res0.5.rds"
)

output_dir <- file.path(
    atlas_root,
    "CRC_DE_ALL_RESULTS",
    "unresolved_clusters",
    "cluster_18",
    "Wilcoxon_annotation_only"
)

dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

cluster_col <- "unresolved_cluster_res0.5"
target_cluster <- "18"
padj_cutoff <- 0.05
log2fc_cutoff <- 0.25
strong_log2fc_cutoff <- 1

if (!file.exists(input_file)) stop("Input not found: ", input_file)

crc <- readRDS(input_file)

if (!cluster_col %in% colnames(crc[[]])) stop("Missing cluster column: ", cluster_col)

Idents(crc) <- cluster_col
DefaultAssay(crc) <- "RNA"

rna_layers <- Layers(crc[["RNA"]])
data_layers <- grep("^data", rna_layers, value=TRUE)

if (length(data_layers) > 1) {
    crc[["RNA"]] <- JoinLayers(crc[["RNA"]])
}

rna_layers_after <- Layers(crc[["RNA"]])

if (!any(grepl("^data$", rna_layers_after))) {
    crc <- NormalizeData(
        crc,
        assay="RNA",
        normalization.method="LogNormalize",
        scale.factor=10000,
        verbose=TRUE
    )
}

cluster_values <- as.character(crc[[cluster_col]][,1])
cluster_counts <- table(cluster_values)

if (!target_cluster %in% names(cluster_counts)) stop("Cluster 18 absent.")

markers <- FindMarkers(
    object=crc,
    ident.1=target_cluster,
    ident.2=NULL,
    assay="RNA",
    test.use="wilcox",
    only.pos=FALSE,
    logfc.threshold=0,
    min.pct=0,
    min.cells.feature=1,
    min.cells.group=1,
    verbose=TRUE
)

markers_dt <- as.data.table(markers, keep.rownames="gene")

if ("avg_log2FC" %in% names(markers_dt)) {
    setnames(markers_dt, "avg_log2FC", "log2FC")
} else if ("avg_logFC" %in% names(markers_dt)) {
    setnames(markers_dt, "avg_logFC", "log2FC")
} else {
    stop("Could not identify Seurat fold-change column.")
}

markers_dt[, `:=`(
    cluster=target_cluster,
    comparison="recluster_18_VS_all_other_new_reclusters",
    method="Seurat_Wilcoxon_cell_level_annotation_only",
    formal_DE=FALSE
)]

setorder(markers_dt, p_val_adj, -log2FC)

up <- markers_dt[!is.na(p_val_adj) & p_val_adj<padj_cutoff & log2FC>=log2fc_cutoff]
down <- markers_dt[!is.na(p_val_adj) & p_val_adj<padj_cutoff & log2FC<=-log2fc_cutoff]
strong <- markers_dt[!is.na(p_val_adj) & p_val_adj<padj_cutoff & log2FC>=strong_log2fc_cutoff]

setorder(up, p_val_adj, -log2FC)
setorder(down, p_val_adj, log2FC)
setorder(strong, p_val_adj, -log2FC)

fwrite(markers_dt, file.path(output_dir,"cluster18_Wilcoxon_all_genes.csv"))
fwrite(up, file.path(output_dir,"cluster18_Wilcoxon_upregulated.csv"))
fwrite(down, file.path(output_dir,"cluster18_Wilcoxon_downregulated.csv"))
fwrite(strong, file.path(output_dir,"cluster18_Wilcoxon_strong_upregulated.csv"))
fwrite(head(up,100), file.path(output_dir,"cluster18_Wilcoxon_top100_upregulated.csv"))
fwrite(head(strong,100), file.path(output_dir,"cluster18_Wilcoxon_top100_strong_upregulated.csv"))

summary_table <- data.table(
    cluster=target_cluster,
    comparison="recluster_18_VS_all_other_new_reclusters",
    method="Seurat_Wilcoxon_cell_level",
    purpose="annotation_only",
    formal_DE=FALSE,
    target_cells=as.integer(cluster_counts[target_cluster]),
    rest_cells=ncol(crc)-as.integer(cluster_counts[target_cluster]),
    genes_tested=nrow(markers_dt),
    upregulated=nrow(up),
    downregulated=nrow(down),
    strong_upregulated=nrow(strong),
    padj_cutoff=padj_cutoff,
    log2fc_cutoff=log2fc_cutoff
)

fwrite(summary_table, file.path(output_dir,"cluster18_Wilcoxon_summary.csv"))
capture.output(sessionInfo(), file=file.path(output_dir,"sessionInfo.txt"))

cat("\nTOP 25 UPREGULATED MARKERS:\n")
print(head(up[,.(gene,log2FC,pct.1,pct.2,p_val_adj)],25))
cat("\nResults saved to:\n",output_dir,"\n",sep="")
