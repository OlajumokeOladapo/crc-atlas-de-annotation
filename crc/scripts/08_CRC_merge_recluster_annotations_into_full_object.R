############################################################
# CRC: MERGE LOCAL UNRESOLVED-RECLUSTER ANNOTATIONS
# BACK INTO THE FULL CRC CCA OBJECT
############################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(data.table)
})

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset="/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

full_crc_file <- file.path(
    atlas_root,
    "Final_Annotated_Seurat_Objects",
    "CRC_CCA_FINAL_ANNOTATED.rds"
)

recluster_file <- file.path(
    atlas_root,
    "Unresolved_Reclustering_By_Condition",
    "CRC",
    "CRC_UNRESOLVED_CCA_RECLUSTERED_res0.5.rds"
)

script_dir <- tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile)),
    error=function(e) getwd()
)

map_candidates <- c(
    file.path(script_dir, "..", "config", "CRC_unresolved_recluster_annotation_map.csv"),
    file.path(getwd(), "config", "CRC_unresolved_recluster_annotation_map.csv"),
    file.path(atlas_root, "crc", "config", "CRC_unresolved_recluster_annotation_map.csv")
)

map_file <- map_candidates[file.exists(map_candidates)][1]
if (is.na(map_file)) stop("Annotation-map CSV not found.")

output_dir <- file.path(
    atlas_root,
    "Final_Annotated_Seurat_Objects"
)
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

output_rds <- file.path(
    output_dir,
    "CRC_CCA_FINAL_REFINED_RENUMBERED.rds"
)

crc <- readRDS(full_crc_file)
crc_unresolved <- readRDS(recluster_file)
map <- fread(map_file, colClasses="character")

original_cluster_col <- "cca_snn_res.0.5"
original_annotation_col <- "CCA_annotation_final"
recluster_col <- "unresolved_cluster_res0.5"

original_unresolved_clusters <- c(
    "3","5","8","9","10","15","16",
    "18","19","20","25","26","27","28"
)

required_map <- c("unresolved_recluster_id","annotation")
if (!all(required_map %in% names(map))) stop("Invalid annotation map.")
if (anyDuplicated(map$unresolved_recluster_id)) stop("Duplicate recluster IDs in map.")

local_ids <- as.character(crc_unresolved[[recluster_col]][,1])
if (length(setdiff(unique(local_ids),map$unresolved_recluster_id))) {
    stop("Annotation map is incomplete.")
}

recluster_annotation <- map$annotation[
    match(local_ids,map$unresolved_recluster_id)
]
if (anyNA(recluster_annotation)) stop("Some local recluster IDs were not mapped.")

full_original_cluster <- as.character(crc[[original_cluster_col]][,1])
full_unresolved_cells <- colnames(crc)[
    full_original_cluster %in% original_unresolved_clusters
]
recluster_cells <- colnames(crc_unresolved)

if (!setequal(full_unresolved_cells,recluster_cells)) {
    stop("Full unresolved-cell set does not match the reclustered object.")
}

crc$CCA_annotation_refined <- as.character(
    crc[[original_annotation_col]][,1]
)
crc$unresolved_recluster_id <- NA_character_
crc$unresolved_recluster_annotation <- NA_character_
crc$annotation_source <- "Original CCA annotation"

idx <- match(recluster_cells,colnames(crc))
if (anyNA(idx)) stop("Some reclustered cells do not match the full CRC object.")

crc$unresolved_recluster_id[idx] <- local_ids
crc$unresolved_recluster_annotation[idx] <- recluster_annotation
crc$CCA_annotation_refined[idx] <- recluster_annotation
crc$annotation_source[idx] <- "Unresolved-cell reclustering"

# Requested cleanup of two original CCA labels.
# Numbering is handled by the final plotting/numbering script.
crc$CCA_annotation_refined[
    crc$CCA_annotation_refined=="Dendritic cells — mixed activated subsets"
] <- "Dendritic cells"

crc$CCA_annotation_refined[
    crc$CCA_annotation_refined=="Donor-limited lineage-divergent epithelial"
] <- "Epithelial cells"

fwrite(
    data.table(
        cell=colnames(crc),
        original_CCA_cluster=as.character(crc[[original_cluster_col]][,1]),
        original_CCA_annotation=as.character(crc[[original_annotation_col]][,1]),
        unresolved_recluster_id=crc$unresolved_recluster_id,
        unresolved_recluster_annotation=crc$unresolved_recluster_annotation,
        annotation_source=crc$annotation_source,
        final_annotation=crc$CCA_annotation_refined
    ),
    file.path(output_dir,"CRC_final_annotation_cell_audit.csv")
)

saveRDS(crc,output_rds,compress=FALSE)
capture.output(sessionInfo(),file=file.path(output_dir,"CRC_final_annotation_sessionInfo.txt"))

cat("\nSaved refined CRC object:\n",output_rds,"\n",sep="")
