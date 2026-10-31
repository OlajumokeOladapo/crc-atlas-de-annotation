suppressPackageStartupMessages({
    library(Seurat)
    library(data.table)
})

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset="/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

input_file <- file.path(
    atlas_root,
    "Final_Annotated_Seurat_Objects",
    "Normal_CCA_FINAL_ANNOTATED.rds"
)

map_file <- file.path(
    atlas_root,
    "normal",
    "config",
    "Normal_final_CCA_annotation_map.csv"
)

normal <- readRDS(input_file)
map <- fread(map_file,colClasses="character")

cluster <- as.character(normal$cca_snn_res.0.5)

normal$CCA_annotation_refined <- map$final_annotation[
    match(cluster,map$cluster)
]
normal$annotation_source <- map$annotation_source[
    match(cluster,map$cluster)
]

if (anyNA(normal$CCA_annotation_refined)) {
    stop("Some Normal clusters are missing from the final annotation map.")
}

out <- file.path(atlas_root,"final_annotated_conditions")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

saveRDS(
    normal,
    file.path(out,"annotated_normal.rds"),
    compress=FALSE
)

fwrite(
    unique(data.table(
        cluster=cluster,
        annotation=normal$CCA_annotation_refined,
        annotation_source=normal$annotation_source
    )),
    file.path(out,"Normal_annotation_summary.csv")
)
