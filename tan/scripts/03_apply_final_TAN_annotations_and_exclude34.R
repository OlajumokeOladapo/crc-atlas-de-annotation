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
    "TAN_CCA_FINAL_ANNOTATED.rds"
)

map_file <- file.path(
    atlas_root,
    "tan",
    "config",
    "TAN_final_CCA_annotation_map.csv"
)

tan <- readRDS(input_file)
map <- fread(map_file,colClasses="character")

cluster <- as.character(tan$cca_snn_res.0.5)

tan$CCA_annotation_refined <- map$final_annotation[
    match(cluster,map$cluster)
]
tan$annotation_source <- map$annotation_source[
    match(cluster,map$cluster)
]
tan$analysis_status <- map$analysis_status[
    match(cluster,map$cluster)
]
tan$exclusion_reason <- NA_character_

tan$exclusion_reason[
    tan$analysis_status=="Excluded"
] <- paste0(
    "CCA cluster 34 restricted to CELLxGENE_Moorman2025; ",
    "insufficient evidence for confident biological annotation"
)

if (
    anyNA(tan$CCA_annotation_refined) ||
    anyNA(tan$analysis_status)
) stop("Some TAN clusters are missing from the final annotation map.")

out <- file.path(atlas_root,"final_annotated_conditions")
dir.create(out,recursive=TRUE,showWarnings=FALSE)

fwrite(
    data.table(
        cell=colnames(tan)[tan$analysis_status=="Excluded"],
        cluster=cluster[tan$analysis_status=="Excluded"],
        annotation=tan$CCA_annotation_refined[tan$analysis_status=="Excluded"],
        exclusion_reason=tan$exclusion_reason[tan$analysis_status=="Excluded"]
    ),
    file.path(out,"TAN_cluster34_excluded_cells.csv")
)

tan_final <- subset(
    tan,
    cells=colnames(tan)[tan$analysis_status=="Included"]
)

saveRDS(
    tan_final,
    file.path(out,"annotated_TAN.rds"),
    compress=FALSE
)

fwrite(
    unique(data.table(
        cluster=cluster,
        annotation=tan$CCA_annotation_refined,
        annotation_source=tan$annotation_source,
        analysis_status=tan$analysis_status
    )),
    file.path(out,"TAN_annotation_summary.csv")
)
