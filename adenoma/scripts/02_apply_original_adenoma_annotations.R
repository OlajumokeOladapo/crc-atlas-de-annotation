#!/usr/bin/env Rscript

# ==============================================================================
# Apply original Adenoma/Polyp CCA annotations after all-cluster DESeq2 analysis
# ==============================================================================
#
# This is the PRE-RECLUSTER annotation map.
# Clusters labeled "Unresolved" here were subsequently extracted for the
# separate unresolved-cell reclustering workflow.
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

atlas_root <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"

cca_file <- file.path(
  atlas_root,
  "Adenoma_Harmony_CCA_with_Chen",
  "Adenoma_CCA_integrated_WITH_CHEN.rds"
)

final_object_dir <- file.path(
  atlas_root,
  "Final_Annotated_Seurat_Objects"
)

dir.create(
  final_object_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

output_file <- file.path(
  final_object_dir,
  "Adenoma_polyp_CCA_FINAL_ANNOTATED.rds"
)

cluster_col <- "cca_snn_res.0.5"

obj <- readRDS(
  cca_file
)

stopifnot(
  cluster_col %in%
    colnames(
      obj[[]]
    )
)

# ==============================================================================
# Original annotation map BEFORE unresolved-cluster reclustering
# ==============================================================================

adenoma_annotation <- c(
  "0"  = "Unresolved",
  "1"  = "Unresolved",
  "2"  = "Unresolved",
  "3"  = "WNT-high stem/progenitor epithelial",
  "4"  = "Unresolved",
  "5"  = "Unresolved",
  "6"  = "Myofibroblast-like stroma",
  "7"  = "Unresolved",
  "8"  = "GPNMB+ / MMP12+ macrophages",
  "9"  = "Inflammatory/regenerative epithelial",
  "10" = "Vascular mural/endothelial-like cells",
  "11" = "Metabolic absorptive epithelial",
  "12" = "Activated antigen-presenting B cells",
  "13" = "Activated T / gamma-delta T–NK-like lymphocytes",
  "14" = "Enteroendocrine cells",
  "15" = "Unresolved",
  "16" = "Enteric glia / neural-like cells",
  "17" = "Ribosome/OXPHOS-high unresolved",
  "18" = "Unresolved",
  "19" = "Goblet-like secretory epithelial",
  "20" = "Mature absorptive enterocyte/colonocyte",
  "21" = "Endothelial cells",
  "22" = "CD163+ tissue-resident-like macrophages",
  "23" = "Unresolved",
  "24" = "Unresolved"
)

cluster_ids <- as.character(
  obj[
    [
      cluster_col
    ]
  ][
    ,
    1
  ]
)

annotation_values <- unname(
  adenoma_annotation[
    cluster_ids
  ]
)

if (
  any(
    is.na(
      annotation_values
    )
  )
) {

  missing_clusters <- sort(
    unique(
      cluster_ids[
        is.na(
          annotation_values
        )
      ]
    )
  )

  stop(
    paste(
      "Missing annotation(s) for cluster(s):",
      paste(
        missing_clusters,
        collapse = ", "
      )
    )
  )
}

obj$CCA_annotation_final <-
  annotation_values

annotation_table <- unique(
  data.table(
    cluster =
      cluster_ids,
    annotation =
      obj$CCA_annotation_final
  )
)

annotation_table[
  ,
  cluster_numeric :=
    as.numeric(
      cluster
    )
]

setorder(
  annotation_table,
  cluster_numeric
)

annotation_table[
  ,
  cluster_numeric := NULL
]

cat(
  "\nOriginal Adenoma CCA annotation map:\n"
)

print(
  annotation_table,
  nrows = Inf
)

cat(
  "\nCell counts by final annotation:\n"
)

print(
  sort(
    table(
      obj$CCA_annotation_final
    ),
    decreasing = TRUE
  )
)

fwrite(
  annotation_table,
  file.path(
    dirname(
      output_file
    ),
    "Adenoma_CCA_original_DESeq2_annotation_map.csv"
  )
)

saveRDS(
  obj,
  output_file,
  compress = FALSE
)

cat(
  "\nAnnotated pre-reclustering Adenoma object saved to:\n",
  output_file,
  "\n"
)
