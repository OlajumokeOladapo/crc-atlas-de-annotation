#!/usr/bin/env Rscript

# ==============================================================================
# 05. Integrate refined unresolved-cell annotations into the full Adenoma atlas
#     and assign final cluster numbers 1-29
# ==============================================================================
#
# FINAL NUMBERING RULE:
#   - Original resolved CCA populations are renumbered first: 1-14.
#   - Original unresolved CCA cluster numbers disappear from the final numbering.
#   - Reclusters 0-14 are then assigned final numbers 15-29.
#
# This is metadata integration by CELL BARCODE. No clustering, CCA integration,
# or UMAP is recomputed here.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
})

atlas_root <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"

full_file <- file.path(
  atlas_root,
  "Final_Annotated_Seurat_Objects",
  "Adenoma_polyp_CCA_FINAL_ANNOTATED.rds"
)

unresolved_file <- file.path(
  atlas_root,
  "Unresolved_Reclustering_By_Condition",
  "Adenoma_polyp",
  "Adenoma_polyp_UNRESOLVED_RECLUSTERED_FINAL_ANNOTATED.rds"
)

output_dir <- file.path(
  atlas_root,
  "Final_Annotated_Seurat_Objects"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

final_object_file <- file.path(
  output_dir,
  "Adenoma_polyp_CCA_FINAL_REFINED_RENUMBERED.rds"
)

# ==============================================================================
# 1. LOAD OBJECTS
# ==============================================================================

adenoma_full <- readRDS(
  full_file
)

unresolved <- readRDS(
  unresolved_file
)

required_full <- c(
  "cca_snn_res.0.5",
  "CCA_annotation_final"
)

required_unresolved <- c(
  "unresolved_cluster_res0.5",
  "Adenoma_recluster_annotation"
)

stopifnot(
  all(
    required_full %in%
      colnames(
        adenoma_full[[]]
      )
  )
)

stopifnot(
  all(
    required_unresolved %in%
      colnames(
        unresolved[[]]
      )
  )
)

cat(
  "Full Adenoma cells:",
  ncol(adenoma_full),
  "\n"
)

cat(
  "Refined unresolved cells:",
  ncol(unresolved),
  "\n"
)

# ==============================================================================
# 2. CELL-BARCODE MATCH
# ==============================================================================

overlap_cells <- intersect(
  colnames(adenoma_full),
  colnames(unresolved)
)

cat(
  "Refined unresolved cells found in full object:",
  length(overlap_cells),
  "/",
  ncol(unresolved),
  "\n"
)

if (
  length(overlap_cells) !=
    ncol(unresolved)
) {
  stop(
    "Not every refined unresolved cell is present in the full Adenoma object."
  )
}

# ==============================================================================
# 3. KEEP ORIGINAL ANNOTATIONS, REPLACE ONLY RECLUSTERED CELLS
# ==============================================================================

adenoma_full$CCA_annotation_refined <- as.character(
  adenoma_full$CCA_annotation_final
)

annotation_lookup <- setNames(
  as.character(
    unresolved$Adenoma_recluster_annotation
  ),
  colnames(unresolved)
)

recluster_lookup <- setNames(
  as.character(
    unresolved$unresolved_cluster_res0.5
  ),
  colnames(unresolved)
)

adenoma_full$CCA_annotation_refined[
  overlap_cells
] <- unname(
  annotation_lookup[
    overlap_cells
  ]
)

adenoma_full$unresolved_recluster_id <- NA_character_

adenoma_full$unresolved_recluster_id[
  overlap_cells
] <- unname(
  recluster_lookup[
    overlap_cells
  ]
)

adenoma_full$annotation_source <- "Original CCA annotation"

adenoma_full$annotation_source[
  overlap_cells
] <- "Unresolved-cell reclustering"

# ==============================================================================
# 4. DEFINE ORIGINAL RESOLVED VS UNRESOLVED CCA CLUSTERS
# ==============================================================================

old_resolved_clusters <- c(
  "3", "6", "8", "9", "10", "11", "12",
  "13", "14", "16", "19", "20", "21", "22"
)

old_unresolved_clusters <- c(
  "0", "1", "2", "4", "5", "7",
  "15", "17", "18", "23", "24"
)

original_cluster_ids <- as.character(
  adenoma_full$cca_snn_res.0.5
)

# Safety: every reclustered cell must originate from an original unresolved CCA
# cluster.
reclustered_original_ids <- original_cluster_ids[
  match(
    overlap_cells,
    colnames(adenoma_full)
  )
]

bad_recluster_sources <- setdiff(
  unique(
    reclustered_original_ids
  ),
  old_unresolved_clusters
)

if (
  length(
    bad_recluster_sources
  ) > 0
) {
  stop(
    paste(
      "Reclustered cells unexpectedly came from original resolved cluster(s):",
      paste(
        bad_recluster_sources,
        collapse = ", "
      )
    )
  )
}

# ==============================================================================
# 5. FINAL NUMBERING MAP
# ==============================================================================

old_number_map <- data.table(
  source_type = "Original resolved CCA",
  source_cluster = old_resolved_clusters,
  final_cluster_number = seq_along(
    old_resolved_clusters
  )
)

new_number_map <- data.table(
  source_type = "Unresolved-cell recluster",
  source_cluster = as.character(
    0:14
  ),
  final_cluster_number = seq(
    from = max(
      old_number_map$final_cluster_number
    ) + 1L,
    length.out = 15L
  )
)

old_number_lookup <- setNames(
  old_number_map$final_cluster_number,
  old_number_map$source_cluster
)

new_number_lookup <- setNames(
  new_number_map$final_cluster_number,
  new_number_map$source_cluster
)

# ==============================================================================
# 6. ASSIGN FINAL CLUSTER NUMBERS
# ==============================================================================

adenoma_full$final_cluster_number <- NA_integer_

resolved_cells <- original_cluster_ids %in%
  old_resolved_clusters

adenoma_full$final_cluster_number[
  resolved_cells
] <- unname(
  old_number_lookup[
    original_cluster_ids[
      resolved_cells
    ]
  ]
)

recluster_ids <- as.character(
  adenoma_full$unresolved_recluster_id
)

reclustered_cells <- !is.na(
  recluster_ids
)

adenoma_full$final_cluster_number[
  reclustered_cells
] <- unname(
  new_number_lookup[
    recluster_ids[
      reclustered_cells
    ]
  ]
)

if (
  any(
    is.na(
      adenoma_full$final_cluster_number
    )
  )
) {
  missing_sources <- sort(
    unique(
      original_cluster_ids[
        is.na(
          adenoma_full$final_cluster_number
        )
      ]
    )
  )

  stop(
    paste(
      "Some cells have no final cluster number. Original clusters:",
      paste(
        missing_sources,
        collapse = ", "
      )
    )
  )
}

# ==============================================================================
# 7. SOURCE-ID / FINAL DISPLAY LABEL
# ==============================================================================

adenoma_full$final_cluster_source_id <- paste0(
  "CCA_",
  original_cluster_ids
)

adenoma_full$final_cluster_source_id[
  reclustered_cells
] <- paste0(
  "RECLUSTER_",
  recluster_ids[
    reclustered_cells
  ]
)

adenoma_full$final_cluster_label <- paste0(
  adenoma_full$final_cluster_number,
  ": ",
  adenoma_full$CCA_annotation_refined
)

# ==============================================================================
# 8. FINAL NUMBER/ANNOTATION TABLE
# ==============================================================================

final_summary <- data.table(
  final_cluster_number =
    adenoma_full$final_cluster_number,
  final_annotation =
    as.character(
      adenoma_full$CCA_annotation_refined
    ),
  annotation_source =
    as.character(
      adenoma_full$annotation_source
    ),
  final_cluster_source_id =
    as.character(
      adenoma_full$final_cluster_source_id
    )
)[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    final_cluster_number,
    final_annotation,
    annotation_source,
    final_cluster_source_id
  )
]

setorder(
  final_summary,
  final_cluster_number
)

cat(
  "\nFINAL ADENOMA CLUSTERS 1-29:\n"
)

print(
  final_summary,
  nrows = Inf
)

if (
  !identical(
    sort(
      unique(
        adenoma_full$final_cluster_number
      )
    ),
    1:29
  )
) {
  stop(
    "Final cluster numbering is not exactly 1:29."
  )
}

fwrite(
  final_summary,
  file.path(
    output_dir,
    "Adenoma_FINAL_cluster_number_annotation_map.csv"
  )
)

# ==============================================================================
# 9. SAVE OLD-TO-NEW CELL FLOW FOR PROVENANCE
# ==============================================================================

mapping_summary <- data.table(
  original_cluster =
    original_cluster_ids,
  original_annotation =
    as.character(
      adenoma_full$CCA_annotation_final
    ),
  recluster_id =
    as.character(
      adenoma_full$unresolved_recluster_id
    ),
  final_cluster_number =
    adenoma_full$final_cluster_number,
  final_annotation =
    as.character(
      adenoma_full$CCA_annotation_refined
    ),
  annotation_source =
    as.character(
      adenoma_full$annotation_source
    )
)[
  ,
  .(
    n_cells = .N
  ),
  by = .(
    original_cluster,
    original_annotation,
    recluster_id,
    final_cluster_number,
    final_annotation,
    annotation_source
  )
]

mapping_summary[
  ,
  original_cluster_numeric :=
    as.numeric(
      original_cluster
    )
]

mapping_summary[
  ,
  recluster_numeric :=
    as.numeric(
      recluster_id
    )
]

setorder(
  mapping_summary,
  original_cluster_numeric,
  recluster_numeric
)

mapping_summary[
  ,
  c(
    "original_cluster_numeric",
    "recluster_numeric"
  ) := NULL
]

fwrite(
  mapping_summary,
  file.path(
    output_dir,
    "Adenoma_original_CCA_to_final_1_29_mapping.csv"
  )
)

# ==============================================================================
# 10. ORDER LEGEND BY FINAL CLUSTER NUMBER
# ==============================================================================

legend_levels <- unique(
  final_summary[
    order(
      final_cluster_number
    ),
    paste0(
      final_cluster_number,
      ": ",
      final_annotation
    )
  ]
)

adenoma_full$final_cluster_label <- factor(
  adenoma_full$final_cluster_label,
  levels = legend_levels
)

# ==============================================================================
# 11. FULL FINAL UMAP — ORIGINAL FULL-ATLAS EMBEDDING
# ==============================================================================

if (
  "umap.cca" %in%
    Reductions(
      adenoma_full
    )
) {
  umap_use <- "umap.cca"
} else if (
  "umap" %in%
    Reductions(
      adenoma_full
    )
) {
  umap_use <- "umap"
} else {
  stop(
    "No full-Adenoma UMAP reduction found."
  )
}

umap_df <- as.data.table(
  Embeddings(
    adenoma_full,
    reduction = umap_use
  ),
  keep.rownames = "cell"
)

setnames(
  umap_df,
  colnames(
    umap_df
  )[2:3],
  c(
    "UMAP_1",
    "UMAP_2"
  )
)

meta_df <- data.table(
  cell = colnames(
    adenoma_full
  ),
  final_cluster_number =
    adenoma_full$final_cluster_number,
  final_annotation =
    as.character(
      adenoma_full$CCA_annotation_refined
    )
)

umap_df <- merge(
  umap_df,
  meta_df,
  by = "cell",
  sort = FALSE
)

cluster_centers <- umap_df[
  ,
  .(
    UMAP_1 = median(
      UMAP_1,
      na.rm = TRUE
    ),
    UMAP_2 = median(
      UMAP_2,
      na.rm = TRUE
    )
  ),
  by = .(
    final_cluster_number,
    final_annotation
  )
]

setorder(
  cluster_centers,
  final_cluster_number
)

p_final <- DimPlot(
  adenoma_full,
  reduction = umap_use,
  group.by = "final_cluster_label",
  raster = TRUE,
  label = FALSE
) +
  geom_text(
    data = cluster_centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = final_cluster_number
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    title = "Adenoma/Polyp",
    x = "UMAP 1",
    y = "UMAP 2",
    color = "Cluster: Annotation"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    legend.text = element_text(
      size = 8
    )
  ) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 4,
        alpha = 1
      )
    )
  )

print(
  p_final
)

ggsave(
  file.path(
    output_dir,
    "Adenoma_FULL_FINAL_RENUMBERED_UMAP.pdf"
  ),
  p_final,
  width = 15,
  height = 10
)

ggsave(
  file.path(
    output_dir,
    "Adenoma_FULL_FINAL_RENUMBERED_UMAP.png"
  ),
  p_final,
  width = 15,
  height = 10,
  dpi = 300
)

# ==============================================================================
# 12. SAVE FINAL SEURAT OBJECT
# ==============================================================================

saveRDS(
  adenoma_full,
  final_object_file,
  compress = FALSE
)

cat(
  "\nFinal Adenoma Seurat object saved to:\n",
  final_object_file,
  "\n"
)

cat(
  "\nFinal UMAP + mapping tables saved to:\n",
  output_dir,
  "\n"
)
