#!/usr/bin/env Rscript

# ==============================================================================
# 04. Differential expression and annotation of reclustered unresolved
#     Adenoma/Polyp populations
# ==============================================================================
#
# Downstream of:
#   03_recluster_original_unresolved_cells_CCA.R
#
# This stage uses:
#   - strict Adenoma raw-count pseudobulk construction,
#   - exclusion of GSE261388 inflammatory polyps P1/P12,
#   - GSE261388 P3/P4 retained,
#   - DESeq2 LRT for new reclusters when biological replication supports it,
#   - descriptive CPM + paired sample-level Wilcoxon for cluster 1,
#   - cell-level Seurat Wilcoxon for sample-specific clusters 9 and 14,
#   - final annotation of new reclusters 0-14.
#
# Cell-level Wilcoxon for clusters 9/14 is marker-discovery evidence only and
# is not treated as independent-biological-replicate inference.
#
# ==============================================================================

# ==============================================================================
# Adenoma/Polyp unresolved-cluster DE and annotation workflow
# ==============================================================================
#
# PURPOSE
#   Reproduce the final DE/annotation workflow used for the already-reclustered
#   unresolved Adenoma/Polyp cells.
#
# IMPORTANT
#   - This script DOES NOT recluster cells.
#   - Formal pseudobulk DE uses ORIGINAL raw RNA counts from study-level objects.
#   - GSE261388 P1/P12 inflammatory polyps are excluded from formal Adenoma DE.
#   - GSE261388 P3/P4 adenomatous polyps are retained.
#   - DESeq2 LRT is the primary formal DE method.
#   - Cluster 1 receives an additional paired-sample Wilcoxon rescue analysis.
#   - Clusters 9 and 14 receive cell-level Seurat Wilcoxon marker tests for
#     annotation only because they are dominated by one biological sample.
#
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(DESeq2)
  library(ggplot2)
})

# ==============================================================================
# 0. USER SETTINGS
# ==============================================================================

atlas_root <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"

adenoma_dir <- file.path(
  atlas_root,
  "Unresolved_Reclustering_By_Condition",
  "Adenoma_polyp"
)

source_dir <- file.path(
  atlas_root,
  "full_seurat_QC_normalized_by_study_condition"
)

# Existing reclustered unresolved Adenoma/Polyp object.
# No reclustering is performed in this workflow.
recluster_candidates <- c(
  file.path(
    adenoma_dir,
    "Adenoma_polyp_UNRESOLVED_RECLUSTERED_FINAL_ANNOTATED.rds"
  ),
  file.path(
    adenoma_dir,
    "Adenoma_polyp_UNRESOLVED_CCA_RECLUSTERED.rds"
  ),
  file.path(
    adenoma_dir,
    "Adenoma_polyp_UNRESOLVED_RECLUSTERED.rds"
  )
)

cluster_col <- "unresolved_cluster_res0.5"

# Primary formal DE thresholds
min_cells_target <- 10
min_cells_rest <- 10
min_count_per_sample <- 10
min_samples_with_counts <- 3

padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Cluster-1 paired Wilcoxon
cluster1_id <- "1"

# Cell-level Wilcoxon annotation-only clusters
annotation_only_wilcox_clusters <- c("9", "14")

# ==============================================================================
# 1. OUTPUT DIRECTORIES
# ==============================================================================

pb_dir <- file.path(
  adenoma_dir,
  "DESeq2_ADENOMA_ONLY_LRT"
)

cluster1_dir <- file.path(
  adenoma_dir,
  "Cluster1_log2CPM_and_paired_Wilcoxon"
)

wilcox_9_14_dir <- file.path(
  adenoma_dir,
  "Wilcoxon_cluster9_cluster14"
)

plot_dir <- file.path(
  adenoma_dir,
  "Final_annotation_plots"
)

dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cluster1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wilcox_9_14_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

find_first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

find_source_file <- function(study, source_dir) {
  files <- list.files(
    source_dir,
    pattern = paste0("^", study, ".*Adenoma_polyp.*\\.rds$"),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) == 0) {
    stop("No Adenoma_polyp source RDS found for ", study, " in ", source_dir)
  }

  preferred <- files[
    grepl("QC_normalized", basename(files), ignore.case = TRUE)
  ]

  if (length(preferred) > 0) return(preferred[1])
  files[1]
}

get_joined_raw_counts <- function(obj, assay = "RNA") {
  DefaultAssay(obj) <- assay
  assay_obj <- obj[[assay]]
  layers <- Layers(assay_obj)

  count_layers <- grep("^counts", layers, value = TRUE)

  if (length(count_layers) == 0) {
    stop("No raw-count layer found in assay ", assay)
  }

  if (length(count_layers) > 1 || !identical(count_layers, "counts")) {
    assay_obj <- JoinLayers(assay_obj)
    obj[[assay]] <- assay_obj
  }

  counts <- GetAssayData(obj, assay = assay, layer = "counts")

  vals <- counts@x
  if (length(vals) > 0) {
    if (any(vals < 0)) stop("Negative values detected in raw counts.")
    if (any(abs(vals - round(vals)) > 1e-8)) {
      stop(
        "Non-integer values detected in the raw-count layer. ",
        "Do not use normalized/integrated values for pseudobulk DE."
      )
    }
  }

  counts
}

aggregate_counts_by_group <- function(counts, cell_table, group_col = "group_id") {
  groups <- unique(cell_table[[group_col]])
  out <- vector("list", length(groups))
  names(out) <- groups

  for (g in groups) {
    cells <- cell_table[get(group_col) == g, cell]
    cells <- intersect(cells, colnames(counts))
    if (length(cells) == 0) next

    out[[g]] <- Matrix::rowSums(
      counts[, cells, drop = FALSE]
    )
  }

  out <- out[!vapply(out, is.null, logical(1))]
  mat <- do.call(cbind, out)
  rownames(mat) <- rownames(counts)

  if (is.null(dim(mat))) {
    mat <- matrix(
      mat,
      ncol = 1,
      dimnames = list(rownames(counts), names(out))
    )
  }

  mat
}

build_target_rest <- function(
    pb_counts,
    pb_meta,
    target_cluster,
    min_target_cells = 10,
    min_rest_cells = 10
) {
  target_cluster <- as.character(target_cluster)
  reps <- unique(pb_meta$replicate)

  target_list <- list()
  rest_list <- list()
  support_list <- list()

  for (rep_id in reps) {
    target_ix <- which(
      pb_meta$replicate == rep_id &
        as.character(pb_meta$cluster) == target_cluster
    )

    rest_ix <- which(
      pb_meta$replicate == rep_id &
        as.character(pb_meta$cluster) != target_cluster
    )

    if (length(target_ix) == 0 || length(rest_ix) == 0) next

    target_cells <- sum(pb_meta$n_cells[target_ix])
    rest_cells <- sum(pb_meta$n_cells[rest_ix])

    if (
      target_cells < min_target_cells ||
      rest_cells < min_rest_cells
    ) next

    target_list[[rep_id]] <- Matrix::rowSums(
      pb_counts[, target_ix, drop = FALSE]
    )

    rest_list[[rep_id]] <- Matrix::rowSums(
      pb_counts[, rest_ix, drop = FALSE]
    )

    support_list[[rep_id]] <- data.table(
      replicate = rep_id,
      target_cells = target_cells,
      rest_cells = rest_cells
    )
  }

  if (length(target_list) == 0) return(NULL)

  target_counts <- do.call(cbind, target_list)
  rest_counts <- do.call(cbind, rest_list)

  rownames(target_counts) <- rownames(pb_counts)
  rownames(rest_counts) <- rownames(pb_counts)

  support <- rbindlist(support_list)

  list(
    target_counts = target_counts,
    rest_counts = rest_counts,
    support = support,
    n_pairs = ncol(target_counts)
  )
}

run_deseq2_lrt_cluster <- function(
    pb_counts,
    pb_meta,
    target_cluster,
    output_dir,
    min_target_cells = 10,
    min_rest_cells = 10,
    min_count = 10,
    min_count_samples = 3,
    padj_cutoff = 0.05,
    log2fc_cutoff = 1
) {
  cat(
    "\n============================================================\n",
    "DESeq2 LRT: cluster ", target_cluster, " vs rest\n",
    "============================================================\n",
    sep = ""
  )

  pair_data <- build_target_rest(
    pb_counts = pb_counts,
    pb_meta = pb_meta,
    target_cluster = target_cluster,
    min_target_cells = min_target_cells,
    min_rest_cells = min_rest_cells
  )

  if (is.null(pair_data) || pair_data$n_pairs < 3) {
    cat("SKIPPED: fewer than 3 eligible paired biological samples.\n")
    return(
      data.table(
        cluster = as.character(target_cluster),
        status = "SKIPPED",
        n_pairs = ifelse(is.null(pair_data), 0, pair_data$n_pairs),
        error = "Fewer than 3 eligible paired biological samples"
      )
    )
  }

  target_counts <- pair_data$target_counts
  rest_counts <- pair_data$rest_counts
  support <- pair_data$support

  keep_gene <- (
    rowSums(target_counts >= min_count) >= min_count_samples &
      rowSums(rest_counts >= min_count) >= min_count_samples
  )

  target_counts_f <- target_counts[keep_gene, , drop = FALSE]
  rest_counts_f <- rest_counts[keep_gene, , drop = FALSE]

  cat("Eligible paired biological samples:", pair_data$n_pairs, "\n")
  cat("Genes entering DESeq2:", nrow(target_counts_f), "\n")

  if (nrow(target_counts_f) < 2) {
    return(
      data.table(
        cluster = as.character(target_cluster),
        status = "FAILED",
        n_pairs = pair_data$n_pairs,
        error = "Too few genes after filtering"
      )
    )
  }

  reps <- colnames(target_counts_f)

  count_matrix <- cbind(
    rest_counts_f,
    target_counts_f
  )

  colnames(count_matrix) <- c(
    paste0(reps, "__rest"),
    paste0(reps, "__target")
  )

  coldata <- data.frame(
    row.names = colnames(count_matrix),
    replicate = factor(c(reps, reps)),
    group = factor(
      c(rep("rest", length(reps)), rep("target", length(reps))),
      levels = c("rest", "target")
    )
  )

  full_mm <- model.matrix(~ replicate + group, data = coldata)
  reduced_mm <- model.matrix(~ replicate, data = coldata)

  cat("Full-model rank:", qr(full_mm)$rank, "/", ncol(full_mm), "\n")
  cat("Reduced-model rank:", qr(reduced_mm)$rank, "/", ncol(reduced_mm), "\n")

  if (
    qr(full_mm)$rank < ncol(full_mm) ||
    qr(reduced_mm)$rank < ncol(reduced_mm)
  ) {
    return(
      data.table(
        cluster = as.character(target_cluster),
        status = "FAILED",
        n_pairs = pair_data$n_pairs,
        error = "Model matrix is not full rank"
      )
    )
  }

  result <- tryCatch({
    dds <- DESeqDataSetFromMatrix(
      countData = round(count_matrix),
      colData = coldata,
      design = ~ replicate + group
    )

    dds <- DESeq(
      dds,
      test = "LRT",
      reduced = ~ replicate,
      fitType = "local",
      quiet = TRUE
    )

    res <- results(
      dds,
      contrast = c("group", "target", "rest")
    )

    res_dt <- as.data.table(
      as.data.frame(res),
      keep.rownames = "gene"
    )

    setnames(res_dt, old = "log2FoldChange", new = "log2FC")

    res_dt[, new_cluster := as.character(target_cluster)]
    res_dt[, n_paired_replicates := pair_data$n_pairs]

    setcolorder(
      res_dt,
      c(
        "new_cluster",
        "gene",
        "n_paired_replicates",
        setdiff(
          names(res_dt),
          c("new_cluster", "gene", "n_paired_replicates")
        )
      )
    )

    setorder(res_dt, padj, -log2FC)

    up <- res_dt[
      !is.na(padj) &
        padj < padj_cutoff &
        log2FC > log2fc_cutoff
    ]

    down <- res_dt[
      !is.na(padj) &
        padj < padj_cutoff &
        log2FC < -log2fc_cutoff
    ]

    fwrite(
      res_dt,
      file.path(
        output_dir,
        paste0(
          "Adenoma_cluster_",
          target_cluster,
          "_VS_rest_DESeq2_LRT_ALL.csv"
        )
      )
    )

    fwrite(
      up,
      file.path(
        output_dir,
        paste0(
          "Adenoma_cluster_",
          target_cluster,
          "_VS_rest_DESeq2_LRT_UP.csv"
        )
      )
    )

    fwrite(
      down,
      file.path(
        output_dir,
        paste0(
          "Adenoma_cluster_",
          target_cluster,
          "_VS_rest_DESeq2_LRT_DOWN.csv"
        )
      )
    )

    fwrite(
      head(up, 100),
      file.path(
        output_dir,
        paste0(
          "Adenoma_cluster_",
          target_cluster,
          "_VS_rest_DESeq2_LRT_TOP100_UP.csv"
        )
      )
    )

    fwrite(
      support,
      file.path(
        output_dir,
        paste0(
          "Adenoma_cluster_",
          target_cluster,
          "_sample_support.csv"
        )
      )
    )

    cat("Valid adjusted P values:", sum(!is.na(res_dt$padj)), "\n")
    cat("Significant positive markers:", nrow(up), "\n")
    cat("Significant negative markers:", nrow(down), "\n")

    data.table(
      cluster = as.character(target_cluster),
      status = "SUCCESS",
      n_pairs = pair_data$n_pairs,
      n_genes_tested = nrow(res_dt),
      n_up = nrow(up),
      n_down = nrow(down),
      error = NA_character_
    )
  }, error = function(e) {
    cat("DESeq2 FAILED:", conditionMessage(e), "\n")
    data.table(
      cluster = as.character(target_cluster),
      status = "FAILED",
      n_pairs = pair_data$n_pairs,
      n_genes_tested = nrow(target_counts_f),
      n_up = NA_integer_,
      n_down = NA_integer_,
      error = conditionMessage(e)
    )
  })

  result
}

# ==============================================================================
# 3. LOAD EXISTING RECLUSTERED OBJECT
# ==============================================================================

recluster_file <- find_first_existing(recluster_candidates)

if (is.null(recluster_file)) {
  stop(
    "Could not find the existing Adenoma/Polyp reclustered unresolved object.\n",
    "Checked:\n",
    paste(recluster_candidates, collapse = "\n")
  )
}

cat("\nLoading existing reclustered object:\n", recluster_file, "\n")

unresolved <- readRDS(recluster_file)

if (!cluster_col %in% colnames(unresolved[[]])) {
  stop("Missing cluster column: ", cluster_col)
}

cat("\nCells in existing unresolved object:", ncol(unresolved), "\n")
cat("\nCluster counts:\n")

print(
  sort(
    table(
      unresolved[[cluster_col]][, 1]
    )
  )
)

cell_cluster_map <- data.table(
  cell = colnames(unresolved),
  cluster = as.character(
    unresolved[[cluster_col]][, 1]
  )
)

# ==============================================================================
# 4. BUILD STRICT ADENOMA RAW-COUNT PSEUDOBULKS
# ==============================================================================
#
# Formal DE uses:
#   GSE201348 Adenoma/Polyp source raw counts
#   GSE261388 P3/P4 only
#
# GSE261388 P1/P12 are inflammatory polyps and are excluded.
#
# CELLxGENE Chen cells remain useful for integrated visualization/annotation but
# are not included in this raw-count pseudobulk construction.
# ==============================================================================

source_files <- list(
  GSE201348 = find_source_file("GSE201348", source_dir),
  GSE261388 = find_source_file("GSE261388", source_dir)
)

cat("\nSource files used for formal pseudobulk DE:\n")
print(unlist(source_files))

study_pb_matrices <- list()
study_pb_meta <- list()

for (study in names(source_files)) {

  cat(
    "\n============================================================\n",
    "Loading raw-count source: ", study, "\n",
    "============================================================\n",
    sep = ""
  )

  src <- readRDS(source_files[[study]])
  counts <- get_joined_raw_counts(src, assay = "RNA")

  src_meta <- as.data.table(
    src[[]],
    keep.rownames = "cell"
  )

  overlap_cells <- intersect(
    colnames(counts),
    cell_cluster_map$cell
  )

  cat(
    "Cells overlapping existing reclustered object:",
    length(overlap_cells),
    "\n"
  )

  if (length(overlap_cells) == 0) {
    stop("No cell overlap for study ", study)
  }

  cell_info <- cell_cluster_map[cell %in% overlap_cells]

  if (study == "GSE201348") {

    if (!"geo_sample_id" %in% colnames(src_meta)) {
      stop("GSE201348 source metadata does not contain geo_sample_id.")
    }

    sample_map <- src_meta[
      cell %in% overlap_cells,
      .(
        cell,
        sample = as.character(geo_sample_id)
      )
    ]

    cell_info <- merge(
      cell_info,
      sample_map,
      by = "cell",
      all.x = TRUE,
      sort = FALSE
    )

    cell_info[, study := "GSE201348"]

  } else if (study == "GSE261388") {

    cell_info[
      ,
      sample := sub(
        "^GSE261388_([^_]+)_.*$",
        "\\1",
        cell
      )
    ]

    cell_info[, study := "GSE261388"]

    cat("\nGSE261388 sample counts before strict Adenoma filtering:\n")
    print(
      cell_info[
        ,
        .N,
        by = sample
      ][
        order(sample)
      ]
    )

    cell_info <- cell_info[
      sample %in% c("P3", "P4")
    ]

    cat("\nGSE261388 retained for formal DE (P3/P4 only):\n")
    print(
      cell_info[
        ,
        .N,
        by = sample
      ][
        order(sample)
      ]
    )
  }

  if (any(is.na(cell_info$sample) | cell_info$sample == "")) {
    stop("Missing biological sample IDs detected for study ", study)
  }

  cell_info[
    ,
    replicate := paste0(study, "__", sample)
  ]

  cell_info[
    ,
    group_id := paste0(
      replicate,
      "__CLUSTER__",
      cluster
    )
  ]

  pb_meta_study <- cell_info[
    ,
    .(
      study = first(study),
      sample = first(sample),
      replicate = first(replicate),
      cluster = first(cluster),
      n_cells = .N
    ),
    by = group_id
  ]

  pb_matrix <- aggregate_counts_by_group(
    counts = counts,
    cell_table = cell_info,
    group_col = "group_id"
  )

  pb_meta_study <- pb_meta_study[
    match(
      colnames(pb_matrix),
      group_id
    )
  ]

  stopifnot(
    identical(
      colnames(pb_matrix),
      pb_meta_study$group_id
    )
  )

  study_pb_matrices[[study]] <- pb_matrix
  study_pb_meta[[study]] <- pb_meta_study

  rm(src, counts, src_meta, cell_info, pb_matrix)
  invisible(gc())
}

# ==============================================================================
# 5. COMBINE STUDIES ON COMMON GENES
# ==============================================================================

common_genes <- Reduce(
  intersect,
  lapply(study_pb_matrices, rownames)
)

cat(
  "\nCommon genes across raw-count source studies:",
  length(common_genes),
  "\n"
)

study_pb_matrices <- lapply(
  study_pb_matrices,
  function(x) {
    x[common_genes, , drop = FALSE]
  }
)

pb_counts <- do.call(
  cbind,
  study_pb_matrices
)

pb_meta <- rbindlist(
  study_pb_meta,
  use.names = TRUE,
  fill = TRUE
)

pb_meta <- pb_meta[
  match(
    colnames(pb_counts),
    group_id
  )
]

stopifnot(
  identical(
    colnames(pb_counts),
    pb_meta$group_id
  )
)

if (
  any(
    pb_meta$study == "GSE261388" &
      pb_meta$sample %in% c("P1", "P12")
  )
) {
  stop(
    "P1/P12 inflammatory polyps are still present in formal DE pseudobulks."
  )
}

cat(
  "\nFinal raw-count pseudobulk matrix:",
  nrow(pb_counts),
  "genes x",
  ncol(pb_counts),
  "pseudobulks\n"
)

cat("\nGSE261388 samples retained:\n")
print(
  unique(
    pb_meta[
      study == "GSE261388",
      sample
    ]
  )
)

saveRDS(
  pb_counts,
  file.path(
    pb_dir,
    "ADENOMA_ONLY_raw_pseudobulk_counts.rds"
  ),
  compress = FALSE
)

fwrite(
  pb_meta,
  file.path(
    pb_dir,
    "ADENOMA_ONLY_pseudobulk_metadata.csv"
  )
)

# ==============================================================================
# 6. PRIMARY FORMAL DE:
#    DESeq2 paired pseudobulk LRT, each cluster vs all other unresolved clusters
# ==============================================================================

clusters <- sort(
  unique(
    as.numeric(
      as.character(
        pb_meta$cluster
      )
    )
  )
)

clusters <- as.character(clusters)

deseq_summary <- rbindlist(
  lapply(
    clusters,
    function(cl) {
      run_deseq2_lrt_cluster(
        pb_counts = pb_counts,
        pb_meta = pb_meta,
        target_cluster = cl,
        output_dir = pb_dir,
        min_target_cells = min_cells_target,
        min_rest_cells = min_cells_rest,
        min_count = min_count_per_sample,
        min_count_samples = min_samples_with_counts,
        padj_cutoff = padj_cutoff,
        log2fc_cutoff = log2fc_cutoff
      )
    }
  ),
  fill = TRUE
)

fwrite(
  deseq_summary,
  file.path(
    pb_dir,
    "ADENOMA_ONLY_DESeq2_LRT_run_summary.csv"
  )
)

cat(
  "\n============================================================\n",
  "DESeq2 RUN SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(deseq_summary)

# ==============================================================================
# 7. CLUSTER 1 RESCUE:
#    A. descriptive pooled log2((CPM_cluster1 + 1)/(CPM_rest + 1))
#    B. paired biological-sample Wilcoxon on log2(CPM + 1)
# ==============================================================================

cluster1_pair <- build_target_rest(
  pb_counts = pb_counts,
  pb_meta = pb_meta,
  target_cluster = cluster1_id,
  min_target_cells = min_cells_target,
  min_rest_cells = min_cells_rest
)

if (is.null(cluster1_pair)) {
  stop("Cluster 1 has no eligible biological-sample pairs.")
}

cat(
  "\nCluster 1 paired biological samples:",
  cluster1_pair$n_pairs,
  "\n"
)

# ------------------------------------------------------------------------------
# 7A. DESCRIPTIVE POOLED CPM RATIO
# ------------------------------------------------------------------------------

cluster1_cols <- which(
  as.character(pb_meta$cluster) == cluster1_id
)

cluster1_rest_cols <- which(
  as.character(pb_meta$cluster) != cluster1_id
)

cluster1_sum <- Matrix::rowSums(
  pb_counts[, cluster1_cols, drop = FALSE]
)

cluster1_rest_sum <- Matrix::rowSums(
  pb_counts[, cluster1_rest_cols, drop = FALSE]
)

cluster1_cpm <- (
  cluster1_sum /
    sum(cluster1_sum)
) * 1e6

cluster1_rest_cpm <- (
  cluster1_rest_sum /
    sum(cluster1_rest_sum)
) * 1e6

cluster1_cpm_ratio <- data.table(
  gene = rownames(pb_counts),
  cluster1_raw_count = as.numeric(cluster1_sum),
  rest_raw_count = as.numeric(cluster1_rest_sum),
  cluster1_CPM = as.numeric(cluster1_cpm),
  rest_CPM = as.numeric(cluster1_rest_cpm)
)

cluster1_cpm_ratio[
  ,
  log2_CPM_ratio :=
    log2(
      (cluster1_CPM + 1) /
        (rest_CPM + 1)
    )
]

cluster1_cpm_ratio <- cluster1_cpm_ratio[
  cluster1_raw_count >= 10
]

setorder(cluster1_cpm_ratio, -log2_CPM_ratio)

fwrite(
  cluster1_cpm_ratio,
  file.path(
    cluster1_dir,
    "Adenoma_cluster1_VS_rest_log2_CPM_ratio_ALL.csv"
  )
)

fwrite(
  head(cluster1_cpm_ratio, 100),
  file.path(
    cluster1_dir,
    "Adenoma_cluster1_VS_rest_log2_CPM_ratio_TOP100.csv"
  )
)

# ------------------------------------------------------------------------------
# 7B. PAIRED SAMPLE-LEVEL WILCOXON
# ------------------------------------------------------------------------------

target_counts <- cluster1_pair$target_counts
rest_counts <- cluster1_pair$rest_counts

keep_gene <- (
  rowSums(target_counts >= min_count_per_sample) >=
    min_samples_with_counts &
  rowSums(rest_counts >= min_count_per_sample) >=
    min_samples_with_counts
)

target_counts <- target_counts[
  keep_gene,
  ,
  drop = FALSE
]

rest_counts <- rest_counts[
  keep_gene,
  ,
  drop = FALSE
]

target_lib <- colSums(target_counts)
rest_lib <- colSums(rest_counts)

target_cpm <- sweep(
  target_counts,
  2,
  target_lib,
  "/"
) * 1e6

rest_cpm <- sweep(
  rest_counts,
  2,
  rest_lib,
  "/"
) * 1e6

target_logcpm <- log2(target_cpm + 1)
rest_logcpm <- log2(rest_cpm + 1)

genes <- rownames(target_logcpm)

cluster1_wilcox <- rbindlist(
  lapply(
    seq_along(genes),
    function(i) {
      target_values <- as.numeric(target_logcpm[i, ])
      rest_values <- as.numeric(rest_logcpm[i, ])
      difference <- target_values - rest_values

      test <- tryCatch(
        wilcox.test(
          target_values,
          rest_values,
          paired = TRUE,
          exact = FALSE
        ),
        error = function(e) NULL
      )

      data.table(
        gene = genes[i],
        n_pairs = length(target_values),
        median_target_log2CPM = median(
          target_values,
          na.rm = TRUE
        ),
        median_rest_log2CPM = median(
          rest_values,
          na.rm = TRUE
        ),
        median_log2CPM_difference = median(
          difference,
          na.rm = TRUE
        ),
        mean_log2CPM_difference = mean(
          difference,
          na.rm = TRUE
        ),
        raw_pvalue = if (is.null(test)) {
          NA_real_
        } else {
          test$p.value
        }
      )
    }
  )
)

cluster1_wilcox[
  ,
  padj := p.adjust(
    raw_pvalue,
    method = "BH"
  )
]

setorder(
  cluster1_wilcox,
  raw_pvalue,
  -median_log2CPM_difference
)

cluster1_wilcox_up <- cluster1_wilcox[
  !is.na(padj) &
    padj < 0.05 &
    median_log2CPM_difference > 0
]

setorder(
  cluster1_wilcox_up,
  -median_log2CPM_difference
)

fwrite(
  cluster1_wilcox,
  file.path(
    cluster1_dir,
    "Adenoma_cluster1_VS_rest_paired_Wilcoxon_ALL.csv"
  )
)

fwrite(
  cluster1_wilcox_up,
  file.path(
    cluster1_dir,
    "Adenoma_cluster1_VS_rest_paired_Wilcoxon_FDR_UP.csv"
  )
)

fwrite(
  cluster1_pair$support,
  file.path(
    cluster1_dir,
    "Adenoma_cluster1_paired_Wilcoxon_sample_support.csv"
  )
)

cat(
  "\nCluster 1 paired-Wilcoxon positive markers, padj < 0.05:",
  nrow(cluster1_wilcox_up),
  "\n"
)

print(
  head(
    cluster1_wilcox_up,
    50
  )
)

# ==============================================================================
# 8. CLUSTERS 9 AND 14:
#    CELL-LEVEL SEURAT WILCOXON FOR ANNOTATION ONLY
# ==============================================================================
#
# These clusters were overwhelmingly represented by GSM6061707.
# Therefore these p-values are NOT interpreted as independent-biological-
# replicate inference. They are used for marker ranking/annotation.
#
# No reclustering is performed here.
# ==============================================================================

DefaultAssay(unresolved) <- "RNA"

unresolved[["RNA"]] <- JoinLayers(
  unresolved[["RNA"]]
)

Idents(unresolved) <- cluster_col

cell_wilcox_results <- list()

for (cl in annotation_only_wilcox_clusters) {

  cat(
    "\n============================================================\n",
    "Cell-level Wilcoxon marker discovery: cluster ", cl, " vs all\n",
    "============================================================\n",
    sep = ""
  )

  res <- FindMarkers(
    unresolved,
    ident.1 = cl,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0,
    only.pos = FALSE,
    verbose = TRUE
  )

  res <- as.data.table(
    res,
    keep.rownames = "gene"
  )

  setorder(
    res,
    -avg_log2FC
  )

  sig_up <- res[
    p_val_adj < 0.05 &
      avg_log2FC > 0.25
  ]

  setorder(
    sig_up,
    -avg_log2FC
  )

  fwrite(
    res,
    file.path(
      wilcox_9_14_dir,
      paste0(
        "Adenoma_cluster",
        cl,
        "_VS_ALL_Wilcoxon_ALL.csv"
      )
    )
  )

  fwrite(
    sig_up,
    file.path(
      wilcox_9_14_dir,
      paste0(
        "Adenoma_cluster",
        cl,
        "_VS_ALL_Wilcoxon_UP.csv"
      )
    )
  )

  cell_wilcox_results[[cl]] <- list(
    all = res,
    up = sig_up
  )

  cat(
    "Significant positive marker rows:",
    nrow(sig_up),
    "\n"
  )

  print(
    head(
      sig_up,
      50
    )
  )
}

# ==============================================================================
# 9. FINAL ANNOTATION MAP
# ==============================================================================

adenoma_recluster_annotation <- c(
  "0" = "LGR5+ stem/progenitor-like epithelial cells",
  "1" = "CA1/CA2+ differentiated absorptive colonocyte-like cells",
  "2" = "Goblet/secretory epithelial cells",
  "3" = "Cycling/proliferating cells",
  "4" = "Inflammatory/regenerative crypt epithelial cells",
  "5" = "Mature absorptive colonocytes",
  "6" = "Activated/plasmablast-like B cells",
  "7" = "Smooth muscle/myofibroblast-like stromal cells",
  "8" = "POU2F3+ tuft-like/tuft-precursor cells",
  "9" = "OLFM4+ intestinal stem/progenitor epithelial cells",
  "10" = "WNT-high intestinal stem/progenitor epithelial cells",
  "11" = "Tuft cells",
  "12" = "BEST4+ OTOP2+ absorptive epithelial cells",
  "13" = "G2/M-high proliferating cells",
  "14" = "Cycling OLFM4+ stem/progenitor epithelial cells"
)

cluster_ids <- as.character(
  unresolved[[cluster_col]][, 1]
)

annotation_values <- unname(
  adenoma_recluster_annotation[
    cluster_ids
  ]
)

if (
  length(annotation_values) != ncol(unresolved)
) {
  stop(
    "Annotation-vector length does not match number of cells."
  )
}

if (any(is.na(annotation_values))) {
  missing_clusters <- unique(
    cluster_ids[
      is.na(annotation_values)
    ]
  )

  stop(
    paste(
      "Missing final annotations for clusters:",
      paste(missing_clusters, collapse = ", ")
    )
  )
}

unresolved$Adenoma_recluster_annotation <-
  annotation_values

unresolved$Adenoma_cluster_label <- paste0(
  cluster_ids,
  ": ",
  unresolved$Adenoma_recluster_annotation
)

annotation_check <- unique(
  data.table(
    cluster = cluster_ids,
    annotation =
      unresolved$Adenoma_recluster_annotation
  )
)

annotation_check[
  ,
  cluster_numeric := as.numeric(cluster)
]

setorder(
  annotation_check,
  cluster_numeric
)

annotation_check[
  ,
  legend_label := paste0(
    cluster,
    ": ",
    annotation
  )
]

legend_levels <- annotation_check$legend_label

unresolved$Adenoma_cluster_label <- factor(
  unresolved$Adenoma_cluster_label,
  levels = legend_levels
)

fwrite(
  annotation_check[
    ,
    .(
      cluster,
      annotation
    )
  ],
  file.path(
    adenoma_dir,
    "Adenoma_FINAL_recluster_annotation_map.csv"
  )
)

cat(
  "\n============================================================\n",
  "FINAL ANNOTATIONS\n",
  "============================================================\n",
  sep = ""
)

print(
  annotation_check[
    ,
    .(
      cluster,
      annotation
    )
  ],
  nrows = Inf
)

# ==============================================================================
# 10. FINAL UMAP
# ==============================================================================
#
# The axes are UMAP dimensions generated from the CCA-integrated representation.
# They are relabeled as UMAP 1 / UMAP 2 for presentation.
# ==============================================================================

available_reductions <- Reductions(unresolved)

if (
  "umap.cca.unresolved" %in%
    available_reductions
) {
  umap_use <- "umap.cca.unresolved"
} else if (
  "umap.cca" %in%
    available_reductions
) {
  umap_use <- "umap.cca"
} else if (
  "umap" %in%
    available_reductions
) {
  umap_use <- "umap"
} else {
  stop("No UMAP reduction found.")
}

cat(
  "\nUsing UMAP reduction:",
  umap_use,
  "\n"
)

umap_df <- as.data.table(
  Embeddings(
    unresolved,
    reduction = umap_use
  ),
  keep.rownames = "cell"
)

setnames(
  umap_df,
  colnames(umap_df)[2:3],
  c("UMAP_1", "UMAP_2")
)

meta_df <- data.table(
  cell = colnames(unresolved),
  cluster = cluster_ids,
  annotation =
    unresolved$Adenoma_recluster_annotation,
  cluster_label =
    unresolved$Adenoma_cluster_label
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
  by = cluster
]

cluster_centers[
  ,
  cluster_numeric := as.numeric(cluster)
]

setorder(
  cluster_centers,
  cluster_numeric
)

p_final <- DimPlot(
  unresolved,
  reduction = umap_use,
  group.by = "Adenoma_cluster_label",
  raster = TRUE,
  label = FALSE
) +
  geom_text(
    data = cluster_centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = cluster
    ),
    inherit.aes = FALSE,
    size = 5,
    fontface = "bold"
  ) +
  labs(
    title = "Adenoma/Polyp Final Annotation",
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

print(p_final)

ggsave(
  filename = file.path(
    plot_dir,
    "Adenoma_FINAL_ANNOTATED_UMAP.pdf"
  ),
  plot = p_final,
  width = 14,
  height = 9
)

ggsave(
  filename = file.path(
    plot_dir,
    "Adenoma_FINAL_ANNOTATED_UMAP.png"
  ),
  plot = p_final,
  width = 14,
  height = 9,
  dpi = 300
)

# ==============================================================================
# 11. SAVE FINAL ANNOTATED OBJECT
# ==============================================================================

final_object_file <- file.path(
  adenoma_dir,
  "Adenoma_polyp_UNRESOLVED_RECLUSTERED_FINAL_ANNOTATED.rds"
)

saveRDS(
  unresolved,
  final_object_file,
  compress = FALSE
)

cat(
  "\nFinal annotated object saved to:\n",
  final_object_file,
  "\n"
)

# ==============================================================================
# 12. SAVE SESSION INFORMATION
# ==============================================================================

capture.output(
  sessionInfo(),
  file = file.path(
    adenoma_dir,
    "Adenoma_DE_annotation_sessionInfo.txt"
  )
)

cat("\nWorkflow complete.\n")
