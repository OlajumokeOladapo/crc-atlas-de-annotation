############################################################
# OSCER R GUI
# CRC PSEUDOBULK: RERUN CLUSTERS 2, 5, 6, 8, 10
#
# Method: voomLmFit + limma eBayes
# Does NOT run estimateDisp() or glmLRT().
#
# Existing normalized checkpoints are reused.
# Missing checkpoints are prepared from the raw cache.
# Every cluster is refitted into a NEW batch folder.
############################################################


# ==========================================================
# 1. PACKAGES AND PATHS
# ==========================================================

required_packages <- c(
  "edgeR",
  "limma",
  "data.table",
  "Matrix"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(data.table)
  library(Matrix)
})

data.table::setDTthreads(1)

if (!"voomLmFit" %in% getNamespaceExports("edgeR")) {
  stop("Your edgeR installation does not provide voomLmFit().")
}

clusters_to_run <- c("2", "5", "6", "8", "10")

output_root <- paste0(
  "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
  "edgeR_pseudobulk_CRC_unresolved_reclusters_res0.5_low_memory"
)

cache_file <- paste0(
  "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
  "DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5/",
  "validated_pseudobulk_cache/",
  "CRC_unresolved_reclusters_res0.5_raw_pseudobulk.rds"
)

settings <- list(
  min_cells_target = 20,
  min_cells_rest = 20,
  min_replicates = 3,
  min_gene_total_count = 10,
  min_count_per_sample = 10,
  min_samples_expressed = 3,
  fdr = 0.05,
  log2fc = 0.25,
  strong_log2fc = 1
)

cat("edgeR:", as.character(packageVersion("edgeR")), "\n")
cat("limma:", as.character(packageVersion("limma")), "\n")


# ==========================================================
# 2. LOGGING AND CHECKPOINT HELPERS
# ==========================================================

log_stage <- function(label, log_file) {

  memory <- gc()

  line <- sprintf(
    "%s | %s | R heap %.1f MiB",
    Sys.time(),
    label,
    sum(memory[, 2])
  )

  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)

  flush.console()
}

save_checkpoint <- function(object, path) {

  temporary_path <- paste0(path, ".tmp")

  saveRDS(
    object,
    temporary_path,
    compress = FALSE
  )

  if (!file.rename(temporary_path, path)) {
    stop("Could not finalize checkpoint: ", path)
  }
}


# ==========================================================
# 3. PREPARE A MISSING NORMALIZED CHECKPOINT
#
# This retains the original:
# - study/donor/replicate pairing
# - target-versus-rest aggregation
# - expression filtering
# - TMM normalization
############################################################

prepare_cluster <- function(
    cluster_id,
    cache_file,
    cluster_dir,
    cfg = settings) {

  dir.create(
    cluster_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  log_file <- file.path(cluster_dir, "progress.log")

  log_stage("Loading raw pseudobulk cache", log_file)

  pb <- readRDS(cache_file)

  if (
    !is.list(pb$parts) ||
    !length(pb$parts) ||
    is.null(names(pb$parts)) ||
    anyDuplicated(names(pb$parts))
  ) {
    stop("Cache must contain uniquely named pb$parts.")
  }

  keys <- c("study", "donor", "replicate")

  metadata_list <- lapply(
    names(pb$parts),
    function(study_id) {

      md <- data.table::copy(
        data.table::as.data.table(
          pb$parts[[study_id]]$metadata
        )
      )

      needed <- c(
        keys,
        "cluster",
        "pseudobulk_id",
        "n_cells"
      )

      if (!all(needed %in% names(md))) {
        stop("Metadata fields missing in ", study_id)
      }

      md <- md[, ..needed]

      identifier_fields <- c(
        keys,
        "cluster",
        "pseudobulk_id"
      )

      for (field in identifier_fields) {
        data.table::set(
          md,
          j = field,
          value = as.character(md[[field]])
        )
      }

      if (
        anyNA(md) ||
        any(md$study != study_id) ||
        any(vapply(
          md[, identifier_fields, with = FALSE],
          function(x) any(!nzchar(x)),
          logical(1)
        ))
      ) {
        stop("Missing or inconsistent metadata in ", study_id)
      }

      if (
        !is.numeric(md$n_cells) ||
        any(!is.finite(md$n_cells)) ||
        any(md$n_cells < 0)
      ) {
        stop("Invalid cell counts in ", study_id)
      }

      if (anyDuplicated(md$pseudobulk_id)) {
        stop("Duplicate pseudobulk IDs in ", study_id)
      }

      count_matrix <- pb$parts[[study_id]]$counts

      if (
        is.null(rownames(count_matrix)) ||
        is.null(colnames(count_matrix)) ||
        anyDuplicated(rownames(count_matrix)) ||
        anyDuplicated(colnames(count_matrix)) ||
        !all(md$pseudobulk_id %in% colnames(count_matrix))
      ) {
        stop("Count/metadata identifier problem in ", study_id)
      }

      md
    }
  )

  md <- data.table::rbindlist(metadata_list)

  rm(metadata_list)

  if (!cluster_id %in% md$cluster) {
    stop("Cluster absent from cache: ", cluster_id)
  }

  eligible <- md[
    ,
    .(
      target_cells = sum(n_cells[cluster == cluster_id]),
      rest_cells = sum(n_cells[cluster != cluster_id])
    ),
    by = keys
  ]

  eligible[
    ,
    eligible := (
      target_cells >= cfg$min_cells_target &
      rest_cells >= cfg$min_cells_rest
    )
  ]

  data.table::fwrite(
    eligible,
    file.path(cluster_dir, "replicate_eligibility.csv")
  )

  eligible <- eligible[eligible == TRUE]

  if (nrow(eligible) < cfg$min_replicates) {
    stop("Insufficient eligible paired replicates.")
  }

  eligible[, pair_id := paste0("pair", seq_len(.N))]

  studies <- unique(eligible$study)

  common_genes <- Reduce(
    intersect,
    lapply(
      studies,
      function(study_id) {
        rownames(pb$parts[[study_id]]$counts)
      }
    )
  )

  if (!length(common_genes)) {
    stop("No common genes across participating studies.")
  }

  writeLines(
    common_genes,
    file.path(cluster_dir, "participating_common_genes.txt")
  )

  cd <- as.data.frame(
    eligible[rep(seq_len(.N), each = 2)]
  )

  cd$eligible <- NULL

  cd$group <- factor(
    rep(c("rest", "target"), nrow(eligible)),
    levels = c("rest", "target")
  )

  cd$n_cells <- ifelse(
    cd$group == "target",
    cd$target_cells,
    cd$rest_cells
  )

  rownames(cd) <- paste(
    cd$pair_id,
    cd$group,
    sep = "_"
  )

  counts <- matrix(
    0,
    nrow = length(common_genes),
    ncol = nrow(cd),
    dimnames = list(
      common_genes,
      rownames(cd)
    )
  )

  log_stage(
    sprintf(
      "Aggregating %d genes across %d pseudobulk libraries",
      nrow(counts),
      ncol(counts)
    ),
    log_file
  )

  # Aggregate in blocks and release each study after use.
  for (study_id in studies) {

    part <- pb$parts[[study_id]]

    gene_rows <- match(
      common_genes,
      rownames(part$counts)
    )

    study_metadata <- md[study == study_id]

    for (j in which(cd$study == study_id)) {

      donor_id <- cd$donor[j]
      replicate_id <- cd$replicate[j]
      is_target <- cd$group[j] == "target"

      ids <- study_metadata[
        donor == donor_id &
          replicate == replicate_id &
          ((cluster == cluster_id) == is_target),
        pseudobulk_id
      ]

      if (!length(ids)) {
        stop("No input pseudobulks for ", rownames(cd)[j])
      }

      columns <- match(
        ids,
        colnames(part$counts)
      )

      for (
        first_row in seq.int(
          1L,
          length(gene_rows),
          by = 1000L
        )
      ) {

        block_rows <- seq.int(
          first_row,
          min(first_row + 999L, length(gene_rows))
        )

        block <- part$counts[
          gene_rows[block_rows],
          columns,
          drop = FALSE
        ]

        values <- if (inherits(block, "sparseMatrix")) {
          block@x
        } else {
          as.vector(block)
        }

        if (
          any(!is.finite(values)) ||
          any(values < 0) ||
          any(abs(values - round(values)) > 1e-7)
        ) {
          stop("Expected finite, nonnegative, raw integer counts.")
        }

        counts[block_rows, j] <- Matrix::rowSums(block)
      }
    }

    pb$parts[[study_id]] <- NULL

    rm(part, study_metadata, block, values)
    invisible(gc())
  }

  rm(pb, md)
  invisible(gc())

  log_stage("Raw cache released; filtering genes", log_file)

  totals <- numeric(nrow(counts))
  expressed <- integer(nrow(counts))

  for (j in seq_len(ncol(counts))) {

    totals <- totals + counts[, j]

    expressed <- expressed + (
      counts[, j] >= cfg$min_count_per_sample
    )
  }

  keep <- (
    totals >= cfg$min_gene_total_count &
    expressed >= cfg$min_samples_expressed
  )

  counts <- counts[keep, , drop = FALSE]

  n_initial <- nrow(counts)

  if (
    !n_initial ||
    any(colSums(counts) <= 0)
  ) {
    stop("No retained genes or a zero-sized library.")
  }

  cd$pair_id <- factor(
    cd$pair_id,
    levels = eligible$pair_id
  )

  design <- model.matrix(
    ~ pair_id + group,
    data = cd
  )

  if (
    qr(design)$rank != ncol(design) ||
    nrow(design) <= ncol(design)
  ) {
    stop("Invalid paired design.")
  }

  data.table::fwrite(
    data.table::data.table(
      pseudobulk_id = rownames(cd),
      cd
    ),
    file.path(cluster_dir, "edgeR_coldata.csv")
  )

  y <- edgeR::DGEList(counts)

  rm(counts, totals, expressed, keep)
  invisible(gc())

  keep <- edgeR::filterByExpr(
    y,
    design = design
  )

  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]

  if (
    !nrow(y) ||
    any(y$samples$lib.size <= 0)
  ) {
    stop("No genes or a zero library after filterByExpr.")
  }

  y <- edgeR::calcNormFactors(
    y,
    method = "TMM"
  )

  info <- list(
    cluster = cluster_id,
    paired_replicates = nrow(eligible),
    n_studies = length(studies),
    target_cells = sum(eligible$target_cells),
    rest_cells = sum(eligible$rest_cells),
    common_genes = length(common_genes),
    genes_after_initial_filter = n_initial,
    genes_after_filterByExpr = nrow(y)
  )

  state <- list(
    stage = "normalized",
    y = y,
    design = design,
    info = info,
    settings = cfg,
    cluster_dir = normalizePath(cluster_dir),
    edgeR_version = as.character(packageVersion("edgeR"))
  )

  save_checkpoint(
    state,
    file.path(cluster_dir, "normalized_checkpoint.rds")
  )

  log_stage("Normalized checkpoint saved", log_file)

  invisible(NULL)
}


# ==========================================================
# 4. FIT LIMMA-VOOM FROM A NORMALIZED CHECKPOINT
# ==========================================================

fit_crc_voom_checkpoint <- function(
    checkpoint_file,
    output_dir,
    expected_cluster) {

  state <- readRDS(checkpoint_file)

  if (!state$stage %in% c("normalized", "dispersions")) {
    stop("Unsupported checkpoint stage.")
  }

  y <- state$y
  design <- state$design
  info <- state$info
  cfg <- state$settings

  if (as.character(info$cluster) != expected_cluster) {
    stop("Checkpoint cluster ID does not match requested cluster.")
  }

  if (!inherits(y, "DGEList")) {
    stop("Checkpoint y must be a DGEList.")
  }

  if (!identical(
    rownames(design),
    colnames(y$counts)
  )) {
    stop("Design rows and count columns do not match.")
  }

  if (
    qr(design)$rank != ncol(design) ||
    nrow(design) <= ncol(design)
  ) {
    stop("Invalid paired design.")
  }

  if (!"grouptarget" %in% colnames(design)) {
    stop("Expected coefficient 'grouptarget' is missing.")
  }

  if (
    any(!is.finite(y$counts)) ||
    any(y$counts < 0)
  ) {
    stop("Invalid counts.")
  }

  effective_library <- (
    y$samples$lib.size *
    y$samples$norm.factors
  )

  if (
    any(!is.finite(effective_library)) ||
    any(effective_library <= 0)
  ) {
    stop("Invalid effective library sizes.")
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  started <- proc.time()[["elapsed"]]
  log_file <- file.path(output_dir, "progress.log")

  progress <- function(label) {

    line <- sprintf(
      "%s | %.1f seconds | %s",
      Sys.time(),
      proc.time()[["elapsed"]] - started,
      label
    )

    cat(line, "\n")
    cat(line, "\n", file = log_file, append = TRUE)

    flush.console()
  }

  capture.output(
    sessionInfo(),
    file = file.path(output_dir, "sessionInfo.txt")
  )

  progress(sprintf(
    "START: cluster %s; %d genes; %d libraries; %d coefficients",
    info$cluster,
    nrow(y),
    ncol(y),
    ncol(design)
  ))

  saveRDS(
    list(
      source_checkpoint = normalizePath(checkpoint_file),
      settings = cfg,
      info = info,
      design = design,
      source_edgeR_version = state$edgeR_version,
      current_edgeR_version = as.character(
        packageVersion("edgeR")
      ),
      current_limma_version = as.character(
        packageVersion("limma")
      ),
      method = "limma_voomLmFit_eBayes"
    ),
    file.path(output_dir, "analysis_settings.rds")
  )

  write.csv(
    data.frame(
      pseudobulk_id = colnames(y$counts),
      library_size = y$samples$lib.size,
      norm_factor = y$samples$norm.factors,
      effective_library_size = effective_library
    ),
    file.path(output_dir, "normalization.csv"),
    row.names = FALSE
  )

  rm(state)
  invisible(gc())

  # --------------------------------------------------------
  # Fit weighted linear models.
  # Existing TMM factors and paired design are retained.
  # --------------------------------------------------------

  progress("Starting voomLmFit")

  arguments <- list(
    counts = y,
    design = design,
    plot = FALSE,
    save.plot = TRUE
  )

  if (
    "keep.EList" %in%
    names(formals(edgeR::voomLmFit))
  ) {
    arguments$keep.EList <- FALSE
  }

  fit <- withCallingHandlers(
    do.call(
      edgeR::voomLmFit,
      arguments
    ),
    warning = function(w) {
      progress(paste(
        "WARNING:",
        conditionMessage(w)
      ))
    }
  )

  rm(y, arguments)
  invisible(gc())

  saveRDS(
    fit,
    file.path(output_dir, "voom_fit_checkpoint.rds"),
    compress = FALSE
  )

  progress("Voom fit checkpoint saved")

  # --------------------------------------------------------
  # Empirical-Bayes moderation
  # --------------------------------------------------------

  progress("Starting empirical-Bayes moderation")

  fit <- limma::eBayes(
    fit,
    trend = FALSE,
    robust = FALSE
  )

  saveRDS(
    fit,
    file.path(output_dir, "moderated_fit.rds"),
    compress = FALSE
  )

  tt <- limma::topTable(
    fit,
    coef = "grouptarget",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )

  results <- data.frame(
    gene = rownames(tt),
    tt,
    row.names = NULL,
    check.names = FALSE
  )

  results$cluster <- as.character(info$cluster)

  results$comparison <- paste0(
    "recluster_",
    info$cluster,
    "_VS_all_other_new_reclusters"
  )

  results$method <- "limma_voomLmFit_eBayes"
  results$FDR <- results$adj.P.Val

  valid <- (
    is.finite(results$FDR) &
    is.finite(results$logFC)
  )

  upregulated <- results[
    valid &
      results$FDR < cfg$fdr &
      results$logFC >= cfg$log2fc,
    ,
    drop = FALSE
  ]

  downregulated <- results[
    valid &
      results$FDR < cfg$fdr &
      results$logFC <= -cfg$log2fc,
    ,
    drop = FALSE
  ]

  strong_upregulated <- results[
    valid &
      results$FDR < cfg$fdr &
      results$logFC >= cfg$strong_log2fc,
    ,
    drop = FALSE
  ]

  tables <- list(
    all_genes = results,
    upregulated = upregulated,
    downregulated = downregulated,
    strong_upregulated = strong_upregulated,
    top100_upregulated = head(upregulated, 100),
    top100_strong_upregulated = head(
      strong_upregulated,
      100
    )
  )

  for (table_name in names(tables)) {

    write.csv(
      tables[[table_name]],
      file.path(
        output_dir,
        paste0(
          "limma_voom_",
          table_name,
          ".csv"
        )
      ),
      row.names = FALSE
    )
  }

  write.csv(
    data.frame(
      gene = rownames(fit$coefficients),
      residual_df = fit$df.residual
    ),
    file.path(output_dir, "gene_residual_df.csv"),
    row.names = FALSE
  )

  # --------------------------------------------------------
  # Save PDF and PNG diagnostics.
  # A graphics error does not discard completed statistics.
  # --------------------------------------------------------

  plot_errors <- character()

  if (
    !is.null(fit$voom.xy) &&
    !is.null(fit$voom.line)
  ) {

    for (format in c("pdf", "png")) {

      plot_error <- tryCatch({

        plot_file <- file.path(
          output_dir,
          paste0("voom_mean_variance.", format)
        )

        if (format == "pdf") {

          grDevices::pdf(
            plot_file,
            width = 8,
            height = 6
          )

        } else {

          grDevices::png(
            plot_file,
            width = 8,
            height = 6,
            units = "in",
            res = 300,
            type = if (capabilities("cairo")) {
              "cairo"
            } else {
              getOption("bitmapType", "Xlib")
            }
          )
        }

        tryCatch({

          plot(
            fit$voom.xy$x,
            fit$voom.xy$y,
            pch = 16,
            cex = 0.3,
            col = "#777777",
            xlab = fit$voom.xy$xlab,
            ylab = fit$voom.xy$ylab,
            main = paste(
              "Cluster",
              info$cluster,
              "- voom mean-variance trend"
            )
          )

          lines(
            fit$voom.line,
            col = "#B23A48",
            lwd = 2
          )

        }, finally = {
          grDevices::dev.off()
        })

        NULL

      }, error = function(e) {

        paste(
          format,
          "plot failed:",
          conditionMessage(e)
        )
      })

      if (!is.null(plot_error)) {

        plot_errors <- c(
          plot_errors,
          plot_error
        )

        progress(plot_error)
      }
    }
  }

  if (length(plot_errors)) {
    writeLines(
      plot_errors,
      file.path(output_dir, "plot_errors.txt")
    )
  }

  summary_table <- data.frame(
    cluster = info$cluster,
    method = "limma_voomLmFit_eBayes",
    genes_reported = nrow(results),
    genes_with_finite_results = sum(valid),
    paired_replicates = info$paired_replicates,
    n_studies = info$n_studies,
    significant_up = nrow(upregulated),
    significant_down = nrow(downregulated),
    strong_up = nrow(strong_upregulated),
    fdr_cutoff = cfg$fdr,
    log2fc_cutoff = cfg$log2fc,
    elapsed_seconds = (
      proc.time()[["elapsed"]] - started
    ),
    plot_errors = length(plot_errors),
    status = "completed"
  )

  write.csv(
    summary_table,
    file.path(
      output_dir,
      "limma_voom_cluster_summary.csv"
    ),
    row.names = FALSE
  )

  writeLines(
    c(
      "Method: edgeR::voomLmFit followed by limma::eBayes.",
      "Uses filtered raw counts and saved TMM normalization factors.",
      "The original paired design is retained.",
      "No estimateDisp or glmLRT is performed.",
      "Positive logFC means higher expression in the target cluster.",
      "BH adjustment is within each cluster comparison.",
      "Fold-change thresholds are post-test filters.",
      "Combined tables do not apply an additional cross-cluster FDR correction."
    ),
    file.path(output_dir, "METHOD.txt")
  )

  progress("COMPLETE")

  print(summary_table)

  cat(
    "\nCluster results:\n",
    output_dir,
    "\n",
    sep = ""
  )

  invisible(list(
    summary = summary_table,
    output_dir = output_dir
  ))
}


# ==========================================================
# 5. RUN ALL FIVE CLUSTERS SEQUENTIALLY
# ==========================================================

run_all_crc_clusters <- function(
    clusters,
    output_root,
    cache_file) {

  clusters <- unique(as.character(clusters))

  if (!length(clusters)) {
    stop("No clusters selected.")
  }

  dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
  )

  checkpoint_paths <- file.path(
    output_root,
    paste0("cluster_", clusters),
    "normalized_checkpoint.rds"
  )

  needs_preparation <- !file.exists(
    checkpoint_paths
  )

  if (
    any(needs_preparation) &&
    !file.exists(cache_file)
  ) {
    stop(
      "Missing normalized checkpoints for clusters: ",
      paste(
        clusters[needs_preparation],
        collapse = ", "
      ),
      "\nThe raw cache was not found at:\n",
      cache_file,
      "\nCorrect cache_file near the top of this script."
    )
  }

  batch_dir <- tempfile(
    pattern = paste0(
      "limma_voom_batch_",
      format(Sys.time(), "%Y%m%d_%H%M%S"),
      "_"
    ),
    tmpdir = output_root
  )

  dir.create(
    batch_dir,
    recursive = TRUE
  )

  capture.output(
    sessionInfo(),
    file = file.path(batch_dir, "sessionInfo.txt")
  )

  cat("\nBatch folder:\n", batch_dir, "\n", sep = "")
  cat("Clusters:", paste(clusters, collapse = ", "), "\n")

  status_rows <- list()
  summaries <- list()
  completed_dirs <- character()

  for (cluster_id in clusters) {

    cat(
      "\n========== CLUSTER ",
      cluster_id,
      " ==========\n",
      sep = ""
    )

    began <- proc.time()[["elapsed"]]

    cluster_dir <- file.path(
      output_root,
      paste0("cluster_", cluster_id)
    )

    dir.create(
      cluster_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    checkpoint <- file.path(
      cluster_dir,
      "normalized_checkpoint.rds"
    )

    results_dir <- file.path(
      batch_dir,
      paste0("cluster_", cluster_id)
    )

    outcome <- tryCatch({

      if (!file.exists(checkpoint)) {

        cat(
          "Preparing missing normalized checkpoint.\n",
          "No edgeR dispersion estimation will be run.\n",
          sep = ""
        )

        prepare_cluster(
          cluster_id = cluster_id,
          cache_file = cache_file,
          cluster_dir = cluster_dir
        )

        invisible(gc())

      } else {

        cat(
          "Reusing normalized checkpoint; ",
          "rerunning limma-voom.\n",
          sep = ""
        )
      }

      fitted <- fit_crc_voom_checkpoint(
        checkpoint_file = checkpoint,
        output_dir = results_dir,
        expected_cluster = cluster_id
      )

      list(
        ok = TRUE,
        summary = fitted$summary,
        error = ""
      )

    }, error = function(e) {

      cat(
        "Cluster ",
        cluster_id,
        " failed: ",
        conditionMessage(e),
        "\n",
        sep = ""
      )

      list(
        ok = FALSE,
        error = conditionMessage(e)
      )
    })

    status_rows[[cluster_id]] <- data.frame(
      cluster = cluster_id,
      status = if (outcome$ok) {
        "completed"
      } else {
        "failed"
      },
      elapsed_seconds = (
        proc.time()[["elapsed"]] - began
      ),
      results_dir = results_dir,
      error = outcome$error
    )

    write.csv(
      do.call(rbind, status_rows),
      file.path(batch_dir, "batch_status.csv"),
      row.names = FALSE
    )

    if (outcome$ok) {

      summaries[[cluster_id]] <- outcome$summary
      completed_dirs[cluster_id] <- results_dir

      write.csv(
        do.call(rbind, summaries),
        file.path(
          batch_dir,
          "combined_cluster_summary.csv"
        ),
        row.names = FALSE
      )
    }

    rm(outcome)
    invisible(gc())
  }

  # --------------------------------------------------------
  # Combine successful results from THIS batch only
  # --------------------------------------------------------

  if (length(completed_dirs)) {

    table_types <- c(
      "all_genes",
      "upregulated",
      "downregulated",
      "strong_upregulated",
      "top100_upregulated",
      "top100_strong_upregulated"
    )

    for (table_type in table_types) {

      paths <- file.path(
        completed_dirs,
        paste0(
          "limma_voom_",
          table_type,
          ".csv"
        )
      )

      combined <- data.table::rbindlist(
        lapply(paths, data.table::fread),
        use.names = TRUE,
        fill = TRUE
      )

      data.table::fwrite(
        combined,
        file.path(
          batch_dir,
          paste0(
            "combined_",
            table_type,
            ".csv"
          )
        )
      )

      rm(combined)
      invisible(gc())
    }
  }

  final_status <- do.call(
    rbind,
    status_rows
  )

  cat("\n========== BATCH SUMMARY ==========\n")

  print(
    final_status[
      ,
      c(
        "cluster",
        "status",
        "elapsed_seconds"
      )
    ],
    row.names = FALSE
  )

  cat(
    "\nCompleted ",
    sum(final_status$status == "completed"),
    " of ",
    length(clusters),
    " clusters.\n",
    sep = ""
  )

  cat(
    "\nAll new results are saved in:\n",
    batch_dir,
    "\n",
    sep = ""
  )

  if (any(final_status$status == "failed")) {
    cat(
      "\nSome clusters failed. ",
      "Check batch_status.csv for the error messages.\n",
      sep = ""
    )
  }

  invisible(list(
    batch_dir = batch_dir,
    status = final_status,
    summaries = summaries
  ))
}


# ==========================================================
# 6. START THE BATCH
# ==========================================================

batch_results <- run_all_crc_clusters(
  clusters = clusters_to_run,
  output_root = output_root,
  cache_file = cache_file
)