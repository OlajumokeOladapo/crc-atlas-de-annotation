############################################################
# CRC ORIGINAL CCA CLUSTER PSEUDOBULK DESeq2
#
# CRC-only version of the validated atlas-wide CCA
# pseudobulk DESeq2 workflow.
#
# Starting clusters:
#   cca_snn_res.0.5 in CRC_CCA_integrated.rds
#
# Comparison:
#   each original CRC CCA cluster vs all other CRC CCA clusters
#   within biological replicate
############################################################

required <- c("Seurat", "SeuratObject", "Matrix", "data.table", "DESeq2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
})

# Do not call library(DESeq2). DESeq2 is used with DESeq2:: only.

atlas_root <- Sys.getenv(
  "CRC_ATLAS_ROOT",
  unset = "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)
source_dir <- file.path(atlas_root, "full_seurat_QC_normalized_by_study_condition")
output_root <- file.path(atlas_root, "DESeq2_pseudobulk_CCA_res0.5_revised")
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

config <- list(
  cluster_column = "cca_snn_res.0.5",
  min_cells_target = 20L,
  min_cells_rest = 20L,
  min_replicates = 3L,
  min_gene_total_count = 10L,
  alpha = 0.05,
  log2fc_cutoff = 0.25,
  min_match_fraction = 0.99,
  block_cells = 25000L,
  integer_tolerance = 1e-8,
  save_dds = TRUE,
  use_cache = TRUE,
  cache_version = "CCA_DESeq2_final_donorfix_v1"
)

conditions <- "CRC"

expected_cells <- c(
  Adenoma_polyp = 306396L,
  CRC = 702029L,
  Normal = 151684L,
  TAN = 287894L
)

cca_files <- c(
  Adenoma_polyp = file.path(
    atlas_root, "Adenoma_Harmony_CCA_with_Chen", "Adenoma_CCA_integrated_WITH_CHEN.rds"
  ),
  CRC = file.path(
    atlas_root, "Integration_by_condition", "CRC", "CRC_CCA_integrated.rds"
  ),
  Normal = file.path(
    atlas_root, "Integration_by_condition", "Normal", "Normal_CCA_integrated.rds"
  ),
  TAN = file.path(
    atlas_root, "Integration_by_condition", "TAN", "TAN_CCA_integrated_clustered.rds"
  )
)

studies_by_condition <- list(
  Adenoma_polyp = c("CELLxGENE_Chen_2021", "GSE201348", "GSE261388"),
  CRC = c(
    "CELLxGENE_Chen_2021", "CELLxGENE_Moorman2025", "GSE132465", "GSE178341",
    "GSE200997", "GSE201348", "GSE221575", "GSE225857", "GSE231559",
    "GSE261388", "GSE294300", "GSE302903", "GSE330797"
  ),
  Normal = c("CELLxGENE_Chen_2021", "GSE201348", "GSE221575", "GSE231559", "GSE261388"),
  TAN = c(
    "CELLxGENE_Moorman2025", "GSE132465", "GSE178341", "GSE200997",
    "GSE225857", "GSE231559", "GSE294300"
  )
)

excluded_nonraw <- c("CELLxGENE_Chen_2021", "CELLxGENE_Moorman2025")

# Optional manual donor-column overrides, for example:
# donor_overrides <- c(GSE132465 = "Patient")
donor_overrides <- character()

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

blank <- function(x) is.na(x) | trimws(as.character(x)) == ""
safe_token <- function(x) gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
save_csv <- function(x, path) fwrite(as.data.table(x), path)

check_names <- function(x, label) {
  x <- as.character(x)
  assert(length(x) > 0, paste0(label, " is empty"))
  assert(!any(blank(x)), paste0(label, " contains missing/blank values"))
  assert(!anyDuplicated(x), paste0(label, " contains duplicates"))
  invisible(TRUE)
}

check_integer_counts <- function(mat) {
  x <- if (inherits(mat, "sparseMatrix")) mat@x else as.numeric(mat)
  if (!length(x)) return(invisible(TRUE))
  assert(all(is.finite(x)), "Count matrix contains non-finite values")
  assert(all(x >= 0), "Count matrix contains negative values")
  assert(
    all(abs(x - round(x)) <= config$integer_tolerance),
    "Count matrix is not raw integer-like counts"
  )
  invisible(TRUE)
}

source_file_for <- function(study, condition) {
  f <- file.path(source_dir, paste0(study, "_", condition, "_QC_normalized_seurat.rds"))
  assert(file.exists(f), paste0("Missing source file: ", f))
  f
}

find_donor_column <- function(md, study) {
  if (length(donor_overrides) && study %in% names(donor_overrides)) {
    d <- as.character(donor_overrides[[study]])
    assert(d %in% colnames(md), paste0(study, ": donor override column absent: ", d))
    assert(any(!blank(md[[d]])), paste0(study, ": donor override column has no usable values"))
    message("  selected donor column by manual override: ", d)
    return(d)
  }

  preferred <- c(
    "geo_patient_id.ch1", "geo_patient_id", "geo_patient", "geo_Patient",
    "patient_id", "Patient_ID", "Patient ID", "patient", "Patient",
    "donor_id", "Donor_ID", "donor", "Donor", "subject_id", "subject",
    "individual_id", "cxg_donor_id", "cxg_Patient"
  )

  for (d in preferred) {
    if (!d %in% colnames(md)) next
    v <- as.character(md[[d]])
    v <- v[!blank(v)]
    if (!length(v)) next
    n_unique <- length(unique(v))
    if (n_unique >= 1 && n_unique < nrow(md)) {
      message("  selected donor column: ", d, " (", n_unique, " unique biological units)")
      return(d)
    }
  }

  pattern_candidates <- colnames(md)[
    grepl(
      "(^|[._])patient([._]|$)|(^|[._])donor([._]|$)|(^|[._])subject([._]|$)|individual",
      colnames(md),
      ignore.case = TRUE
    )
  ]
  pattern_candidates <- pattern_candidates[
    !grepl(
      "type|status|stage|score|count|tumor|disease|sex|age|number$",
      pattern_candidates,
      ignore.case = TRUE
    )
  ]

  for (d in pattern_candidates) {
    v <- as.character(md[[d]])
    v <- v[!blank(v)]
    if (!length(v)) next
    n_unique <- length(unique(v))
    if (n_unique >= 1 && n_unique < nrow(md)) {
      message("  selected donor column by pattern: ", d, " (", n_unique, " unique biological units)")
      return(d)
    }
  }

  sample_candidates <- c(
    "geo_sample_id", "cxg_Sample.ID", "sample_id", "Sample_ID", "Sample ID",
    "sample", "Sample", "sample_name", "sample_accession", "library_id", "orig.ident"
  )

  for (d in sample_candidates) {
    if (!d %in% colnames(md)) next
    v <- as.character(md[[d]])
    v <- v[!blank(v)]
    if (!length(v)) next
    n_unique <- length(unique(v))
    if (n_unique >= 1 && n_unique < nrow(md)) {
      warning("No usable patient/donor identifier found. Falling back to sample-level column: ", d)
      return(d)
    }
  }

  stop(
    paste0(study, ": no donor/sample column found. Columns: ", paste(colnames(md), collapse = "; ")),
    call. = FALSE
  )
}

validate_cca <- function(condition) {
  f <- unname(cca_files[condition])
  assert(file.exists(f), paste0("Missing CCA file: ", f))
  message("Validating CCA object: ", f)

  obj <- readRDS(f)
  assert(inherits(obj, "Seurat"), paste0(condition, ": CCA object is not Seurat"))
  assert(
    ncol(obj) == unname(expected_cells[condition]),
    paste0(condition, ": CCA cell count differs from expected full atlas count")
  )

  md <- obj@meta.data
  check_names(rownames(md), paste0(condition, " CCA cell IDs"))
  assert(identical(rownames(md), colnames(obj)), paste0(condition, ": CCA metadata order mismatch"))
  assert(config$cluster_column %in% colnames(md), paste0(condition, ": CCA res0.5 column missing"))
  assert(!any(blank(md[[config$cluster_column]])), paste0(condition, ": missing cluster labels"))

  out <- md[, config$cluster_column, drop = FALSE]
  rm(obj, md)
  invisible(gc())
  out
}

aggregate_study <- function(study, condition, cca_md, out_dir) {
  message("Aggregating ", study)
  f <- source_file_for(study, condition)
  obj <- readRDS(f)

  assert(inherits(obj, "Seurat"), paste0(study, ": source object is not Seurat"))
  assert("RNA" %in% names(obj@assays), paste0(study, ": RNA assay absent"))
  message("  RNA assay confirmed")

  md <- obj@meta.data
  check_names(rownames(md), paste0(study, " metadata cell IDs"))
  assert(identical(rownames(md), colnames(obj)), paste0(study, ": metadata/cell order mismatch"))

  donor_col <- find_donor_column(md, study)
  message("  donor column: ", donor_col)

  rna <- obj@assays$RNA
  layers <- as.character(SeuratObject::Layers(rna))
  message("  RNA layers: ", paste(layers, collapse = ", "))

  count_layers <- layers[grepl("^counts($|\\.)", layers)]
  assert(length(count_layers) > 0, paste0(study, ": no counts layer found"))
  message("  count layers: ", paste(count_layers, collapse = ", "))

  layer_features <- list()
  source_cells <- character()
  seen_cells <- character()

  for (layer_name in count_layers) {
    mat <- SeuratObject::LayerData(object = rna, layer = layer_name)
    mat <- as(mat, "CsparseMatrix")

    genes <- as.character(rownames(mat))
    cells <- as.character(colnames(mat))
    check_names(genes, paste0(study, " genes in ", layer_name))
    check_names(cells, paste0(study, " cells in ", layer_name))

    overlap <- base::intersect(seen_cells, cells)
    assert(!length(overlap), paste0(study, ": cells occur in multiple counts layers"))

    layer_features[layer_name] <- list(genes)
    source_cells <- c(source_cells, cells)
    seen_cells <- c(seen_cells, cells)

    rm(mat)
    invisible(gc())
  }

  check_names(source_cells, paste0(study, " source cells"))
  assert(all(source_cells %in% rownames(md)), paste0(study, ": count-layer cells absent from metadata"))
  message("  source cells in counts layers: ", format(length(source_cells), big.mark = ","))

  cca_cells <- rownames(cca_md)
  cca_lookup <- stats::setNames(cca_cells, cca_cells)
  mapped_cca <- unname(cca_lookup[source_cells])
  matched <- !is.na(mapped_cca)
  match_fraction <- sum(matched) / length(source_cells)

  message(
    "  cells matched to CCA: ", format(sum(matched), big.mark = ","), " / ",
    format(length(source_cells), big.mark = ","), " (", sprintf("%.2f%%", 100 * match_fraction), ")"
  )

  assert(
    match_fraction >= config$min_match_fraction,
    paste0(study, ": CCA match fraction too low: ", sprintf("%.2f%%", 100 * match_fraction))
  )

  donor <- as.character(md[source_cells, donor_col, drop = TRUE])
  donor_missing <- blank(donor)
  assert(
    !any(matched & donor_missing),
    paste0(study, ": matched cells have missing donor values in ", donor_col)
  )

  cluster_lookup <- stats::setNames(
    as.character(cca_md[[config$cluster_column]]),
    rownames(cca_md)
  )
  cluster <- rep(NA_character_, length(source_cells))
  cluster[matched] <- unname(cluster_lookup[mapped_cca[matched]])

  use <- matched & !donor_missing & !blank(cluster)
  assert(any(use), paste0(study, ": no cells remain after CCA/donor filtering"))

  cells_dt <- data.table(
    source_cell = as.character(source_cells[use]),
    cca_cell = as.character(mapped_cca[use]),
    donor = as.character(donor[use]),
    cluster = as.character(cluster[use])
  )
  check_names(cells_dt$source_cell, paste0(study, " selected source cells"))
  check_names(cells_dt$cca_cell, paste0(study, " selected CCA cells"))

  groups <- unique(cells_dt[, .(donor, cluster)])
  groups[, key := paste(donor, cluster, sep = "\r")]
  cells_dt[, key := paste(donor, cluster, sep = "\r")]

  group_lookup <- stats::setNames(seq_len(nrow(groups)), groups$key)
  cells_dt[, group_index := as.integer(group_lookup[key])]
  assert(!anyNA(cells_dt$group_index), paste0(study, ": donor-cluster indexing failed"))

  group_n <- tabulate(cells_dt$group_index, nbins = nrow(groups))

  all_genes <- unique(as.character(unlist(layer_features, use.names = FALSE)))
  check_names(all_genes, paste0(study, " feature union"))

  gene_lookup <- stats::setNames(seq_along(all_genes), all_genes)
  cell_lookup <- stats::setNames(cells_dt$group_index, cells_dt$source_cell)

  pb <- Matrix::Matrix(
    0,
    nrow = length(all_genes),
    ncol = nrow(groups),
    sparse = TRUE,
    dimnames = list(all_genes, NULL)
  )

  expected_totals <- numeric(nrow(groups))
  represented_cells <- character()

  for (layer_name in count_layers) {
    message("  aggregating ", layer_name)
    mat <- SeuratObject::LayerData(object = rna, layer = layer_name)
    mat <- as(mat, "CsparseMatrix")

    layer_genes <- as.character(rownames(mat))
    assert(
      identical(layer_genes, getElement(layer_features, layer_name)),
      paste0(study, ": feature order changed between passes")
    )

    row_i <- as.integer(gene_lookup[layer_genes])
    assert(!anyNA(row_i), paste0(study, ": gene indexing failed"))

    row_map <- Matrix::sparseMatrix(
      i = row_i,
      j = seq_along(layer_genes),
      x = 1,
      dims = c(length(all_genes), length(layer_genes))
    )

    for (start in seq.int(1L, ncol(mat), by = config$block_cells)) {
      ix <- seq.int(start, min(ncol(mat), start + config$block_cells - 1L))
      block <- mat[, ix, drop = FALSE]
      check_integer_counts(block)

      block_cells <- as.character(colnames(block))
      gi <- as.integer(cell_lookup[block_cells])
      keep <- !is.na(gi)
      if (!any(keep)) next

      block <- block[, keep, drop = FALSE]
      gi <- gi[keep]
      block@x <- round(block@x)

      membership <- Matrix::sparseMatrix(
        i = seq_along(gi),
        j = gi,
        x = 1,
        dims = c(length(gi), nrow(groups))
      )

      pb <- pb + row_map %*% (block %*% membership)
      expected_totals <- expected_totals + as.numeric(
        Matrix::crossprod(membership, Matrix::colSums(block))
      )
      represented_cells <- c(represented_cells, as.character(colnames(block)))
    }

    rm(mat, row_map)
    invisible(gc())
  }

  check_names(represented_cells, paste0(study, " represented cells"))
  assert(
    base::setequal(represented_cells, cells_dt$source_cell),
    paste0(study, ": not every selected cell was represented exactly once")
  )

  pb_totals <- as.numeric(Matrix::colSums(pb))
  assert(all(pb_totals == expected_totals), paste0(study, ": pseudobulk total counts changed"))
  assert(all(pb_totals > 0), paste0(study, ": zero-library pseudobulk found"))

  pids <- paste0("S", safe_token(study), "_PB", seq_len(nrow(groups)))
  colnames(pb) <- pids

  donor_levels <- unique(as.character(groups$donor))
  donor_number <- as.integer(factor(as.character(groups$donor), levels = donor_levels))

  pb_meta <- data.table(
    pseudobulk_id = pids,
    study = study,
    donor = as.character(groups$donor),
    replicate = paste0("S", safe_token(study), "_D", donor_number),
    cluster = as.character(groups$cluster),
    n_cells = as.integer(group_n),
    donor_column = donor_col
  )

  audit <- data.table(
    condition = condition,
    study = study,
    status = "included",
    source_file = f,
    source_cells = length(source_cells),
    matched_cells = sum(matched),
    match_fraction = match_fraction,
    used_cells = sum(use),
    missing_donor = sum(matched & donor_missing),
    donors = uniqueN(cells_dt$donor),
    genes = nrow(pb),
    pseudobulks = ncol(pb),
    layers = paste(count_layers, collapse = ";"),
    donor_column = donor_col
  )

  save_csv(audit, file.path(out_dir, paste0("source_audit_", safe_token(study), ".csv")))

  rm(obj, md, rna)
  invisible(gc())

  list(counts = pb, metadata = pb_meta, cca_cells = cells_dt$cca_cell, audit = audit)
}

build_pseudobulk <- function(condition, cca_md, out_dir) {
  studies <- getElement(studies_by_condition, condition)
  parts <- list()
  audits <- list()
  seen <- character()

  for (study in studies) {
    f <- source_file_for(study, condition)

    if (study %in% excluded_nonraw) {
      message("Excluding ", study, " from DESeq2 raw-count pseudobulk")
      audits[study] <- list(
        data.table(
          condition = condition,
          study = study,
          status = "excluded",
          source_file = f,
          reason = "Current atlas matrix is processed/normalized rather than raw integer counts"
        )
      )
      next
    }

    p <- aggregate_study(study, condition, cca_md, out_dir)
    overlap <- base::intersect(seen, p$cca_cells)
    assert(!length(overlap), paste0(condition, ": CCA cells represented by multiple sources"))

    seen <- c(seen, p$cca_cells)
    parts[study] <- list(list(counts = p$counts, metadata = p$metadata))
    audits[study] <- list(p$audit)

    rm(p)
    invisible(gc())
  }

  assert(length(parts) > 0, paste0(condition, ": no eligible raw-count studies"))

  audit <- rbindlist(audits, use.names = TRUE, fill = TRUE)
  save_csv(audit, file.path(out_dir, "source_audit.csv"))

  clusters <- as.character(cca_md[[config$cluster_column]])
  coverage <- data.table(
    cell = rownames(cca_md),
    cluster = clusters,
    used = rownames(cca_md) %in% seen
  )[
    ,
    .(
      cca_cells = .N,
      used_cells = sum(used),
      coverage_fraction = sum(used) / .N
    ),
    by = cluster
  ]

  save_csv(coverage, file.path(out_dir, "CCA_cluster_coverage.csv"))

  pmeta <- rbindlist(lapply(parts, function(x) x$metadata), use.names = TRUE, fill = TRUE)
  save_csv(pmeta, file.path(out_dir, "pseudobulk_metadata.csv"))

  list(parts = parts, all_clusters = unique(clusters), audit = audit, coverage = coverage)
}

sort_cluster_ids <- function(x) {
  x <- unique(as.character(x))
  nx <- suppressWarnings(as.numeric(x))
  if (!anyNA(nx)) return(x[order(nx)])
  sort(x)
}

run_cluster <- function(pb, cluster_id, condition, out_dir) {
  d <- file.path(out_dir, paste0("cluster_", safe_token(cluster_id)))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(cluster_id), file.path(d, "cluster_id.txt"))

  md <- rbindlist(lapply(pb$parts, function(x) x$metadata), use.names = TRUE, fill = TRUE)

  rep_summary <- md[
    ,
    .(
      target_cells = sum(n_cells[cluster == cluster_id]),
      rest_cells = sum(n_cells[cluster != cluster_id])
    ),
    by = .(study, donor, replicate)
  ]

  rep_summary[
    ,
    eligible := target_cells >= config$min_cells_target &
      rest_cells >= config$min_cells_rest
  ]
  save_csv(rep_summary, file.path(d, "replicate_eligibility.csv"))

  eligible <- rep_summary[eligible == TRUE]

  summary_row <- data.table(
    condition = condition,
    method = "CCA",
    resolution = 0.5,
    cluster = cluster_id,
    paired_donors = nrow(eligible),
    studies = uniqueN(eligible$study),
    target_cells = sum(eligible$target_cells),
    rest_cells = sum(eligible$rest_cells),
    input_genes = NA_integer_,
    genes_tested = NA_integer_,
    genes_with_padj = NA_integer_,
    significant_up = 0L,
    significant_down = 0L,
    status = ""
  )

  skip <- function(reason) {
    summary_row[, status := reason]
    save_csv(summary_row, file.path(d, "cluster_summary.csv"))
    list(summary = summary_row, top = NULL)
  }

  if (nrow(eligible) < config$min_replicates) {
    return(skip("skipped: insufficient paired replicate units"))
  }

  participating_studies <- unique(as.character(eligible$study))
  gene_sets <- lapply(pb$parts[participating_studies], function(x) rownames(x$counts))
  genes <- Reduce(
    function(a, b) base::intersect(as.character(a), as.character(b)),
    gene_sets
  )

  if (!length(genes)) return(skip("skipped: no shared genes across participating studies"))
  writeLines(genes, file.path(d, "participating_common_genes.txt"))

  columns <- list()
  coldata <- list()

  for (i in seq_len(nrow(eligible))) {
    er <- eligible[i]
    study_id <- as.character(er$study)
    part <- getElement(pb$parts, study_id)
    pm <- part$metadata

    for (group_name in c("rest", "target")) {
      if (group_name == "target") {
        take <- pm$replicate == er$replicate & pm$cluster == cluster_id
      } else {
        take <- pm$replicate == er$replicate & pm$cluster != cluster_id
      }

      ids <- pm$pseudobulk_id[take]
      assert(length(ids) > 0, paste0(condition, " cluster ", cluster_id, ": empty ", group_name, " group"))

      id <- paste0("pair", i, "_", group_name)
      columns[id] <- list(Matrix::rowSums(part$counts[genes, ids, drop = FALSE]))
      coldata[id] <- list(
        data.frame(
          study = study_id,
          donor = as.character(er$donor),
          replicate = as.character(er$replicate),
          group = group_name,
          n_cells = if (group_name == "target") er$target_cells else er$rest_cells,
          row.names = id,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  x <- do.call(cbind, columns)
  rownames(x) <- genes
  cd <- do.call(rbind, unname(coldata))
  cd <- cd[colnames(x), , drop = FALSE]

  check_integer_counts(x)
  x <- round(x)
  assert(max(x) <= .Machine$integer.max, "Pseudobulk exceeds R integer range")

  summary_row[, input_genes := nrow(x)]
  keep <- rowSums(x) >= config$min_gene_total_count
  x <- x[keep, , drop = FALSE]

  if (!nrow(x)) return(skip("skipped: no genes pass total-count filter"))
  if (any(colSums(x) == 0)) return(skip("skipped: zero library after gene filtering"))

  storage.mode(x) <- "integer"
  cd$replicate <- factor(cd$replicate)
  cd$group <- factor(cd$group, levels = c("rest", "target"))

  design_matrix <- model.matrix(~ replicate + group, data = cd)
  assert(qr(design_matrix)$rank == ncol(design_matrix), "DESeq2 design is not full rank")
  assert(nrow(design_matrix) > ncol(design_matrix), "No residual degrees of freedom")

  save_csv(data.table(pseudobulk_id = rownames(cd), cd), file.path(d, "DESeq2_coldata.csv"))

  message("  cluster ", cluster_id, ": ", nrow(eligible), " paired replicate units; ", nrow(x), " genes")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = x,
    colData = cd,
    design = ~ replicate + group
  )

  dds <- DESeq2::DESeq(
    dds,
    sfType = "poscounts",
    minReplicatesForReplace = Inf,
    parallel = FALSE,
    quiet = TRUE
  )

  res <- DESeq2::results(
    dds,
    contrast = c("group", "target", "rest"),
    alpha = config$alpha
  )

  tab <- as.data.table(as.data.frame(res), keep.rownames = "gene")
  tab[
    ,
    `:=`(
      condition = condition,
      method = "CCA",
      resolution = 0.5,
      cluster = cluster_id,
      paired_donors = nrow(eligible),
      n_studies = length(participating_studies)
    )
  ]

  setcolorder(
    tab,
    c(
      "condition", "method", "resolution", "cluster", "paired_donors", "n_studies",
      "gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
    )
  )
  setorderv(tab, c("padj", "pvalue"), na.last = TRUE)

  up <- tab[
    !is.na(padj) & padj < config$alpha & log2FoldChange >= config$log2fc_cutoff
  ]
  down <- tab[
    !is.na(padj) & padj < config$alpha & log2FoldChange <= -config$log2fc_cutoff
  ]
  top <- head(up[order(padj, -log2FoldChange)], 50)
  if (nrow(top)) top[, rank := seq_len(.N)]

  save_csv(tab, file.path(d, "DESeq2_all_genes.csv"))
  save_csv(up, file.path(d, "DESeq2_upregulated.csv"))
  save_csv(down, file.path(d, "DESeq2_downregulated.csv"))
  save_csv(top, file.path(d, "DESeq2_top50_upregulated.csv"))
  save_csv(
    data.table(pseudobulk_id = colnames(dds), size_factor = DESeq2::sizeFactors(dds)),
    file.path(d, "size_factors.csv")
  )

  if (config$save_dds) saveRDS(dds, file.path(d, "DESeq2_dds.rds"), compress = FALSE)

  summary_row[
    ,
    `:=`(
      genes_tested = sum(!is.na(tab$pvalue)),
      genes_with_padj = sum(!is.na(tab$padj)),
      significant_up = nrow(up),
      significant_down = nrow(down),
      status = "completed"
    )
  ]
  save_csv(summary_row, file.path(d, "cluster_summary.csv"))

  rm(dds, res, tab, up, down, x, cd, design_matrix)
  invisible(gc())
  list(summary = summary_row, top = top)
}

run_condition <- function(condition, run_dir) {
  message("")
  message("Condition: ", condition)

  out_dir <- file.path(run_dir, condition)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  cca_md <- validate_cca(condition)
  cluster_values <- as.character(cca_md[[config$cluster_column]])

  message("  CCA cells: ", format(nrow(cca_md), big.mark = ","))
  message("  CCA clusters: ", length(unique(cluster_values)))

  save_csv(
    data.table(cell = rownames(cca_md), cluster = cluster_values),
    file.path(out_dir, "CCA_cluster_assignments.csv")
  )

  cache_file <- file.path(
    output_root,
    "validated_pseudobulk_cache",
    paste0(safe_token(condition), "_", config$cache_version, ".rds")
  )
  dir.create(dirname(cache_file), recursive = TRUE, showWarnings = FALSE)

  if (config$use_cache && file.exists(cache_file)) {
    message("  Loading pseudobulk cache: ", cache_file)
    pb <- readRDS(cache_file)
    save_csv(pb$audit, file.path(out_dir, "source_audit.csv"))
    save_csv(pb$coverage, file.path(out_dir, "CCA_cluster_coverage.csv"))
    pmeta <- rbindlist(lapply(pb$parts, function(x) x$metadata), use.names = TRUE, fill = TRUE)
    save_csv(pmeta, file.path(out_dir, "pseudobulk_metadata.csv"))
  } else {
    pb <- build_pseudobulk(condition, cca_md, out_dir)
    tmp <- tempfile(tmpdir = dirname(cache_file), fileext = ".rds")
    saveRDS(pb, tmp, compress = FALSE)
    assert(file.rename(tmp, cache_file), "Could not finalize pseudobulk cache")
  }

  summaries <- list()
  tops <- list()

  for (cluster_id in sort_cluster_ids(pb$all_clusters)) {
    message(condition, " CCA cluster ", cluster_id, " vs all other CCA clusters")

    result <- tryCatch(
      run_cluster(pb, cluster_id, condition, out_dir),
      error = function(e) {
        message("  Cluster failed: ", conditionMessage(e))
        failure <- data.table(
          condition = condition,
          method = "CCA",
          resolution = 0.5,
          cluster = cluster_id,
          paired_donors = NA_integer_,
          studies = NA_integer_,
          target_cells = NA_integer_,
          rest_cells = NA_integer_,
          input_genes = NA_integer_,
          genes_tested = NA_integer_,
          genes_with_padj = NA_integer_,
          significant_up = NA_integer_,
          significant_down = NA_integer_,
          status = paste0("ERROR: ", conditionMessage(e))
        )
        d <- file.path(out_dir, paste0("cluster_", safe_token(cluster_id)))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        save_csv(failure, file.path(d, "cluster_summary.csv"))
        list(summary = failure, top = NULL)
      }
    )

    summaries[cluster_id] <- list(result$summary)
    if (!is.null(result$top) && nrow(result$top)) tops[cluster_id] <- list(result$top)

    save_csv(
      rbindlist(summaries, use.names = TRUE, fill = TRUE),
      file.path(out_dir, "CCA_DESeq2_cluster_summary.csv")
    )
  }

  summary_table <- rbindlist(summaries, use.names = TRUE, fill = TRUE)
  save_csv(summary_table, file.path(out_dir, "CCA_DESeq2_cluster_summary.csv"))

  top_table <- NULL
  if (length(tops)) {
    top_table <- rbindlist(tops, use.names = TRUE, fill = TRUE)
    save_csv(top_table, file.path(out_dir, "CCA_DESeq2_top50_upregulated_all_clusters.csv"))
  }

  list(summary = summary_table, top = top_table)
}

run_workflow <- function() {
  run_dir <- file.path(output_root, paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

  writeLines(capture.output(sessionInfo()), file.path(run_dir, "sessionInfo.txt"))

  source_manifest <- rbindlist(
    lapply(conditions, function(condition) {
      rbindlist(lapply(getElement(studies_by_condition, condition), function(study) {
        data.table(
          condition = condition,
          study = study,
          source_file = source_file_for(study, condition),
          include_in_DESeq2 = !study %in% excluded_nonraw
        )
      }))
    })
  )
  save_csv(source_manifest, file.path(run_dir, "source_manifest.csv"))

  condition_manifest <- data.table(
    condition = conditions,
    cca_file = unname(cca_files[conditions]),
    expected_cells = as.integer(expected_cells[conditions]),
    cluster_column = config$cluster_column
  )
  save_csv(condition_manifest, file.path(run_dir, "condition_manifest.csv"))

  all_summaries <- list()
  all_tops <- list()
  statuses <- list()

  for (condition in conditions) {
    result <- tryCatch(
      run_condition(condition, run_dir),
      error = function(e) {
        message("Condition failed: ", condition, " - ", conditionMessage(e))
        d <- file.path(run_dir, condition)
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        writeLines(conditionMessage(e), file.path(d, "CONDITION_ERROR.txt"))
        list(error = conditionMessage(e))
      }
    )

    if (!is.null(result$error)) {
      statuses[condition] <- list(data.table(condition = condition, status = "failed", message = result$error))
      next
    }

    all_summaries[condition] <- list(result$summary)
    if (!is.null(result$top) && nrow(result$top)) all_tops[condition] <- list(result$top)
    statuses[condition] <- list(data.table(condition = condition, status = "completed", message = ""))
  }

  if (length(all_summaries)) {
    save_csv(
      rbindlist(all_summaries, use.names = TRUE, fill = TRUE),
      file.path(run_dir, "CCA_DESeq2_ALL_conditions_cluster_summary.csv")
    )
  }

  if (length(all_tops)) {
    save_csv(
      rbindlist(all_tops, use.names = TRUE, fill = TRUE),
      file.path(run_dir, "CCA_DESeq2_ALL_conditions_top50_upregulated_per_cluster.csv")
    )
  }

  save_csv(
    rbindlist(statuses, use.names = TRUE, fill = TRUE),
    file.path(run_dir, "condition_status.csv")
  )

  writeLines(normalizePath(run_dir, mustWork = TRUE), file.path(output_root, "LATEST_RUN.txt"))
  message("")
  message("CCA DESeq2 workflow finished")
  message("Run directory: ", run_dir)
  invisible(run_dir)
}

run_dir <- run_workflow()
