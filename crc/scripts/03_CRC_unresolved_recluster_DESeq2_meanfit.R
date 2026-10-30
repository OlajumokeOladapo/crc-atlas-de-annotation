############################################################
# CRC UNRESOLVED RECLUSTER PSEUDOBULK DESeq2
#
# Primary replicate-aware differential-expression workflow
# for CRC unresolved reclusters 0-18.
#
# Comparison:
#   each new recluster vs all other new reclusters
#
# Design:
#   ~ replicate + group
#
# Input:
#   existing validated raw-count pseudobulk cache
#
# Key settings:
#   min target cells = 20
#   min rest cells   = 20
#   min paired units = 3
#   FDR              = 0.05
#   |log2FC|         = 0.25
#   strong UP        = log2FC >= 1
#   DESeq2 fitType   = "mean"
############################################################


# ============================================================
# 1. PACKAGES
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
    library(data.table)
    library(Matrix)
})


# ============================================================
# 2. PATHS
# ============================================================

cache_file <- paste0(
    "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
    "DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5/",
    "validated_pseudobulk_cache/",
    "CRC_unresolved_reclusters_res0.5_raw_pseudobulk.rds"
)

output_root <- paste0(
    "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
    "DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5_meanfit"
)

dir.create(
    output_root,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# 3. SETTINGS
# ============================================================

clusters_to_run <- as.character(0:18)

min_cells_target <- 20
min_cells_rest <- 20
min_replicates <- 3

min_gene_total_count <- 10
min_count_per_sample <- 10
min_samples_expressed <- 3

alpha_cutoff <- 0.05
log2fc_cutoff <- 0.25
strong_log2fc_cutoff <- 1


# ============================================================
# 4. HELPERS
# ============================================================

assert <- function(condition, message) {
    if (!isTRUE(condition)) {
        stop(message, call. = FALSE)
    }
}

save_csv <- function(x, file) {
    fwrite(as.data.table(x), file)
}

safe_token <- function(x) {
    gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}


# ============================================================
# 5. LOAD VALIDATED RAW PSEUDOBULK CACHE
# ============================================================

assert(
    file.exists(cache_file),
    paste0("Cache file not found:\n", cache_file)
)

pb <- readRDS(cache_file)

assert(
    is.list(pb) && "parts" %in% names(pb) && length(pb$parts) > 0,
    "Invalid pseudobulk cache: expected non-empty pb$parts."
)

cat(
    "\nLoaded pseudobulk cache with ",
    length(pb$parts),
    " studies.\n",
    sep = ""
)


# ============================================================
# 6. COMBINE CACHE METADATA
# ============================================================

all_metadata <- rbindlist(
    lapply(
        names(pb$parts),
        function(study_id) {

            md <- as.data.table(
                pb$parts[[study_id]]$metadata
            )

            required <- c(
                "study",
                "donor",
                "replicate",
                "cluster",
                "pseudobulk_id",
                "n_cells"
            )

            assert(
                all(required %in% names(md)),
                paste0(
                    study_id,
                    ": missing required cache metadata."
                )
            )

            md[, study := as.character(study)]
            md[, donor := as.character(donor)]
            md[, replicate := as.character(replicate)]
            md[, cluster := as.character(cluster)]
            md[, pseudobulk_id := as.character(pseudobulk_id)]

            md
        }
    ),
    use.names = TRUE,
    fill = TRUE
)


# ============================================================
# 7. RESULT STORAGE
# ============================================================

summaries <- list()
tops <- list()
strong_tops <- list()


# ============================================================
# 8. RUN ONE RECLUSTER
# ============================================================

run_cluster <- function(cluster_id) {

    cat(
        "\n\n############################################################\n",
        "STARTING RECLUSTER ",
        cluster_id,
        "\n",
        "############################################################\n",
        sep = ""
    )

    d <- file.path(
        output_root,
        paste0("cluster_", safe_token(cluster_id))
    )

    dir.create(
        d,
        recursive = TRUE,
        showWarnings = FALSE
    )

    # --------------------------------------------------------
    # Eligible paired replicate units
    # --------------------------------------------------------

    rep_summary <- all_metadata[
        ,
        .(
            target_cells = sum(
                n_cells[cluster == cluster_id]
            ),
            rest_cells = sum(
                n_cells[cluster != cluster_id]
            )
        ),
        by = .(
            study,
            donor,
            replicate
        )
    ]

    rep_summary[
        ,
        eligible :=
            target_cells >= min_cells_target &
            rest_cells >= min_cells_rest
    ]

    save_csv(
        rep_summary,
        file.path(d, "replicate_eligibility.csv")
    )

    eligible <- rep_summary[
        eligible == TRUE
    ]

    summary_row <- data.table(
        cluster = cluster_id,
        comparison = paste0(
            "recluster_",
            cluster_id,
            "_VS_all_other_new_reclusters"
        ),
        recluster_cells = NA_integer_,
        raw_count_cells = NA_integer_,
        paired_replicates = nrow(eligible),
        n_studies = uniqueN(eligible$study),
        target_cells_in_model = sum(eligible$target_cells),
        rest_cells_in_model = sum(eligible$rest_cells),
        common_genes = NA_integer_,
        genes_after_total_filter = NA_integer_,
        genes_after_sparse_filter = NA_integer_,
        genes_tested = NA_integer_,
        genes_with_padj = NA_integer_,
        significant_up = 0L,
        significant_down = 0L,
        strong_up = 0L,
        status = ""
    )

    skip <- function(reason) {
        summary_row[, status := reason]
        save_csv(
            summary_row,
            file.path(d, "cluster_summary.csv")
        )
        return(
            list(
                summary = summary_row,
                top = NULL,
                strong_top = NULL
            )
        )
    }

    if (nrow(eligible) < min_replicates) {
        cat(
            "SKIPPED: insufficient paired replicate units: ",
            nrow(eligible),
            "\n",
            sep = ""
        )
        return(
            skip(
                paste0(
                    "insufficient paired replicate units: ",
                    nrow(eligible)
                )
            )
        )
    }

    # --------------------------------------------------------
    # Common genes across participating studies
    # --------------------------------------------------------

    participating_studies <- unique(
        as.character(eligible$study)
    )

    gene_sets <- lapply(
        pb$parts[participating_studies],
        function(x) rownames(x$counts)
    )

    genes <- Reduce(intersect, gene_sets)
    genes <- as.character(genes)

    if (!length(genes)) {
        return(
            skip(
                "skipped: no shared genes across participating studies"
            )
        )
    }

    summary_row[, common_genes := length(genes)]

    writeLines(
        genes,
        file.path(d, "participating_common_genes.txt")
    )

    # --------------------------------------------------------
    # Construct paired target/rest pseudobulk matrix
    # --------------------------------------------------------

    columns <- list()
    coldata <- list()

    for (i in seq_len(nrow(eligible))) {

        er <- eligible[i]
        study_id <- as.character(er$study)
        part <- pb$parts[[study_id]]

        pm <- as.data.table(part$metadata)
        pm[, replicate := as.character(replicate)]
        pm[, cluster := as.character(cluster)]

        for (group_name in c("rest", "target")) {

            if (group_name == "target") {
                take <- (
                    pm$replicate == as.character(er$replicate) &
                    pm$cluster == cluster_id
                )
            } else {
                take <- (
                    pm$replicate == as.character(er$replicate) &
                    pm$cluster != cluster_id
                )
            }

            ids <- pm$pseudobulk_id[take]

            assert(
                length(ids) > 0,
                paste0(
                    "Cluster ",
                    cluster_id,
                    ": empty ",
                    group_name,
                    " group for pair ",
                    i
                )
            )

            id <- paste0(
                "pair",
                i,
                "_",
                group_name
            )

            columns[[id]] <- Matrix::rowSums(
                part$counts[
                    genes,
                    ids,
                    drop = FALSE
                ]
            )

            coldata[[id]] <- data.frame(
                study = study_id,
                donor = as.character(er$donor),
                replicate = paste0("pair", i),
                source_replicate = as.character(er$replicate),
                group = group_name,
                n_cells = if (
                    group_name == "target"
                ) {
                    er$target_cells
                } else {
                    er$rest_cells
                },
                row.names = id,
                stringsAsFactors = FALSE
            )
        }
    }

    x <- do.call(cbind, columns)
    rownames(x) <- genes

    cd <- do.call(
        rbind,
        unname(coldata)
    )

    cd <- cd[
        colnames(x),
        ,
        drop = FALSE
    ]

    assert(
        all(is.finite(x)) &&
        all(x >= 0),
        paste0(
            "Cluster ",
            cluster_id,
            ": invalid pseudobulk counts."
        )
    )

    x <- round(x)
    storage.mode(x) <- "integer"

    # --------------------------------------------------------
    # Gene filtering
    # --------------------------------------------------------

    keep_total <- rowSums(x) >= min_gene_total_count

    summary_row[
        ,
        genes_after_total_filter :=
            sum(keep_total)
    ]

    x <- x[
        keep_total,
        ,
        drop = FALSE
    ]

    if (!nrow(x)) {
        return(
            skip(
                "skipped: no genes pass total-count filter"
            )
        )
    }

    keep_sparse <- rowSums(
        x >= min_count_per_sample
    ) >= min_samples_expressed

    x <- x[
        keep_sparse,
        ,
        drop = FALSE
    ]

    summary_row[
        ,
        genes_after_sparse_filter :=
            nrow(x)
    ]

    if (!nrow(x)) {
        return(
            skip(
                "skipped: no genes pass sparse-expression filter"
            )
        )
    }

    if (any(colSums(x) == 0)) {
        return(
            skip(
                "skipped: zero library after gene filtering"
            )
        )
    }

    # --------------------------------------------------------
    # Paired design
    # --------------------------------------------------------

    cd$replicate <- factor(cd$replicate)

    cd$group <- factor(
        cd$group,
        levels = c(
            "rest",
            "target"
        )
    )

    design_matrix <- model.matrix(
        ~ replicate + group,
        data = cd
    )

    assert(
        qr(design_matrix)$rank ==
        ncol(design_matrix),
        paste0(
            "Cluster ",
            cluster_id,
            ": DESeq2 design is not full rank."
        )
    )

    assert(
        nrow(design_matrix) >
        ncol(design_matrix),
        paste0(
            "Cluster ",
            cluster_id,
            ": no residual degrees of freedom."
        )
    )

    save_csv(
        data.table(
            pseudobulk_id = rownames(cd),
            cd
        ),
        file.path(
            d,
            "DESeq2_coldata.csv"
        )
    )

    cat(
        "\n============================================================\n",
        "NEW CRC RECLUSTER ",
        cluster_id,
        " VS ALL OTHER NEW RECLUSTERS\n",
        "============================================================\n",
        sep = ""
    )

    cat(
        "Eligible paired replicate units: ",
        nrow(eligible),
        "\n",
        sep = ""
    )

    cat(
        "Studies: ",
        length(participating_studies),
        "\n",
        sep = ""
    )

    cat(
        "Target cells: ",
        sum(eligible$target_cells),
        "\n",
        sep = ""
    )

    cat(
        "Rest cells: ",
        sum(eligible$rest_cells),
        "\n",
        sep = ""
    )

    cat(
        "Common genes: ",
        length(genes),
        "\n",
        sep = ""
    )

    cat(
        "Genes after total-count filter: ",
        summary_row$genes_after_total_filter,
        "\n",
        sep = ""
    )

    cat(
        "Genes after sparse-expression filter: ",
        nrow(x),
        "\n",
        sep = ""
    )

    cat(
        "Design: ",
        nrow(design_matrix),
        " x ",
        ncol(design_matrix),
        "\n",
        sep = ""
    )

    cat(
        "Design rank: ",
        qr(design_matrix)$rank,
        " / ",
        ncol(design_matrix),
        "\n",
        sep = ""
    )

    cat(
        "Running DESeq2 with fitType='mean'...\n"
    )

    # --------------------------------------------------------
    # DESeq2
    # --------------------------------------------------------

    dds <- DESeqDataSetFromMatrix(
        countData = x,
        colData = cd,
        design = ~ replicate + group
    )

    dds <- DESeq(
        dds,
        sfType = "poscounts",
        fitType = "mean",
        minReplicatesForReplace = Inf,
        parallel = FALSE,
        quiet = TRUE
    )

    res <- results(
        dds,
        contrast = c(
            "group",
            "target",
            "rest"
        ),
        alpha = alpha_cutoff
    )

    tab <- as.data.table(
        as.data.frame(res),
        keep.rownames = "gene"
    )

    tab[
        ,
        `:=`(
            cluster = cluster_id,
            comparison = paste0(
                "recluster_",
                cluster_id,
                "_VS_all_other_new_reclusters"
            ),
            paired_replicates = nrow(eligible),
            n_studies = length(participating_studies)
        )
    ]

    setcolorder(
        tab,
        c(
            "cluster",
            "comparison",
            "paired_replicates",
            "n_studies",
            "gene",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    )

    setorderv(
        tab,
        c(
            "padj",
            "pvalue"
        ),
        na.last = TRUE
    )

    up <- tab[
        !is.na(padj) &
        padj < alpha_cutoff &
        log2FoldChange >= log2fc_cutoff
    ]

    down <- tab[
        !is.na(padj) &
        padj < alpha_cutoff &
        log2FoldChange <= -log2fc_cutoff
    ]

    strong_up_tbl <- tab[
        !is.na(padj) &
        padj < alpha_cutoff &
        log2FoldChange >= strong_log2fc_cutoff
    ]

    setorder(
        up,
        padj,
        -log2FoldChange
    )

    setorder(
        down,
        padj,
        log2FoldChange
    )

    setorder(
        strong_up_tbl,
        padj,
        -log2FoldChange
    )

    top100 <- head(up, 100)
    strong_top100 <- head(
        strong_up_tbl,
        100
    )

    # --------------------------------------------------------
    # Save outputs
    # --------------------------------------------------------

    save_csv(
        tab,
        file.path(
            d,
            "DESeq2_all_genes.csv"
        )
    )

    save_csv(
        up,
        file.path(
            d,
            "DESeq2_upregulated.csv"
        )
    )

    save_csv(
        down,
        file.path(
            d,
            "DESeq2_downregulated.csv"
        )
    )

    save_csv(
        strong_up_tbl,
        file.path(
            d,
            "DESeq2_strong_upregulated.csv"
        )
    )

    save_csv(
        top100,
        file.path(
            d,
            "DESeq2_top100_upregulated.csv"
        )
    )

    save_csv(
        strong_top100,
        file.path(
            d,
            "DESeq2_top100_strong_upregulated.csv"
        )
    )

    save_csv(
        data.table(
            pseudobulk_id = colnames(dds),
            size_factor = sizeFactors(dds)
        ),
        file.path(
            d,
            "size_factors.csv"
        )
    )

    summary_row[
        ,
        `:=`(
            genes_tested = sum(!is.na(tab$pvalue)),
            genes_with_padj = sum(!is.na(tab$padj)),
            significant_up = nrow(up),
            significant_down = nrow(down),
            strong_up = nrow(strong_up_tbl),
            status = "completed"
        )
    ]

    save_csv(
        summary_row,
        file.path(
            d,
            "cluster_summary.csv"
        )
    )

    cat(
        "\nSignificant UP: ",
        nrow(up),
        "\n",
        sep = ""
    )

    cat(
        "Significant DOWN: ",
        nrow(down),
        "\n",
        sep = ""
    )

    cat(
        "Strong UP: ",
        nrow(strong_up_tbl),
        "\n",
        sep = ""
    )

    cat(
        "\nTOP 25 UPREGULATED GENES:\n"
    )

    print(
        head(
            up[
                ,
                .(
                    gene,
                    log2FoldChange,
                    padj
                )
            ],
            25
        )
    )

    rm(
        dds,
        res,
        tab,
        up,
        down,
        strong_up_tbl,
        x,
        cd,
        design_matrix,
        columns,
        coldata
    )

    invisible(gc())

    list(
        summary = summary_row,
        top = top100,
        strong_top = strong_top100
    )
}


# ============================================================
# 9. RUN ALL CRC NEW RECLUSTERS
# ============================================================

for (cluster_id in clusters_to_run) {

    out <- tryCatch(
        run_cluster(cluster_id),

        error = function(e) {

            cat(
                "\n============================================================\n",
                "ERROR — RECLUSTER ",
                cluster_id,
                "\n",
                conditionMessage(e),
                "\n",
                "============================================================\n",
                sep = ""
            )

            d <- file.path(
                output_root,
                paste0(
                    "cluster_",
                    safe_token(cluster_id)
                )
            )

            dir.create(
                d,
                recursive = TRUE,
                showWarnings = FALSE
            )

            error_row <- data.table(
                cluster = cluster_id,
                comparison = paste0(
                    "recluster_",
                    cluster_id,
                    "_VS_all_other_new_reclusters"
                ),
                status = paste0(
                    "ERROR: ",
                    conditionMessage(e)
                )
            )

            save_csv(
                error_row,
                file.path(
                    d,
                    "cluster_summary.csv"
                )
            )

            list(
                summary = error_row,
                top = NULL,
                strong_top = NULL
            )
        }
    )

    summaries[[cluster_id]] <- out$summary

    if (!is.null(out$top)) {
        tops[[cluster_id]] <- out$top
    }

    if (!is.null(out$strong_top)) {
        strong_tops[[cluster_id]] <- out$strong_top
    }

    invisible(gc())
}


# ============================================================
# 10. FINAL SUMMARY
# ============================================================

summary_table <- rbindlist(
    summaries,
    use.names = TRUE,
    fill = TRUE
)

summary_table[
    ,
    cluster_numeric :=
        suppressWarnings(
            as.numeric(cluster)
        )
]

setorder(
    summary_table,
    cluster_numeric
)

summary_table[
    ,
    cluster_numeric := NULL
]

save_csv(
    summary_table,
    file.path(
        output_root,
        "CRC_unresolved_recluster_DESeq2_summary.csv"
    )
)

cat(
    "\n\n============================================================\n",
    "FINAL CORRECTED DESeq2 SUMMARY\n",
    "============================================================\n",
    sep = ""
)

print(
    summary_table,
    nrows = Inf
)


# ============================================================
# 11. COMBINED TOP 100
# ============================================================

if (length(tops) > 0) {

    top_table <- rbindlist(
        tops,
        use.names = TRUE,
        fill = TRUE
    )

    top_table[
        ,
        cluster_numeric :=
            as.numeric(cluster)
    ]

    setorder(
        top_table,
        cluster_numeric,
        padj,
        -log2FoldChange
    )

    top_table[
        ,
        cluster_numeric := NULL
    ]

    save_csv(
        top_table,
        file.path(
            output_root,
            "CRC_unresolved_recluster_DESeq2_top100_all_clusters.csv"
        )
    )
}


# ============================================================
# 12. COMBINED STRONG TOP 100
# ============================================================

if (length(strong_tops) > 0) {

    strong_table <- rbindlist(
        strong_tops,
        use.names = TRUE,
        fill = TRUE
    )

    strong_table[
        ,
        cluster_numeric :=
            as.numeric(cluster)
    ]

    setorder(
        strong_table,
        cluster_numeric,
        padj,
        -log2FoldChange
    )

    strong_table[
        ,
        cluster_numeric := NULL
    ]

    save_csv(
        strong_table,
        file.path(
            output_root,
            "CRC_unresolved_recluster_DESeq2_top100_strong_all_clusters.csv"
        )
    )
}


# ============================================================
# 13. SESSION INFO
# ============================================================

capture.output(
    sessionInfo(),
    file = file.path(
        output_root,
        "sessionInfo.txt"
    )
)


# ============================================================
# 14. DONE
# ============================================================

cat(
    "\n============================================================\n",
    "CORRECTED CRC RECLUSTER DESeq2 COMPLETE\n",
    "============================================================\n",
    sep = ""
)

cat(
    "\nResults:\n",
    output_root,
    "\n",
    sep = ""
)

cat(
    "\nExisting pseudobulk cache reused:\n",
    cache_file,
    "\n",
    sep = ""
)
