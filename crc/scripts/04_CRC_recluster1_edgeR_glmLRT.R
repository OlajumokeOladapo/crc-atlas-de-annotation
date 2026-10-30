############################################################
# CRC RECLUSTER 1: edgeR GLM-LRT FALLBACK
#
# This script is used for recluster 1 only.
# The primary DESeq2 model failed numerically for this cluster.
#
# Comparison:
#   recluster 1 vs all other new CRC reclusters
#
# Design:
#   ~ replicate + group
#
# Method:
#   edgeR DGEList
#   filterByExpr
#   TMM normalization
#   estimateDisp(robust = FALSE)
#   glmFit
#   glmLRT
#
# Positive logFC = higher in recluster 1.
############################################################


suppressPackageStartupMessages({
    library(edgeR)
    library(data.table)
    library(Matrix)
})


# ============================================================
# 1. SETTINGS
# ============================================================

cluster_id <- "1"

cache_file <- paste0(
    "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
    "DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5/",
    "validated_pseudobulk_cache/",
    "CRC_unresolved_reclusters_res0.5_raw_pseudobulk.rds"
)

output_root <- paste0(
    "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/",
    "edgeR_pseudobulk_CRC_unresolved_reclusters_res0.5_fallback"
)

cluster_dir <- file.path(
    output_root,
    paste0("cluster_", cluster_id)
)

dir.create(
    cluster_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

min_cells_target <- 20
min_cells_rest <- 20
min_replicates <- 3

min_gene_total_count <- 10
min_count_per_sample <- 10
min_samples_expressed <- 3

fdr_cutoff <- 0.05
log2fc_cutoff <- 0.25
strong_log2fc_cutoff <- 1


# ============================================================
# 2. LOAD CACHE
# ============================================================

pb <- readRDS(cache_file)

all_metadata <- rbindlist(
    lapply(
        names(pb$parts),
        function(study_id) {

            md <- as.data.table(
                pb$parts[[study_id]]$metadata
            )

            md[, study := as.character(study)]
            md[, donor := as.character(donor)]
            md[, replicate := as.character(replicate)]
            md[, cluster := as.character(cluster)]

            md
        }
    ),
    use.names = TRUE,
    fill = TRUE
)


# ============================================================
# 3. ELIGIBLE PAIRED REPLICATES
# ============================================================

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

eligible <- rep_summary[
    eligible == TRUE
]

if (nrow(eligible) < min_replicates) {
    stop("Insufficient paired biological replicates.")
}

fwrite(
    rep_summary,
    file.path(
        cluster_dir,
        "replicate_eligibility.csv"
    )
)


# ============================================================
# 4. COMMON GENES
# ============================================================

participating_studies <- unique(
    as.character(eligible$study)
)

gene_sets <- lapply(
    pb$parts[participating_studies],
    function(part) rownames(part$counts)
)

common_genes <- Reduce(
    intersect,
    gene_sets
)

common_genes <- as.character(common_genes)


# ============================================================
# 5. BUILD TARGET/REST MATRIX
# ============================================================

columns <- list()
coldata <- list()

for (i in seq_len(nrow(eligible))) {

    er <- eligible[i]
    study_id <- as.character(er$study)

    part <- pb$parts[[study_id]]

    pm <- as.data.table(part$metadata)

    pm[, replicate := as.character(replicate)]
    pm[, cluster := as.character(cluster)]

    rest_ids <- pm[
        replicate == as.character(er$replicate) &
        cluster != cluster_id,
        pseudobulk_id
    ]

    target_ids <- pm[
        replicate == as.character(er$replicate) &
        cluster == cluster_id,
        pseudobulk_id
    ]

    if (!length(rest_ids) || !length(target_ids)) {
        stop(
            paste0(
                "Missing target/rest pseudobulk IDs for pair ",
                i
            )
        )
    }

    rest_name <- paste0(
        "pair",
        i,
        "_rest"
    )

    target_name <- paste0(
        "pair",
        i,
        "_target"
    )

    columns[[rest_name]] <- Matrix::rowSums(
        part$counts[
            common_genes,
            rest_ids,
            drop = FALSE
        ]
    )

    columns[[target_name]] <- Matrix::rowSums(
        part$counts[
            common_genes,
            target_ids,
            drop = FALSE
        ]
    )

    coldata[[rest_name]] <- data.frame(
        study = study_id,
        donor = as.character(er$donor),
        replicate = paste0("pair", i),
        source_replicate = as.character(er$replicate),
        group = "rest",
        n_cells = er$rest_cells,
        row.names = rest_name,
        stringsAsFactors = FALSE
    )

    coldata[[target_name]] <- data.frame(
        study = study_id,
        donor = as.character(er$donor),
        replicate = paste0("pair", i),
        source_replicate = as.character(er$replicate),
        group = "target",
        n_cells = er$target_cells,
        row.names = target_name,
        stringsAsFactors = FALSE
    )
}

counts <- do.call(
    cbind,
    columns
)

rownames(counts) <- common_genes

cd <- do.call(
    rbind,
    unname(coldata)
)

cd <- cd[
    colnames(counts),
    ,
    drop = FALSE
]

counts <- round(counts)
storage.mode(counts) <- "integer"


# ============================================================
# 6. FILTER GENES
# ============================================================

keep_total <- rowSums(
    counts
) >= min_gene_total_count

keep_sparse <- rowSums(
    counts >= min_count_per_sample
) >= min_samples_expressed

counts_initial <- counts[
    keep_total & keep_sparse,
    ,
    drop = FALSE
]


# ============================================================
# 7. DESIGN
# ============================================================

cd$replicate <- factor(cd$replicate)

cd$group <- factor(
    cd$group,
    levels = c(
        "rest",
        "target"
    )
)

design <- model.matrix(
    ~ replicate + group,
    data = cd
)

if (
    qr(design)$rank != ncol(design) ||
    nrow(design) <= ncol(design)
) {
    stop("Invalid paired design.")
}

fwrite(
    data.table(
        pseudobulk_id = rownames(cd),
        cd
    ),
    file.path(
        cluster_dir,
        "edgeR_coldata.csv"
    )
)


# ============================================================
# 8. edgeR
# ============================================================

y <- DGEList(
    counts = counts_initial
)

keep <- filterByExpr(
    y,
    design = design
)

y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
]

y <- calcNormFactors(
    y,
    method = "TMM"
)

cat("Starting estimateDisp...\n")

y <- estimateDisp(
    y,
    design,
    robust = FALSE
)

cat("estimateDisp complete\n")

fit <- glmFit(
    y,
    design
)

cat("glmFit complete\n")

lrt <- glmLRT(
    fit,
    coef = "grouptarget"
)

cat("glmLRT complete\n")


# ============================================================
# 9. RESULTS
# ============================================================

tab <- as.data.table(
    topTags(
        lrt,
        n = Inf,
        sort.by = "PValue"
    )$table,
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
        n_studies = uniqueN(eligible$study),
        method = "edgeR_glmLRT"
    )
]

up <- tab[
    !is.na(FDR) &
    FDR < fdr_cutoff &
    logFC >= log2fc_cutoff
]

down <- tab[
    !is.na(FDR) &
    FDR < fdr_cutoff &
    logFC <= -log2fc_cutoff
]

strong_up <- tab[
    !is.na(FDR) &
    FDR < fdr_cutoff &
    logFC >= strong_log2fc_cutoff
]

setorder(
    up,
    FDR,
    -logFC
)

setorder(
    down,
    FDR,
    logFC
)

setorder(
    strong_up,
    FDR,
    -logFC
)

top100 <- head(
    up,
    100
)

strong_top100 <- head(
    strong_up,
    100
)


# ============================================================
# 10. SAVE
# ============================================================

fwrite(
    tab,
    file.path(
        cluster_dir,
        "edgeR_all_genes.csv"
    )
)

fwrite(
    up,
    file.path(
        cluster_dir,
        "edgeR_upregulated.csv"
    )
)

fwrite(
    down,
    file.path(
        cluster_dir,
        "edgeR_downregulated.csv"
    )
)

fwrite(
    strong_up,
    file.path(
        cluster_dir,
        "edgeR_strong_upregulated.csv"
    )
)

fwrite(
    top100,
    file.path(
        cluster_dir,
        "edgeR_top100_upregulated.csv"
    )
)

fwrite(
    strong_top100,
    file.path(
        cluster_dir,
        "edgeR_top100_strong_upregulated.csv"
    )
)

fwrite(
    data.table(
        pseudobulk_id = colnames(y),
        library_size = y$samples$lib.size,
        norm_factor = y$samples$norm.factors,
        effective_library_size = (
            y$samples$lib.size *
            y$samples$norm.factors
        )
    ),
    file.path(
        cluster_dir,
        "edgeR_TMM_normalization.csv"
    )
)

summary_row <- data.table(
    cluster = cluster_id,
    method = "edgeR_glmLRT",
    paired_replicates = nrow(eligible),
    n_studies = uniqueN(eligible$study),
    target_cells = sum(eligible$target_cells),
    rest_cells = sum(eligible$rest_cells),
    common_genes = length(common_genes),
    genes_after_initial_filter = nrow(counts_initial),
    genes_after_filterByExpr = nrow(y),
    genes_tested = nrow(tab),
    significant_up = nrow(up),
    significant_down = nrow(down),
    strong_up = nrow(strong_up),
    common_dispersion = y$common.dispersion,
    status = "completed"
)

fwrite(
    summary_row,
    file.path(
        cluster_dir,
        "edgeR_cluster_summary.csv"
    )
)

capture.output(
    sessionInfo(),
    file.path(
        cluster_dir,
        "sessionInfo.txt"
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
    nrow(strong_up),
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
                logFC,
                FDR
            )
        ],
        25
    )
)
