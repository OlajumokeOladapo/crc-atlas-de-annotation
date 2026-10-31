############################################################
# DIRECT LIMMA-VOOM REFINEMENT OF ORIGINAL UNRESOLVED CCA CLUSTERS
############################################################

suppressPackageStartupMessages({
    library(edgeR)
    library(limma)
    library(data.table)
    library(Matrix)
})

data.table::setDTthreads(1)

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset="/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

condition <- "TAN"
clusters_to_run <- c("17","26")

cache_dir <- file.path(
    atlas_root,
    "DESeq2_pseudobulk_CCA_res0.5_revised",
    "validated_pseudobulk_cache"
)

cache_candidates <- list.files(
    cache_dir,
    pattern=paste0("^",condition,"_.*\\.rds$"),
    full.names=TRUE
)

if (!length(cache_candidates)) {
    stop("No validated pseudobulk cache found for ",condition,
         ". Run the condition-specific 01_* pseudobulk script first.")
}

if (length(cache_candidates) > 1) {
    print(cache_candidates)
    stop("Multiple validated cache candidates found for ",condition,
         ". Set cache_file explicitly after confirming the correct cache.")
}

cache_file <- cache_candidates[1]

output_root <- file.path(
    atlas_root,
    paste0("limma_voom_",condition,"_unresolved_original_CCA")
)

dir.create(output_root,recursive=TRUE,showWarnings=FALSE)

settings <- list(
    min_cells_target=20,
    min_cells_rest=20,
    min_replicates=3,
    min_gene_total_count=10,
    min_count_per_sample=10,
    min_samples_expressed=3,
    fdr=0.05,
    log2fc=0.25,
    strong_log2fc=1
)

pb <- readRDS(cache_file)

if (!is.list(pb$parts) || !length(pb$parts)) {
    stop("Invalid pseudobulk cache: missing pb$parts.")
}

md <- rbindlist(
    lapply(
        names(pb$parts),
        function(study_id) {
            x <- as.data.table(pb$parts[[study_id]]$metadata)
            x[,study:=as.character(study)]
            x[,donor:=as.character(donor)]
            x[,replicate:=as.character(replicate)]
            x[,cluster:=as.character(cluster)]
            x
        }
    ),
    use.names=TRUE,
    fill=TRUE
)

all_summaries <- list()

run_one <- function(cluster_id) {

    cat("\n============================================================\n")
    cat(condition," ORIGINAL CCA CLUSTER ",cluster_id,
        " VS ALL OTHER CCA CLUSTERS\n",sep="")
    cat("============================================================\n")

    out <- file.path(output_root,paste0("cluster_",cluster_id))
    dir.create(out,recursive=TRUE,showWarnings=FALSE)

    eligible <- md[
        ,
        .(
            target_cells=sum(n_cells[cluster==cluster_id]),
            rest_cells=sum(n_cells[cluster!=cluster_id])
        ),
        by=.(study,donor,replicate)
    ]

    eligible[
        ,
        eligible := (
            target_cells >= settings$min_cells_target &
            rest_cells >= settings$min_cells_rest
        )
    ]

    fwrite(eligible,file.path(out,"replicate_eligibility.csv"))
    eligible <- eligible[eligible==TRUE]

    if (nrow(eligible) < settings$min_replicates) {
        summary <- data.table(
            condition=condition,
            cluster=cluster_id,
            method="limma_voomLmFit_eBayes",
            paired_replicates=nrow(eligible),
            n_studies=uniqueN(eligible$study),
            status="skipped_insufficient_paired_replicates"
        )
        fwrite(summary,file.path(out,"cluster_summary.csv"))
        return(summary)
    }

    studies <- unique(as.character(eligible$study))

    genes <- Reduce(
        intersect,
        lapply(pb$parts[studies],function(part) rownames(part$counts))
    )

    columns <- list()
    coldata <- list()

    for (i in seq_len(nrow(eligible))) {
        er <- eligible[i]
        study_id <- as.character(er$study)
        part <- pb$parts[[study_id]]

        pm <- as.data.table(part$metadata)
        pm[,replicate:=as.character(replicate)]
        pm[,cluster:=as.character(cluster)]

        rest_ids <- pm[
            replicate==as.character(er$replicate) & cluster!=cluster_id,
            pseudobulk_id
        ]
        target_ids <- pm[
            replicate==as.character(er$replicate) & cluster==cluster_id,
            pseudobulk_id
        ]

        if (!length(rest_ids) || !length(target_ids)) {
            stop("Empty target/rest pseudobulk for pair ",i)
        }

        rn <- paste0("pair",i,"_rest")
        tn <- paste0("pair",i,"_target")

        columns[[rn]] <- Matrix::rowSums(
            part$counts[genes,rest_ids,drop=FALSE]
        )
        columns[[tn]] <- Matrix::rowSums(
            part$counts[genes,target_ids,drop=FALSE]
        )

        coldata[[rn]] <- data.frame(
            pair=paste0("pair",i),
            study=study_id,
            donor=as.character(er$donor),
            source_replicate=as.character(er$replicate),
            group="rest",
            row.names=rn,
            stringsAsFactors=FALSE
        )
        coldata[[tn]] <- data.frame(
            pair=paste0("pair",i),
            study=study_id,
            donor=as.character(er$donor),
            source_replicate=as.character(er$replicate),
            group="target",
            row.names=tn,
            stringsAsFactors=FALSE
        )
    }

    counts <- do.call(cbind,columns)
    rownames(counts) <- genes

    cd <- do.call(rbind,unname(coldata))
    cd <- cd[colnames(counts),,drop=FALSE]

    counts <- round(counts)
    storage.mode(counts) <- "integer"

    keep <- (
        rowSums(counts) >= settings$min_gene_total_count &
        rowSums(counts >= settings$min_count_per_sample) >=
            settings$min_samples_expressed
    )
    counts <- counts[keep,,drop=FALSE]

    cd$pair <- factor(cd$pair)
    cd$group <- factor(cd$group,levels=c("rest","target"))
    design <- model.matrix(~pair+group,data=cd)

    if (
        qr(design)$rank != ncol(design) ||
        nrow(design) <= ncol(design)
    ) stop("Invalid paired design.")

    y <- edgeR::DGEList(counts)
    keep2 <- edgeR::filterByExpr(y,design=design)
    y <- y[keep2,,keep.lib.sizes=FALSE]
    y <- edgeR::calcNormFactors(y,method="TMM")

    fit <- edgeR::voomLmFit(
        counts=y,
        design=design,
        plot=FALSE,
        save.plot=TRUE
    )

    fit <- limma::eBayes(
        fit,
        trend=FALSE,
        robust=FALSE
    )

    tt <- limma::topTable(
        fit,
        coef="grouptarget",
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    results <- data.table(
        gene=rownames(tt),
        tt,
        check.names=FALSE
    )

    results[
        ,
        `:=`(
            condition=condition,
            cluster=cluster_id,
            comparison=paste0(
                "CCA_cluster_",cluster_id,
                "_VS_all_other_CCA_clusters"
            ),
            method="limma_voomLmFit_eBayes",
            FDR=adj.P.Val
        )
    ]

    valid <- is.finite(results$FDR) & is.finite(results$logFC)

    up <- results[
        valid & FDR < settings$fdr & logFC >= settings$log2fc
    ]
    down <- results[
        valid & FDR < settings$fdr & logFC <= -settings$log2fc
    ]
    strong <- results[
        valid & FDR < settings$fdr & logFC >= settings$strong_log2fc
    ]

    setorder(up,FDR,-logFC)
    setorder(down,FDR,logFC)
    setorder(strong,FDR,-logFC)

    outputs <- list(
        all_genes=results,
        upregulated=up,
        downregulated=down,
        strong_upregulated=strong,
        top100_upregulated=head(up,100),
        top100_strong_upregulated=head(strong,100)
    )

    for (nm in names(outputs)) {
        fwrite(
            outputs[[nm]],
            file.path(out,paste0("limma_voom_",nm,".csv"))
        )
    }

    summary <- data.table(
        condition=condition,
        cluster=cluster_id,
        method="limma_voomLmFit_eBayes",
        paired_replicates=nrow(eligible),
        n_studies=uniqueN(eligible$study),
        genes_tested=nrow(results),
        significant_up=nrow(up),
        significant_down=nrow(down),
        strong_up=nrow(strong),
        status="completed"
    )

    fwrite(summary,file.path(out,"cluster_summary.csv"))

    cat("\nTOP 25 UPREGULATED:\n")
    print(head(up[,.(gene,logFC,FDR)],25))

    summary
}

for (cluster_id in clusters_to_run) {
    all_summaries[[cluster_id]] <- tryCatch(
        run_one(cluster_id),
        error=function(e) data.table(
            condition=condition,
            cluster=cluster_id,
            method="limma_voomLmFit_eBayes",
            status="failed",
            error=conditionMessage(e)
        )
    )
}

summary <- rbindlist(all_summaries,use.names=TRUE,fill=TRUE)

fwrite(
    summary,
    file.path(
        output_root,
        paste0(condition,"_unresolved_CCA_limma_voom_summary.csv")
    )
)

capture.output(sessionInfo(),file=file.path(output_root,"sessionInfo.txt"))
print(summary)
