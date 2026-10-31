# Complete OSCER R GUI workflow: saved subclusters of parents 9 and 12.
# If needed, set verified_map_file below to the EXISTING historical donor-recovery CSV.
# Missing historical multiplexing mappings cannot be inferred from pooled GEO titles.
# Reads only whitelisted function definitions from the old script; never sources it.
required <- c("SeuratObject", "Matrix", "data.table", "edgeR", "statmod")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse=", "))
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1)
ql_method <- getS3method("glmQLFit", "default", envir=asNamespace("edgeR"))
if (!"legacy" %in% names(formals(ql_method))) stop("This script requires edgeR's legacy=FALSE QL implementation.")

atlas <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
recluster_dir <- file.path(atlas, "CRC_individual_reclustering", "run_20261030_155306_1a9ab35069091b")
mapping_script <- file.path(atlas, "crc/scripts/01_CRC_original_CCA_pseudobulk_DESeq2.R")
stopifnot(file.exists(mapping_script), dir.exists(recluster_dir))
out <- tempfile(pattern="edgeR_verified_cells_only_", tmpdir=recluster_dir)
dir.create(out)
cat("Output directory: ", out, "\n", sep="")

# Import function DEFINITIONS only into an isolated environment.
helpers <- new.env(parent=globalenv())
wanted <- c("assert", "blank", "safe_token", "save_csv", "check_names",
            "check_integer_counts", "source_file_for", "find_donor_column", "aggregate_study")
expressions <- parse(mapping_script)
loaded <- character()
for (e in expressions) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && length(e)==3L &&
        is.symbol(e[[2]]) && as.character(e[[2]]) %in% wanted &&
        is.call(e[[3]]) && identical(e[[3]][[1]], as.name("function"))) {
        eval(e, envir=helpers)
        loaded <- c(loaded, as.character(e[[2]]))
    }
}
if (!setequal(loaded, wanted)) stop("Required helper definitions missing from mapping script.")
helpers$source_dir <- file.path(atlas, "full_seurat_QC_normalized_by_study_condition")
# Optional EXISTING cell-level donor recovery table. Leave NULL if saved
# verified fields or the explicit rules below cover all selected cells.
# Required CSV columns: study, cell, verified_donor_id,
# donor_source_field, donor_source_value.
verified_map_file <- file.path(recluster_dir,
    "GSE225857_donor_recovery_1a9ab375aa375c",
    "GSE225857_verified_donor_map.csv")
# If restoration was rerun, replace this with the newly printed CSV path.
verified_map <- NULL
if (!is.null(verified_map_file)) {
    verified_map <- fread(verified_map_file, colClasses="character")
    required_map <- c("study","cell","verified_donor_id",
                      "donor_source_field","donor_source_value")
    if (!all(required_map %in% names(verified_map))) stop("Donor map lacks required columns.")
    if (anyDuplicated(verified_map[,.(study,cell)])) stop("Duplicate study/cell donor mappings.")
    if (anyNA(verified_map[,..required_map]) ||
        any(vapply(verified_map[,..required_map],function(z) any(!nzchar(trimws(z))),logical(1)))) {
        stop("Missing IDs or provenance in donor map.")
    }
}
raw_studies <- c("GSE132465","GSE178341","GSE200997","GSE201348",
                 "GSE221575","GSE225857","GSE231559","GSE261388",
                 "GSE294300","GSE302903","GSE330797")
helpers$donor_overrides <- setNames(rep("verified_donor_id",length(raw_studies)),raw_studies)
helpers$config <- list(cluster_column="new_cluster", min_match_fraction=0,
                       block_cells=2000L, integer_tolerance=1e-8)
helpers$blank <- function(x) {
    z <- trimws(as.character(x))
    is.na(z) | tolower(z) %in% c("", "na", "n/a", "nan", "null", "unknown", "not reported", "not available")
}
# No generic column selection or sample/library fallback.
helpers$find_donor_column <- function(md, study) {
    if (!study %in% names(helpers$donor_overrides)) stop("No donor override for ",study)
    column <- helpers$donor_overrides[[study]]
    if (!column %in% names(md)) stop(study,": verified donor field missing.")
    column
}

recover_donors <- function(md, study, chosen_cells, external=NULL) {
    cells <- rownames(md)
    if (is.null(cells) || anyDuplicated(cells)) stop(study,": invalid cell names.")
    ix <- match(chosen_cells,cells)
    if (anyNA(ix)) stop(study,": selected cell IDs absent from raw-source metadata.")
    m <- md[ix,,drop=FALSE]
    getcol <- function(nm) {
        if (nm %in% names(m)) as.character(m[[nm]]) else rep(NA_character_,nrow(m))
    }
    blank <- helpers$blank
    audit <- data.table(study=study,cell=chosen_cells,
                        verified_donor_id=getcol("verified_donor_id"),
                        donor_source_field=getcol("donor_source_field"),
                        donor_source_value=getcol("donor_source_value"),
                        donor_rule=rep("previous_verified_field",nrow(m)))
    # Preserve existing verified identities. If original provenance is absent,
    # say so explicitly rather than inventing it.
    prior <- !blank(audit$verified_donor_id)
    audit[prior & blank(donor_source_field),donor_source_field:="existing_verified_donor_id"]
    audit[prior & blank(donor_source_value),donor_source_value:=verified_donor_id]
    audit[prior & donor_source_field=="existing_verified_donor_id",
          donor_rule:="previous_verified_field_original_provenance_unavailable"]
    audit[!prior,verified_donor_id:=NA_character_]
    assign_values <- function(ids, values, field, rule, check_conflict=FALSE) {
        ids <- as.character(ids)
        ok <- !blank(ids)
        conflict <- ok & !blank(audit$verified_donor_id) & audit$verified_donor_id!=ids
        if (check_conflict && any(conflict,na.rm=TRUE)) stop(study,": conflicting verified donor mappings.")
        use <- ok & blank(audit$verified_donor_id)
        audit[use, `:=`(verified_donor_id=ids[use], donor_source_field=field,
                        donor_source_value=as.character(values)[use],donor_rule=rule)]
    }
    # An audited external mapping takes precedence over new parsing, but cannot
    # contradict a saved verified donor.
    if (!is.null(external)) {
        ext <- external[external$study==study,]
        k <- match(chosen_cells,ext$cell)
        ids <- ext$verified_donor_id[k]
        conflict <- !blank(ids) & !blank(audit$verified_donor_id) & ids!=audit$verified_donor_id
        if (any(conflict,na.rm=TRUE)) stop(study,": saved and external donor IDs disagree.")
        use <- !blank(ids)
        audit[use, `:=`(verified_donor_id=ids[use],
                        donor_source_field=ext$donor_source_field[k[use]],
                        donor_source_value=ext$donor_source_value[k[use]],
                        donor_rule="external_verified_cell_map")]
    }
    capture <- function(x, pattern, replacement="\\1") {
        answer <- rep(NA_character_,length(x))
        ok <- !is.na(x) & grepl(pattern,x,perl=TRUE)
        answer[ok] <- sub(pattern,replacement,x[ok],perl=TRUE)
        answer
    }
    # The supplied preprocessing workflow saves GEO patients under this name.
    # Preserve that restored field rather than falling back to titles/library IDs.
    if (study=="GSE225857") {
        v <- getcol("gse225857_patients")
        assign_values(v,v,"gse225857_patients",
                      "GEO_patients_restored_by_preprocessing",check_conflict=TRUE)
    }
    # First use the explicit source fields specified for these studies.
    if (study=="GSE132465") {
        v <- getcol("geo_patient_id.ch1")
        assign_values(capture(v,"^(SMC[0-9]+)$"),v,"geo_patient_id.ch1","patient_metadata")
    }
    if (study=="GSE178341") {
        v <- getcol("geo_patienttypeid.ch1")
        assign_values(capture(v,"^(C[0-9]+)(?:[^0-9].*)?$"),v,
                      "geo_patienttypeid.ch1","leading_patient_code")
    }
    title <- trimws(getcol("geo_title"))
    switch(study,
        GSE132465=assign_values(capture(chosen_cells,"^(SMC[0-9]{2})-T_.+$"),
                               chosen_cells,"cell_id","SMC_patient_prefix"),
        GSE178341=assign_values(capture(chosen_cells,"^(C[0-9]+)_.+$"),
                               chosen_cells,"cell_id","C_patient_prefix"),
        GSE200997=assign_values(capture(chosen_cells,"^T_(cac[0-9]+)_.+$"),
                               chosen_cells,"cell_id","cac_patient_prefix"),
        GSE201348=assign_values(capture(title,"^([A-Za-z0-9_-]+),.*$"),
                               title,"geo_title","user_confirmed_subject_code_before_comma"),
        # GSE221575 title numbers do not establish donor identity.
        # Require saved verified IDs or the original cleaned donor mapping.
        GSE221575=NULL,
        # GSE225857 immune/nonimmune titles do NOT identify patients.
        # Only saved verified IDs or the external recovery map can resolve it.
        GSE225857=NULL,
        GSE231559=assign_values(capture(title,"^([CL][0-9]+)[NT]$"),
                               title,"geo_title","user_confirmed_code_before_terminal_N_T"),
        GSE261388=assign_values(capture(chosen_cells,"^GSE261388_(T[0-9]+)_.+$"),
                               chosen_cells,"cell_id","user_confirmed_source_subject_token"),
        GSE294300=assign_values(capture(title,"^.*\\b(patient_[0-9]+)\\b.*$"),
                               title,"geo_title","explicit_patient_token"),
        GSE302903=assign_values(capture(title,"^(T[0-9]+)$"),
                               title,"geo_title","user_confirmed_source_patient_token"),
        GSE330797=assign_values(capture(title,"^.*\\bpatient ([0-9]+)\\b.*$"),
                               title,"geo_title","explicit_patient_number")
    )
    audit[, `:=`(source_title=title,source_sample=getcol("geo_sample_id"),
                  donor_valid=!blank(verified_donor_id))]
    audit
}

# This is the readRDS used only by the imported aggregation helper.
# It updates the object in memory, never overwrites the source RDS.
helpers$readRDS <- function(file, ...) {
    object <- base::readRDS(file,...)
    study <- sub("_CRC_QC_normalized_seurat\\.rds$","",basename(file))
    if (!study %in% raw_studies) return(object)
    selected <- assignments$cell[assignments$study==study]
    audit <- recover_donors(object[[]],study,selected,verified_map)
    fwrite(audit,file.path(out,paste0(study,"_verified_donor_audit.csv")))
    if (any(!audit$donor_valid)) {
        missing <- audit[donor_valid==FALSE]
        fwrite(missing,file.path(out,paste0(study,"_unresolved_donor_cells.csv")))
        stop(study,": ",nrow(missing)," selected cells lack verified donor mapping. See donor audit.")
    }
    m <- object[[]]
    for (field in c("verified_donor_id","donor_source_field","donor_source_value")) {
        v <- if(field %in% names(m)) as.character(m[[field]]) else rep(NA_character_,nrow(m))
        v[match(audit$cell,rownames(m))] <- audit[[field]]
        object[[field]] <- v
    }
    object
}

assignments <- rbindlist(lapply(c("9", "12"), function(parent) {
    f <- file.path(recluster_dir, paste0("parent_cluster_", parent), "cell_subcluster_assignments.csv")
    a <- fread(f, colClasses="character")
    stopifnot(all(c("cell", "parent_cluster", "subcluster", "study") %in% names(a)),
              all(a$parent_cluster==parent))
    a[, new_cluster := paste(parent_cluster, subcluster, sep=":")]
    a
}))
if (anyDuplicated(assignments$cell) || anyNA(assignments) ||
    any(!nzchar(assignments$cell)) || any(!nzchar(assignments$study))) stop("Invalid assignments.")

# Explicit final exclusion: ONLY the previously audited 117 unresolved cells.
original_assignments <- copy(assignments)
blocked_run <- file.path(recluster_dir,"edgeR_1a9ab34daec974")
expected_exclusions <- c(GSE221575=65L,GSE231559=52L)
exclusions <- rbindlist(lapply(names(expected_exclusions),function(s) {
    f <- file.path(blocked_run,paste0(s,"_unresolved_donor_cells.csv"))
    if(!file.exists(f)) stop("Missing exclusion audit: ",f)
    a <- fread(f,colClasses="character")
    if(!all(c("study","cell") %in% names(a))) stop("Invalid exclusion audit: ",f)
    if(anyNA(a$study) || anyNA(a$cell) || !all(a$study==s) ||
       nrow(a)!=expected_exclusions[[s]]) stop("Unexpected exclusions: ",s)
    a[,exclusion_reason:="unresolved_biological_donor"]
    a
}),fill=TRUE)
if(anyDuplicated(exclusions$cell)) stop("Duplicate excluded cells.")
k <- match(exclusions$cell,original_assignments$cell)
if(anyNA(k) || any(exclusions$study!=original_assignments$study[k]))
    stop("Exclusions do not match original assignments.")
exclusions[,`:=`(parent_cluster=original_assignments$parent_cluster[k],
                 subcluster=original_assignments$subcluster[k])]
assignments <- original_assignments[!cell %in% exclusions$cell]
stopifnot(nrow(original_assignments)-nrow(assignments)==117L)
fwrite(exclusions,file.path(out,"excluded_unresolved_donor_cells.csv"))
fwrite(exclusions[,.(excluded_cells=.N),by=.(study,parent_cluster,subcluster)],
       file.path(out,"excluded_unresolved_donor_summary.csv"))
# No additional exclusions based on anatomical site are made.

fwrite(assignments, file.path(out, "input_assignments.csv"))
excluded <- c("CELLxGENE_Chen_2021", "CELLxGENE_Moorman2025")
fwrite(assignments[study %in% excluded, .N, by=.(study,parent_cluster,subcluster)],
       file.path(out,"excluded_processed_sources.csv"))
studies <- setdiff(unique(assignments$study), excluded)
if(any(!studies %in% raw_studies)) stop("Unconfigured raw-count study.")
# Per-study caches survive interrupted runs. Cache reuse requires identical
# assignments, raw-file size/mtime, mapping-table hash and recovery rules.
cache_dir <- file.path(recluster_dir,"edgeR_verified_donor_pseudobulk_cache")
dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
map_hash <- if(is.null(verified_map_file)) "none" else unname(tools::md5sum(verified_map_file))
parts <- list()
source_status <- list()
for (s in studies) {
    cat("\nBUILDING RAW PSEUDOBULK: ",s,"\n",sep=""); flush.console()
    chosen <- assignments[study==s]
    f <- file.path(helpers$source_dir,paste0(s,"_CRC_QC_normalized_seurat.rds"))
    cached <- file.path(cache_dir,paste0(s,".rds"))
    signature <- list(version=if(s %in% c("GSE221575","GSE225857")) "verified_donors_v3" else "verified_donors_v2",assignment=as.data.frame(chosen),
                      source=f,size=file.info(f)$size,mtime=as.character(file.info(f)$mtime),
                      map_hash=if(is.null(verified_map) || !any(verified_map$study==s)) "none" else map_hash,helper_hash=unname(tools::md5sum(mapping_script)))
    result <- tryCatch({
        reuse <- FALSE
        if(file.exists(cached)) {
            z <- readRDS(cached)
            if(identical(z$signature,signature)) {
                parts[[s]] <- z$part
                fwrite(z$donor_audit,file.path(out,paste0(s,"_verified_donor_audit.csv")))
                fwrite(z$source_audit,file.path(out,paste0("source_audit_",s,".csv")))
                reuse <- TRUE
                cat("Reused verified pseudobulk cache.\n")
            }
            rm(z)
        }
        if(!reuse) {
            labels <- data.frame(new_cluster=chosen$new_cluster,row.names=chosen$cell)
            p <- helpers$aggregate_study(s,"CRC",labels,out)
            if(!setequal(p$cca_cells,chosen$cell)) {
                fwrite(chosen[!cell %in% p$cca_cells],file.path(out,paste0("unmatched_",s,".csv")))
                stop("Selected-cell barcode match incomplete for ",s)
            }
            p$metadata[,parent:=sub(":.*$","",cluster)]
            p$metadata[,subcluster:=sub("^[^:]+:","",cluster)]
            part <- list(counts=p$counts,metadata=p$metadata)
            record <- list(signature=signature,part=part,source_audit=p$audit,
                           donor_audit=fread(file.path(out,paste0(s,"_verified_donor_audit.csv")),colClasses="character"))
            temp <- tempfile(tmpdir=cache_dir,fileext=".rds")
            saveRDS(record,temp,compress=FALSE)
            if(!file.rename(temp,cached)) stop("Could not save verified cache for ",s)
            parts[[s]] <- part
            rm(p,part,record)
        }
        data.table(study=s,status=if(reuse) "cached" else "aggregated",error="")
    },error=function(e) data.table(study=s,status="blocked",error=conditionMessage(e)))
    source_status[[s]] <- result
    fwrite(rbindlist(source_status),file.path(out,"source_status.csv"))
    print(result)
    invisible(gc())
}
source_status <- rbindlist(source_status)
if(any(source_status$status=="blocked")) {
    stop("DE has not started: one or more sources need verified mapping or count fixes. ",
         "Other study caches were preserved. See ",file.path(out,"source_status.csv"))
}
if (!length(parts)) stop("No raw-count studies available.")
md <- rbindlist(lapply(parts, function(p) p$metadata))
fwrite(md, file.path(out,"pseudobulk_metadata.csv"))
saveRDS(list(parts=parts, assignments=assignments, source_script=mapping_script),
        file.path(out,"new_subcluster_pseudobulk_cache.rds"), compress=FALSE)

save_diagnostic <- function(draw, stem) {
    for (fmt in c("pdf","png")) {
        tryCatch(local({
            if (fmt=="pdf") grDevices::pdf(paste0(stem,".pdf"),width=8,height=6)
            else {
                if (!capabilities("cairo")) stop("Cairo unavailable; PNG skipped")
                grDevices::png(paste0(stem,".png"),width=8,height=6,units="in",res=300,type="cairo")
            }
            on.exit(grDevices::dev.off())
            draw()
        }),error=function(e) cat(fmt,conditionMessage(e),"\n",
                               file=paste0(stem,"_plot_errors.txt"),append=TRUE))
    }
}

fit_one <- function(parent_id, sub_id, folder) {
    m <- md[parent==parent_id]
    eligible <- m[, .(target_cells=sum(n_cells[subcluster==sub_id]),
                      rest_cells=sum(n_cells[subcluster!=sub_id])), by=.(study,donor)]
    eligible[, eligible := target_cells>=20 & rest_cells>=20]
    fwrite(eligible,file.path(folder,"replicate_eligibility.csv"))
    eligible <- eligible[eligible==TRUE]
    base <- data.frame(parent=parent_id,subcluster=sub_id,paired_donors=nrow(eligible),
                       studies=uniqueN(eligible$study),genes_tested=NA_integer_,
                       significant_up=NA_integer_,significant_down=NA_integer_,status="",error="")
    if (nrow(eligible)<3L) {
        base$status <- "skipped_insufficient_paired_donors"
        return(base)
    }
    included <- unique(eligible$study)
    genes <- Reduce(intersect,lapply(parts[included],function(p) rownames(p$counts)))
    if (!length(genes)) stop("No shared genes across participating studies.")
    writeLines(genes,file.path(folder,"tested_gene_universe_before_filtering.txt"))
    cd <- as.data.frame(eligible[rep(seq_len(.N),each=2L)])
    cd$pair <- factor(rep(seq_len(nrow(eligible)),each=2L))
    cd$group <- factor(rep(c("rest","target"),nrow(eligible)),levels=c("rest","target"))
    rownames(cd) <- paste0("pair",cd$pair,"_",cd$group)
    counts <- matrix(0,length(genes),nrow(cd),dimnames=list(genes,rownames(cd)))
    for (j in seq_len(nrow(cd))) {
        p <- parts[[cd$study[j]]]
        pm <- p$metadata
        take <- pm$parent==parent_id & pm$donor==cd$donor[j] &
                ((pm$subcluster==sub_id)==(cd$group[j]=="target"))
        ids <- pm$pseudobulk_id[take]
        if (!length(ids)) stop("Empty target/rest group.")
        counts[,j] <- Matrix::rowSums(p$counts[genes,ids,drop=FALSE])
    }
    design <- model.matrix(~pair+group,cd)
    if (qr(design)$rank!=ncol(design) || nrow(design)<=ncol(design)) stop("Invalid paired design.")
    if (any(colSums(counts)<=0)) stop("Zero library.")
    y <- edgeR::DGEList(counts)
    keep <- edgeR::filterByExpr(y,design=design) & rowSums(counts)>=10
    y <- y[keep,,keep.lib.sizes=FALSE]
    if (!nrow(y) || any(y$samples$lib.size<=0)) stop("No genes or empty library after filtering.")
    y <- edgeR::calcNormFactors(y,method="TMM")
    rm(counts); invisible(gc())
    fwrite(data.table(library=rownames(cd),cd),file.path(folder,"model_metadata.csv"))
    fwrite(data.table(library=rownames(y$samples),y$samples),file.path(folder,"library_normalization.csv"))
    saveRDS(list(y=y,design=design,metadata=cd),file.path(folder,"normalized_checkpoint.rds"),compress=FALSE)
    cat("  edgeR QL: ",nrow(y)," genes, ",ncol(y)," libraries\n",sep=""); flush.console()
    fit <- edgeR::glmQLFit(y,design=design,legacy=FALSE,robust=TRUE)
    saveRDS(fit,file.path(folder,"edgeR_QL_fit.rds"),compress=FALSE)
    test <- edgeR::glmQLFTest(fit,coef="grouptarget")
    tab <- as.data.table(edgeR::topTags(test,n=Inf,sort.by="PValue")$table,keep.rownames="gene")
    tab[, `:=`(parent=parent_id,subcluster=sub_id,method="edgeR_QLF",
               comparison=paste0("parent",parent_id,"_sub",sub_id,"_vs_other_subclusters_within_parent"))]
    good <- is.finite(tab$FDR) & is.finite(tab$logFC)
    up <- tab[good & FDR<0.05 & logFC>=0.25]
    down <- tab[good & FDR<0.05 & logFC<= -0.25]
    strong <- tab[good & FDR<0.05 & logFC>=1]
    outputs <- list(all_genes=tab,upregulated=up,downregulated=down,strong_upregulated=strong,
                    top100_upregulated=head(up,100),top100_strong_upregulated=head(strong,100))
    for (nm in names(outputs)) fwrite(outputs[[nm]],file.path(folder,paste0("edgeR_",nm,".csv")))
    save_diagnostic(function() edgeR::plotQLDisp(fit),file.path(folder,"QL_dispersion"))
    save_diagnostic(function() edgeR::plotMD(test),file.path(folder,"MD_plot"))
    base$genes_tested <- nrow(tab)
    base$significant_up <- nrow(up)
    base$significant_down <- nrow(down)
    base$status <- "completed"
    base
}

summaries <- list()
comparisons <- unique(original_assignments[,.(parent_cluster,subcluster)])
comparisons <- comparisons[order(as.integer(parent_cluster),as.integer(subcluster))]
for (i in seq_len(nrow(comparisons))) {
    parent_id <- comparisons$parent_cluster[i]
    sub_id <- comparisons$subcluster[i]
    folder <- file.path(out,paste0("parent_cluster_",parent_id),paste0("subcluster_",sub_id),"edgeR")
    dir.create(folder,recursive=TRUE,showWarnings=FALSE)
    cat("\nPARENT ",parent_id," / SUBCLUSTER ",sub_id,"\n",sep=""); flush.console()
    ans <- tryCatch(withCallingHandlers(fit_one(parent_id,sub_id,folder),warning=function(w) {
        cat(conditionMessage(w),"\n",file=file.path(folder,"warnings.log"),append=TRUE)
    }),error=function(e) data.frame(parent=parent_id,subcluster=sub_id,status="failed",error=conditionMessage(e)))
    fwrite(ans,file.path(folder,"cluster_summary.csv"))
    summaries[[i]] <- ans
    fwrite(rbindlist(summaries,fill=TRUE),file.path(out,"edgeR_all_subclusters_summary.csv"))
    print(ans)
    invisible(gc())
}
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
writeLines(c("edgeR QL F-test, legacy=FALSE, robust=TRUE; paired design ~ pair + group.",
             "Each target subcluster is compared only to other subclusters within the same parent.",
             "Minimum 20 cells in each group per donor; minimum 3 paired donors.",
             "TMM; filterByExpr; BH correction separately within each comparison.",
             "FDR<0.05 and absolute logFC>=0.25 are post-test filters; strong-up logFC>=1.",
             "CELLxGENE sources excluded following the previous raw-count workflow.",
             "117 audited unresolved-donor cells excluded: GSE221575=65; GSE231559=52.",
             "No additional tissue-site exclusions. Libraries pooled by study and biological donor.",
             "Marker discovery uses clusters learned from these data and is exploratory."),
           file.path(out,"METHOD.txt"))
cat("\nResults: ",out,"\n",sep="")
print(rbindlist(summaries,fill=TRUE))
