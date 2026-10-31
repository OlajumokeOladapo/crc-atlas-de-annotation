# Restore selected GSE225857 cells from the ORIGINAL GEO cell metadata.
# No DE, clustering, downloads, or source-object overwrites are performed.
atlas_dir <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
recluster_dir <- file.path(atlas_dir,"CRC_individual_reclustering",
                          "run_20261030_155306_1a9ab35069091b")

# Original GEO metadata already saved on OSCER.
metadata_restore_dir <- file.path(atlas_dir,"GSE225857_original_metadata")
metadata_files <- c(
    GSM7058754=file.path(metadata_restore_dir,"GSM7058754_immune_meta.txt.gz"),
    GSM7058755=file.path(metadata_restore_dir,"GSM7058755_non_immune_meta.txt.gz")
)
print(data.frame(library=names(metadata_files),file=unname(metadata_files),
                 exists=file.exists(metadata_files)))
if(!all(file.exists(metadata_files))) stop("Original metadata missing at the explicit paths above. No download attempted.")

clean_cell <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^GSM[0-9]+_", "", x)
    gsub("(SampleTag[0-9]+)\\.", "\\1-", x, perl=TRUE)
}
cell_variants <- function(x) {
    x <- clean_cell(x)
    variants <- x
    # Seurat merge suffixes are considered only as explicit alternative keys.
    # Multiple matching metadata rows are an error, never resolved arbitrarily.
    while (grepl("_[0-9]+$", x)) {
        x <- sub("_[0-9]+$", "", x)
        variants <- c(variants,x)
    }
    unique(variants)
}
match_cells <- function(cells, metadata_ids, library) {
    keys <- paste(library,clean_cell(metadata_ids),sep="__")
    vapply(cells,function(cell) {
        query <- paste(library,cell_variants(cell),sep="__")
        hits <- unique(na.omit(match(query,keys)))
        if(length(hits)>1L) return(-1L)
        if(!length(hits)) return(NA_integer_)
        as.integer(hits)
    },integer(1),USE.NAMES=FALSE)
}


selected <- do.call(rbind,lapply(c("9","12"),function(parent) {
    read.csv(file.path(recluster_dir,paste0("parent_cluster_",parent),
                       "cell_subcluster_assignments.csv"),colClasses="character",check.names=FALSE)
}))
stopifnot(all(c("cell","study") %in% names(selected)))
selected <- selected[selected$study=="GSE225857",,drop=FALSE]
if(!nrow(selected) || anyNA(selected$cell) || anyDuplicated(selected$cell)) stop("Invalid selected cell IDs.")
selected$source_library <- sub("_.*$","",selected$cell)
if(any(!selected$source_library %in% names(metadata_files))) stop("Unexpected source-library prefix.")

recovery_dir <- tempfile(pattern="GSE225857_donor_recovery_",tmpdir=recluster_dir)
dir.create(recovery_dir)
recovered <- list()
for(lib in names(metadata_files)) {
    cells <- selected$cell[selected$source_library==lib]
    if(!length(cells)) next
    con <- gzfile(metadata_files[[lib]],open="rt")
    md <- tryCatch(read.delim(con,header=TRUE,row.names=NULL,check.names=FALSE,
                             colClasses="character",quote='"',comment.char=""),finally=close(con))
    patient_col <- which(tolower(trimws(names(md)))=="patients")
    if(length(patient_col)!=1L) stop(lib,": expected exactly one patients column. Columns: ",paste(names(md),collapse=", "))

    # Identify the barcode column by actual matches, not its spelling.
    matches <- lapply(md,function(values) match_cells(cells,values,lib))
    coverage <- vapply(matches,function(ix) sum(!is.na(ix) & ix>0L),integer(1))
    best <- which(coverage==max(coverage))
    if(!length(best) || max(coverage)==0L) stop(lib,": no metadata column matches selected cell IDs.")
    first <- best[1]
    if(length(best)>1L && !all(vapply(best,function(j)
        identical(clean_cell(md[[first]]),clean_cell(md[[j]])),logical(1)))) {
        stop(lib,": ambiguous barcode-column identification.")
    }
    ids <- clean_cell(md[[first]])
    if(anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(paste(lib,ids,sep="__")))
        stop(lib,": missing or duplicate normalized metadata match keys.")
    ix <- matches[[first]]
    if(anyNA(ix) || any(ix<1L)) {
        write.csv(data.frame(cell=cells,metadata_row=ix),file.path(recovery_dir,paste0(lib,"_match_audit.csv")),row.names=FALSE)
        stop(lib,": unmatched or ambiguous selected cells; see match audit.")
    }
    if(anyDuplicated(ix)) stop(lib,": multiple selected cells map to the same source metadata row.")
    patients <- trimws(md[[patient_col]][ix])
    if(anyNA(patients) || any(tolower(patients) %in% c("","na","n/a","unknown","null")))
        stop(lib,": missing patients values in matched rows.")

    result <- data.frame(study="GSE225857",cell=cells,verified_donor_id=patients,
                         donor_source_field="patients",donor_source_value=patients,
                         source_library=lib,source_metadata_file=metadata_files[[lib]],
                         source_cell_id=md[[first]][ix],
                         match_key=paste(lib,ids[ix],sep="__"),
                         matching_rule=ifelse(clean_cell(cells)==ids[ix],
                                              "GSM_prefix_and_SampleTag_punctuation",
                                              "also_removed_Seurat_merge_suffix"),
                         stringsAsFactors=FALSE)
    for(field in c("sampletag","organs","patients_organ")) {
        j <- which(tolower(trimws(names(md)))==field)
        result[[field]] <- if(length(j)==1L) md[[j]][ix] else NA_character_
    }
    recovered[[lib]] <- result
    cat(lib,": matched ",nrow(result)," cells; ",length(unique(patients))," patients\n",sep="")
}
result <- do.call(rbind,recovered)
stopifnot(nrow(result)==nrow(selected),!anyDuplicated(result$cell),setequal(result$cell,selected$cell))
map_file <- file.path(recovery_dir,"GSE225857_verified_donor_map.csv")
write.csv(result,map_file,row.names=FALSE)
write.csv(data.frame(source_library=names(metadata_files),file=unname(metadata_files),
                     md5=unname(tools::md5sum(metadata_files))),
          file.path(recovery_dir,"source_files.csv"),row.names=FALSE)
cat("\nVerified donor map:\n",map_file,"\n",sep="")
print(with(result,table(source_library,verified_donor_id)))
cat("\nUse this CSV as verified_map_file in the edgeR script.\n")
cat("GSE221575 and GSE231559 are not changed by this recovery script.\n")
