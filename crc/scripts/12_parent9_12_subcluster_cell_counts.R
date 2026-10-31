# Standalone count report: reads saved assignments and the exact exclusion audits.
library(data.table)
atlas <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
recluster_dir <- file.path(atlas,"CRC_individual_reclustering",
                          "run_20261030_155306_1a9ab35069091b")
blocked_run <- file.path(recluster_dir,"edgeR_1a9ab34daec974")
a <- rbindlist(lapply(c("9","12"),function(p) fread(file.path(
    recluster_dir,paste0("parent_cluster_",p),"cell_subcluster_assignments.csv"),
    colClasses="character")),fill=TRUE)
e <- rbindlist(lapply(c("GSE221575","GSE231559"),function(s) fread(file.path(
    blocked_run,paste0(s,"_unresolved_donor_cells.csv")),colClasses="character")),fill=TRUE)
stopifnot(!anyDuplicated(a$cell),!anyDuplicated(e$cell),nrow(e)==117L,
          all(e$cell %in% a$cell))
excluded_sources <- c("CELLxGENE_Chen_2021","CELLxGENE_Moorman2025")
a[,donor_excluded:=cell %in% e$cell]
a[,available_for_DE:=!donor_excluded & !study %in% excluded_sources]
counts <- a[,.(original_cells=.N,donor_excluded=sum(donor_excluded),
               after_donor_exclusion=sum(!donor_excluded),
               cells_available_for_DE=sum(available_for_DE),
               studies=uniqueN(study[available_for_DE])),
            by=.(parent_cluster,subcluster)]
counts <- counts[order(as.integer(parent_cluster),as.integer(subcluster))]
print(counts)
cat("\nTOTALS BY PARENT:\n")
print(counts[,.(original_cells=sum(original_cells),
               donor_excluded=sum(donor_excluded),
               cells_available_for_DE=sum(cells_available_for_DE)),by=parent_cluster])
report <- file.path(recluster_dir,"subcluster_cell_counts.csv")
fwrite(counts,report)
cat("\nSaved: ",report,"\n",sep="")
cat("Available cells are before per-comparison donor eligibility filtering.\n")
