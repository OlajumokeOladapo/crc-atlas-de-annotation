############################################################
# CRC FINAL NUMBERING, UMAP, AND FINAL CONDITION EXPORT
############################################################

suppressPackageStartupMessages({
    library(Seurat)
    library(data.table)
    library(ggplot2)
    library(grid)
})

atlas_root <- Sys.getenv(
    "CRC_ATLAS_ROOT",
    unset="/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
)

input_file <- file.path(
    atlas_root,
    "Final_Annotated_Seurat_Objects",
    "CRC_CCA_FINAL_REFINED_RENUMBERED.rds"
)

final_dir <- file.path(
    atlas_root,
    "final_annotated_conditions"
)
dir.create(final_dir,recursive=TRUE,showWarnings=FALSE)

crc <- readRDS(input_file)

if (!"CCA_annotation_refined" %in% colnames(crc[[]])) stop("CCA_annotation_refined missing.")
if (!"umap.cca" %in% Reductions(crc)) stop("umap.cca missing.")

# Biology-first order for the recluster-derived populations.
preferred <- c(
    "LGR5/OLFM4+ epithelial stem/progenitor-like cells",
    "LGR5+ absorptive progenitor-like epithelial cells",
    "Cycling epithelial cells",
    "DKK1+/MUC5AC+ stress-associated epithelial cells",
    "Activated antigen-presenting B cells",
    "MHC-II-high B cells",
    "B-cell / plasma-cell population",
    "Plasma/plasmablast-like secretory B cells",
    "Plasma cells / plasmablast-like cells",
    "Cycling B-lineage / plasmablast-like cells",
    "GZMM+ regulatory T cells",
    "Activated regulatory T cells",
    "Gamma-delta T / NK-like cytotoxic lymphocytes",
    "Cycling gamma-delta T/NK-like cytotoxic lymphocytes",
    "ACKR1+/MADCAM1+ venous endothelial cells",
    "RGS5+/PDGFRB+ mural/pericyte cells",
    "COL2A1/COL11A2+ mesenchymal stromal cells",
    "INHBA+/BMP2+ vascular/perivascular stromal-like cells"
)

all_annotations <- unique(as.character(crc$CCA_annotation_refined))
all_annotations <- all_annotations[!is.na(all_annotations) & nzchar(all_annotations)]
order_final <- c(preferred[preferred %in% all_annotations],
                 sort(setdiff(all_annotations,preferred)))

number_map <- data.table(
    final_cluster_number=0:(length(order_final)-1),
    annotation=order_final
)
number_map[,legend_label:=paste0(final_cluster_number,": ",annotation)]

crc$final_cluster_number <- number_map$final_cluster_number[
    match(crc$CCA_annotation_refined,number_map$annotation)
]
if (anyNA(crc$final_cluster_number)) stop("Some cells lack final numbers.")

crc$final_cluster_label <- number_map$legend_label[
    match(crc$CCA_annotation_refined,number_map$annotation)
]

um <- Embeddings(crc,"umap.cca")
plot_dt <- data.table(
    cell=rownames(um),
    UMAP_1=um[,1],
    UMAP_2=um[,2],
    final_cluster_number=as.integer(crc$final_cluster_number),
    final_cluster_label=as.character(crc$final_cluster_label)
)
plot_dt[,final_cluster_label:=factor(
    final_cluster_label,
    levels=number_map$legend_label
)]

labels <- plot_dt[,.(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),
                  by=.(final_cluster_number,final_cluster_label)]

p_crc <- ggplot(plot_dt,aes(UMAP_1,UMAP_2,color=final_cluster_label)) +
    geom_point(size=0.15,alpha=0.9,stroke=0) +
    geom_text(data=labels,aes(UMAP_1,UMAP_2,label=final_cluster_number),
              inherit.aes=FALSE,color="black",fontface="bold",size=5) +
    labs(title="CRC",x="UMAP 1",y="UMAP 2",color="Cell type") +
    theme_classic(base_size=14) +
    theme(
        plot.title=element_text(hjust=0.5,face="bold",size=20),
        axis.title=element_text(face="bold",size=15),
        axis.text=element_text(size=11,color="black"),
        legend.position="right",
        legend.title=element_text(face="bold",size=14),
        legend.text=element_text(size=9),
        legend.key.height=unit(0.42,"cm"),
        legend.key.width=unit(0.42,"cm")
    ) +
    guides(color=guide_legend(override.aes=list(size=4,alpha=1),ncol=2))

fwrite(number_map,file.path(final_dir,"annotated_CRC_annotation_key.csv"))
saveRDS(crc,file.path(final_dir,"annotated_CRC.rds"),compress=FALSE)

ggsave(file.path(final_dir,"annotated_CRC_UMAP.png"),
       p_crc,width=19,height=11,dpi=400,bg="white")
ggsave(file.path(final_dir,"annotated_CRC_UMAP.pdf"),
       p_crc,width=19,height=11,bg="white")

capture.output(sessionInfo(),file=file.path(final_dir,"annotated_CRC_sessionInfo.txt"))

cat("\nFinal CRC outputs saved to:\n",final_dir,"\n",sep="")
