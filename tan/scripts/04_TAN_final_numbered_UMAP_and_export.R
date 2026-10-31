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

condition <- "TAN"
title_text <- "TAN"
object_file <- "annotated_TAN.rds"

out <- file.path(atlas_root,"final_annotated_conditions")
obj <- readRDS(file.path(out,object_file))

reduction <- if ("umap.cca" %in% Reductions(obj)) {
    "umap.cca"
} else if ("umap" %in% Reductions(obj)) {
    "umap"
} else {
    stop("No UMAP reduction found.")
}

annotation <- as.character(obj$CCA_annotation_refined)
final_order <- sort(unique(annotation[!is.na(annotation) & nzchar(annotation)]))

number_map <- data.table(
    final_cluster_number=0:(length(final_order)-1),
    annotation=final_order
)
number_map[,legend_label:=paste0(final_cluster_number,": ",annotation)]

obj$final_cluster_number <- number_map$final_cluster_number[
    match(annotation,number_map$annotation)
]
obj$final_cluster_label <- number_map$legend_label[
    match(annotation,number_map$annotation)
]

um <- Embeddings(obj,reduction)
plot_dt <- data.table(
    UMAP_1=um[,1],
    UMAP_2=um[,2],
    final_cluster_number=obj$final_cluster_number,
    final_cluster_label=obj$final_cluster_label
)

plot_dt[,final_cluster_label:=factor(
    final_cluster_label,
    levels=number_map$legend_label
)]

label_dt <- plot_dt[
    ,
    .(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),
    by=.(final_cluster_number,final_cluster_label)
]

cols <- grDevices::hcl.colors(
    length(levels(plot_dt$final_cluster_label)),
    palette="Dark 3",
    alpha=1
)
names(cols) <- levels(plot_dt$final_cluster_label)

p <- ggplot(
    plot_dt,
    aes(UMAP_1,UMAP_2,color=final_cluster_label)
) +
    geom_point(size=0.20,alpha=1,stroke=0) +
    geom_text(
        data=label_dt,
        aes(UMAP_1,UMAP_2,label=final_cluster_number),
        inherit.aes=FALSE,
        color="black",
        fontface="bold",
        size=5
    ) +
    scale_color_manual(values=cols,drop=FALSE) +
    labs(
        title=title_text,
        x="UMAP 1",
        y="UMAP 2",
        color="Cell type"
    ) +
    theme_classic(base_size=14) +
    theme(
        plot.title=element_text(hjust=0.5,face="bold",size=20),
        axis.title=element_text(face="bold",size=15),
        axis.text=element_text(size=11,color="black"),
        legend.position="right",
        legend.title=element_text(face="bold",size=14),
        legend.text=element_text(size=9,color="black"),
        legend.key.height=unit(0.42,"cm"),
        legend.key.width=unit(0.42,"cm")
    ) +
    guides(
        color=guide_legend(
            override.aes=list(size=4,alpha=1),
            ncol=2
        )
    )

print(p)

fwrite(
    number_map,
    file.path(out,paste0("annotated_",condition,"_annotation_key.csv"))
)

saveRDS(
    obj,
    file.path(out,object_file),
    compress=FALSE
)

ggsave(
    file.path(out,paste0("annotated_",condition,"_UMAP.png")),
    p,width=19,height=11,dpi=400,bg="white"
)

ggsave(
    file.path(out,paste0("annotated_",condition,"_UMAP.pdf")),
    p,width=19,height=11,bg="white"
)
