# Adenoma/Polyp Differential Expression and Cell-Type Annotation

This repository contains the reproducible workflow used to perform differential-expression analysis and finalize annotations for the **already-reclustered unresolved Adenoma/Polyp cell populations** in the colorectal cancer single-cell atlas.

> **Important:** this workflow does **not** perform reclustering. It starts from the existing reclustered Adenoma/Polyp unresolved-cell Seurat object and performs DE, marker inspection, annotation, and final UMAP visualization.

## Workflow overview

The analysis separates **formal biological-replicate differential expression** from **annotation-oriented marker discovery**.

1. Load the existing reclustered Adenoma/Polyp unresolved-cell Seurat object.
2. Recover raw RNA counts from the original study-level Seurat objects.
3. Restrict the formal Adenoma analysis to appropriate biological samples.
4. Build sample-by-cluster raw-count pseudobulks.
5. Run paired pseudobulk DESeq2 likelihood-ratio tests for each cluster versus the remaining unresolved populations.
6. Use a sample-level paired Wilcoxon analysis for cluster 1 because its DESeq2/edgeR marker profile was weak/unstable.
7. Use cell-level Seurat Wilcoxon marker discovery for clusters 9 and 14 because these populations are overwhelmingly contributed by a single biological sample.
8. Apply the final annotation map.
9. Generate the final annotated UMAP and save the annotated Seurat object.

## Repository structure

```text
adenoma_de_annotation_workflow/
├── README.md
└── scripts/
    └── 01_adenoma_de_annotation_workflow.R
```

## Required R packages

```r
library(Seurat)
library(Matrix)
library(data.table)
library(DESeq2)
library(ggplot2)
```

The script saves `sessionInfo()` at the end of the run so package versions can be documented with the analysis.

## Input data

The default atlas root is:

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas
```

The existing Adenoma/Polyp reclustered object is expected under:

```text
Unresolved_Reclustering_By_Condition/Adenoma_polyp/
```

The script searches these filenames in order:

```text
Adenoma_polyp_UNRESOLVED_RECLUSTERED_FINAL_ANNOTATED.rds
Adenoma_polyp_UNRESOLVED_CCA_RECLUSTERED.rds
Adenoma_polyp_UNRESOLVED_RECLUSTERED.rds
```

The recluster identity column is:

```text
unresolved_cluster_res0.5
```

### Raw-count source objects

Formal pseudobulk DE does **not** use the integrated Seurat object's count layer. Raw counts are recovered from study-level objects under:

```text
full_seurat_QC_normalized_by_study_condition/
```

The formal pseudobulk analysis uses:

- `GSE201348`
- `GSE261388`

## GSE261388 sample filtering

GSE261388 contains both inflammatory and adenomatous polyps.

| Sample | Type | Formal Adenoma DE |
|---|---|---|
| P1 | Inflammatory polyp | Excluded |
| P12 | Inflammatory polyp | Excluded |
| P3 | Adenomatous polyp | Included |
| P4 | Adenomatous polyp | Included |

For GSE261388, sample identity is recovered from the cell barcode, for example:

```text
GSE261388_P3_AAAC...
GSE261388_P4_AAAC...
```

## Pseudobulk construction

Raw counts are summed for every:

```text
biological sample × new cluster
```

A replicate is defined as:

```r
replicate = paste0(study, "__", sample)
```

and the pseudobulk identifier is:

```r
group_id = paste0(replicate, "__CLUSTER__", cluster)
```

Generated core files:

```text
DESeq2_ADENOMA_ONLY_LRT/
├── ADENOMA_ONLY_raw_pseudobulk_counts.rds
└── ADENOMA_ONLY_pseudobulk_metadata.csv
```

## Primary formal DE: DESeq2 paired pseudobulk LRT

For each reclustered population, the target cluster is compared with all other unresolved populations from the **same biological samples**.

Samples are retained when they contain at least:

```text
10 target-cluster cells
10 rest cells
```

Genes are retained for testing when they have at least:

```text
10 raw counts in >= 3 target pseudobulks
and
10 raw counts in >= 3 rest pseudobulks
```

### Model

Full model:

```r
~ replicate + group
```

Reduced model:

```r
~ replicate
```

Likelihood-ratio test:

```r
DESeq(
    dds,
    test = "LRT",
    reduced = ~ replicate,
    fitType = "local"
)
```

The contrast is `target` versus `rest`.

The primary positive-marker threshold used for annotation is:

```text
adjusted P < 0.05
log2FC > 1
```

## Why pseudobulk?

Cells from the same specimen are not independent biological replicates. Pseudobulk aggregation preserves the biological sample as the replicate and avoids treating thousands of individual cells from one sample as thousands of independent observations.

## Cluster 1: additional paired Wilcoxon analysis

Cluster 1 was broadly represented across biological samples but did not yield a useful positive marker program with the initial DESeq2/edgeR rescue attempts.

Two additional analyses were used.

### Descriptive CPM enrichment

```text
log2((cluster1 CPM + 1) / (rest CPM + 1))
```

This analysis is descriptive and does not produce a p-value.

Representative enriched genes included:

```text
CA1
CA2
PCK1
GPT
SLC4A4
SLC26A2
SELENBP1
SCNN1B
CES3
NR5A2
MAOA
```

### Paired biological-sample Wilcoxon

A paired Wilcoxon signed-rank test is performed across eligible biological samples using:

```text
log2(CPM + 1)
```

for cluster 1 and the rest of the unresolved populations within each sample.

Final cluster-1 annotation:

```text
CA1/CA2+ differentiated absorptive colonocyte-like cells
```

## Clusters 9 and 14: annotation-only Wilcoxon

Clusters 9 and 14 are overwhelmingly contributed by `GSM6061707`.

Because they effectively lack independent biological replication, DESeq2 or edgeR pseudobulk p-values are **not appropriate for formal inferential DE** for these populations.

Instead, Seurat cell-level Wilcoxon tests are used strictly for marker ranking and cell-type annotation:

```r
FindMarkers(
    unresolved,
    ident.1 = cluster,
    ident.2 = NULL,
    assay = "RNA",
    slot = "data",
    test.use = "wilcox",
    logfc.threshold = 0,
    min.pct = 0
)
```

These cell-level p-values must **not** be interpreted as independent-sample DE evidence.

### Cluster 9

Representative markers:

```text
OLFM4
DACH1
FOXP2
FGF13
BMPR1B
SEMA3A
MAGI2
ENTPD1
```

Final annotation:

```text
OLFM4+ intestinal stem/progenitor epithelial cells
```

### Cluster 14

Representative markers:

```text
OLFM4
DIAPH3
ARHGAP11B
KIF23
UBE2C
NUSAP1
KNL1
PCLAF
CENPP
TPX2
CDCA2
PTTG1
TOP2A
CENPF
BUB1B
```

Final annotation:

```text
Cycling OLFM4+ stem/progenitor epithelial cells
```

## Final annotation map

| Cluster | Final annotation |
|---:|---|
| 0 | LGR5+ stem/progenitor-like epithelial cells |
| 1 | CA1/CA2+ differentiated absorptive colonocyte-like cells |
| 2 | Goblet/secretory epithelial cells |
| 3 | Cycling/proliferating cells |
| 4 | Inflammatory/regenerative crypt epithelial cells |
| 5 | Mature absorptive colonocytes |
| 6 | Activated/plasmablast-like B cells |
| 7 | Smooth muscle/myofibroblast-like stromal cells |
| 8 | POU2F3+ tuft-like/tuft-precursor cells |
| 9 | OLFM4+ intestinal stem/progenitor epithelial cells |
| 10 | WNT-high intestinal stem/progenitor epithelial cells |
| 11 | Tuft cells |
| 12 | BEST4+ OTOP2+ absorptive epithelial cells |
| 13 | G2/M-high proliferating cells |
| 14 | Cycling OLFM4+ stem/progenitor epithelial cells |

## Representative marker programs for the remaining clusters

| Cluster | Representative markers |
|---:|---|
| 0 | LGR5, EDAR, ALDH1B1, SAA1 |
| 2 | MUC2, CLCA1, ITLN1, REG4, FCGBP, RETNLB, B3GNT6 |
| 3 | MKI67, TOP2A, ASPM, BUB1, BUB1B, CENPE, KIF20A |
| 4 | OLFM4, REG3A, CLDN2, DEFA5, LCN2, S100P, CEACAM6 |
| 5 | SLC6A19, MS4A12, SLC15A1, SLC9A3, CA4, CEACAM7, ALDOB |
| 6 | FCRL5, MZB1, IRF4, CD38, DERL3, IGHA2, IGKC, IGHG1 |
| 7 | CARMN, LMOD1, CALD1, COL12A1, POSTN, VCAN, DCN |
| 8 | POU2F3, PSTPIP2, ST18, PLCG2 |
| 10 | NOTUM, SMOC2, AXIN2, EPHB3, PTK7, PROX1 |
| 11 | POU2F3, TRPM5, HPGDS, RGS13, PLCB2, KIT |
| 12 | BEST4, OTOP2, CA4, GUCA2B, LYPD8 |
| 13 | KIF18B, KIFC1, KIF2C, AURKB, BIRC5, CDC20, NDC80 |

## Major output files

```text
DESeq2_ADENOMA_ONLY_LRT/
    ADENOMA_ONLY_raw_pseudobulk_counts.rds
    ADENOMA_ONLY_pseudobulk_metadata.csv
    ADENOMA_ONLY_DESeq2_LRT_run_summary.csv
    Adenoma_cluster_<N>_VS_rest_DESeq2_LRT_ALL.csv
    Adenoma_cluster_<N>_VS_rest_DESeq2_LRT_UP.csv
    Adenoma_cluster_<N>_VS_rest_DESeq2_LRT_DOWN.csv
    Adenoma_cluster_<N>_VS_rest_DESeq2_LRT_TOP100_UP.csv

Cluster1_log2CPM_and_paired_Wilcoxon/
    Adenoma_cluster1_VS_rest_log2_CPM_ratio_ALL.csv
    Adenoma_cluster1_VS_rest_log2_CPM_ratio_TOP100.csv
    Adenoma_cluster1_VS_rest_paired_Wilcoxon_ALL.csv
    Adenoma_cluster1_VS_rest_paired_Wilcoxon_FDR_UP.csv

Wilcoxon_cluster9_cluster14/
    Adenoma_cluster9_VS_ALL_Wilcoxon_ALL.csv
    Adenoma_cluster9_VS_ALL_Wilcoxon_UP.csv
    Adenoma_cluster14_VS_ALL_Wilcoxon_ALL.csv
    Adenoma_cluster14_VS_ALL_Wilcoxon_UP.csv

Final_annotation_plots/
    Adenoma_FINAL_ANNOTATED_UMAP.pdf
    Adenoma_FINAL_ANNOTATED_UMAP.png
```

The final annotated Seurat object is saved as:

```text
Adenoma_polyp_UNRESOLVED_RECLUSTERED_FINAL_ANNOTATED.rds
```

## UMAP interpretation

The final plot uses the existing UMAP generated from the CCA-integrated representation. The axes are labeled:

```text
UMAP 1
UMAP 2
```

These are embedding coordinates and do not represent physical biological units.

## Statistical interpretation

The workflow intentionally distinguishes three levels of evidence:

1. **Formal biological-replicate DE** — paired pseudobulk DESeq2 LRT.
2. **Sample-level rescue analysis** — paired Wilcoxon across biological samples for cluster 1.
3. **Annotation-only cell-level marker discovery** — Seurat Wilcoxon for clusters 9 and 14.

## Running the workflow

From R:

```r
source("scripts/01_adenoma_de_annotation_workflow.R")
```

or from the shell:

```bash
Rscript scripts/01_adenoma_de_annotation_workflow.R
```

Before running on another system, edit:

```r
atlas_root <- "/ourdisk/hpc/moussalab/jumi/dont_archive/atlas"
```

to match the local project path.

## Reproducibility

The script writes:

```text
Adenoma_DE_annotation_sessionInfo.txt
```

at the end of the workflow. Commit that file with analysis results if exact package-version tracking is needed.

## Project note

This workflow was developed for the Adenoma/Polyp component of a colorectal cancer single-cell atlas analysis. When adapting the code to another dataset, preserve biological-replicate structure and reassess study-specific sample filtering rather than copying sample IDs or thresholds blindly.
