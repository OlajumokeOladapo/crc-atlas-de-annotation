# CRC Atlas Differential-Expression and Unresolved-Reclustering Workflow

This directory contains the **full CRC cluster differential-expression workflow through the current unresolved-recluster analysis stage**.

It is intentionally organized as a sibling of `adenoma/` in the existing `crc-atlas-de-annotation` repository.

## Folder structure

```text
crc/
├── README.md
└── scripts/
    ├── 01_CRC_original_CCA_pseudobulk_DESeq2.R
    ├── 02_CRC_unresolved_CCA_reclustering.R
    ├── 03_CRC_unresolved_recluster_DESeq2_meanfit.R
    ├── 04_CRC_recluster1_edgeR_glmLRT.R
    └── 05_CRC_reclusters_2_5_6_8_10_limma_voom.R
```

## What the workflow contains

The CRC analysis has **two distinct clustering levels**.

### Stage 1: Original global CRC CCA clusters

`01_CRC_original_CCA_pseudobulk_DESeq2.R`

The starting CRC CCA object is:

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
Integration_by_condition/CRC/CRC_CCA_integrated.rds
```

The original cluster column is:

```text
cca_snn_res.0.5
```

The CRC CCA object contains approximately:

```text
702,029 cells
```

Raw RNA counts from eligible source datasets are mapped back to these CCA cluster assignments and aggregated by biological replicate and CCA cluster.

For each original CCA cluster:

```text
target = selected CCA cluster
rest   = all other CCA clusters
```

The paired DESeq2 design is:

```r
~ replicate + group
```

The original validated workflow excludes:

```text
CELLxGENE_Chen_2021
CELLxGENE_Moorman2025
```

from formal raw-count DE because the atlas matrices available for those cohorts are processed/normalized rather than valid raw integer-count input.

The remaining raw-count studies are used for formal pseudobulk testing.

## Stage 2: Identify original unresolved CRC CCA clusters

After first-pass marker review and annotation of the original CCA clusters, the unresolved CRC clusters retained for a second round of clustering were:

```text
3
5
8
9
10
15
16
18
19
20
25
26
27
28
```

These numbers are **original global CCA cluster IDs**.

Together they contain:

```text
270,276 cells
```

## Stage 3: CRC-only unresolved-cell CCA reclustering

`02_CRC_unresolved_CCA_reclustering.R`

The reclustering workflow:

1. Loads the annotated full CRC CCA object.
2. Selects only cells belonging to the original unresolved CCA clusters.
3. Reuses the RNA variable features already defined on the full object.
4. Does **not** rerun `FindVariableFeatures()`/VST on the subset.
5. Scales the existing HVGs.
6. Runs a 20-PC PCA named:

```text
pca.unresolved
```

7. Reintegrates the unresolved subset with Seurat `CCAIntegration` using:

```text
integrated.cca.unresolved
```

8. Runs UMAP:

```text
umap.cca.unresolved
```

9. Builds graphs:

```text
CRC_unresolved_nn
CRC_unresolved_snn
```

10. Clusters at resolution:

```text
0.5
```

with the local cluster column:

```text
unresolved_cluster_res0.5
```

The saved object is:

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
Unresolved_Reclustering_By_Condition/CRC/
CRC_UNRESOLVED_CCA_RECLUSTERED_res0.5.rds
```

The resulting **local** unresolved recluster IDs are:

```text
0-18
```

These local IDs must not be confused with the original global CCA IDs.

The recorded recluster sizes were:

| Local recluster | Cells |
|---:|---:|
| 0 | 54,966 |
| 1 | 33,077 |
| 2 | 32,982 |
| 3 | 26,742 |
| 4 | 20,194 |
| 5 | 19,218 |
| 6 | 15,951 |
| 7 | 14,404 |
| 8 | 12,935 |
| 9 | 11,984 |
| 10 | 9,237 |
| 11 | 6,249 |
| 12 | 5,278 |
| 13 | 1,953 |
| 14 | 1,928 |
| 15 | 1,502 |
| 16 | 815 |
| 17 | 561 |
| 18 | 300 |

Total:

```text
270,276 cells
```

## Stage 4: Pseudobulk DE for the new unresolved reclusters

`03_CRC_unresolved_recluster_DESeq2_meanfit.R`

The validated raw-count pseudobulk cache is:

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5/
validated_pseudobulk_cache/
CRC_unresolved_reclusters_res0.5_raw_pseudobulk.rds
```

For each local recluster:

```text
target = selected local unresolved recluster
rest   = all other local unresolved reclusters
```

within each eligible biological replicate.

The paired design remains:

```r
~ replicate + group
```

The corrected DESeq2 rerun uses:

```text
sfType  = "poscounts"
fitType = "mean"
```

### Eligibility thresholds

Each replicate must contain at least:

```text
20 target cells
20 rest cells
```

and each formal comparison requires at least:

```text
3 paired replicate units
```

### Gene filtering

Genes must have:

```text
total count >= 10
```

and:

```text
>= 10 counts in at least 3 pseudobulk samples
```

### Differential-expression thresholds

```text
adjusted P / FDR < 0.05
|log2FC| >= 0.25
```

Strong upregulation:

```text
log2FC >= 1
```

Positive fold change always means higher expression in the target cluster.

## Stage 5: Fallback for local recluster 1

`04_CRC_recluster1_edgeR_glmLRT.R`

The corrected DESeq2 model was numerically unstable for recluster 1.

The successful fallback preserved the same paired target/rest construction and used:

```text
edgeR::DGEList
filterByExpr
TMM normalization
estimateDisp(robust = FALSE)
glmFit
glmLRT
```

The coefficient tested is:

```text
grouptarget
```

Positive `logFC` therefore means higher expression in recluster 1.

## Stage 6: Low-memory limma-voom fallback

`05_CRC_reclusters_2_5_6_8_10_limma_voom.R`

The remaining unstable comparisons were:

```text
2
5
6
8
10
```

For these clusters, the workflow retains:

- raw-count pseudobulk construction;
- biological-replicate pairing;
- expression filtering;
- edgeR `filterByExpr`;
- TMM normalization;
- the same paired `~ pair_id + group` design.

Differential expression is then performed with:

```text
edgeR::voomLmFit
+
limma::eBayes
```

This fallback intentionally does **not** run:

```text
estimateDisp()
glmLRT()
```

The script was written to reduce memory pressure on OSCER by:

- processing clusters sequentially;
- storing normalized checkpoints;
- releasing the large raw cache before fitting;
- writing model checkpoints;
- recording progress logs;
- and combining only successful results from the current batch.

## Local recluster 18

Recluster 18 had only:

```text
1 eligible paired replicate
```

and was therefore not used for formal replicate-aware differential expression.

Its eventual biological annotation must be treated as exploratory rather than as a formal pseudobulk DE result.

## Output locations

### Original global CCA DE

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
DESeq2_pseudobulk_CCA_res0.5_revised/
```

### Unresolved reclustering

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
Unresolved_Reclustering_By_Condition/CRC/
```

### Corrected unresolved DESeq2

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5_meanfit/
```

### Recluser-1 edgeR fallback

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
edgeR_pseudobulk_CRC_unresolved_reclusters_res0.5_fallback/
```

### Low-memory limma-voom fallback

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas/
edgeR_pseudobulk_CRC_unresolved_reclusters_res0.5_low_memory/
```

## Recommended execution order

The scripts represent the analysis history and should be read/run in this order:

```bash
Rscript scripts/01_CRC_original_CCA_pseudobulk_DESeq2.R
Rscript scripts/02_CRC_unresolved_CCA_reclustering.R
Rscript scripts/03_CRC_unresolved_recluster_DESeq2_meanfit.R
Rscript scripts/04_CRC_recluster1_edgeR_glmLRT.R
Rscript scripts/05_CRC_reclusters_2_5_6_8_10_limma_voom.R
```

The fallback scripts are not meant to replace successful primary results. They are used only for the reclusters for which the earlier model-fitting route was numerically unstable.

## Reproducibility notes

- Individual cells are never treated as independent biological replicates for formal DE.
- All formal comparisons are target-vs-rest within eligible replicate units.
- Raw-count pseudobulk is used for count-based formal testing.
- Processed CELLxGENE matrices are excluded from formal raw-count DE.
- Original global CCA IDs and local unresolved-recluster IDs are explicitly kept separate.
- Large Seurat objects, pseudobulk caches, DESeq2 model objects, checkpoints, and result tables should remain on OSCER `/ourdisk`; they should not be committed to GitHub.
- `sessionInfo()` is saved by the analysis scripts to document the software environment.

## Repository placement

This directory belongs in the existing repository as:

```text
crc-atlas-de-annotation/
├── adenoma/
└── crc/
```

---

**Author:** Olajumoke B. Oladapo  
Ph.D. Biomedical Engineering  
University of Oklahoma  
Moussa Lab
