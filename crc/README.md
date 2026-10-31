# CRC DE and Annotation Workflow

This folder contains the CRC-specific differential-expression, unresolved-cluster refinement, annotation, and final-output workflow for the colorectal cancer single-cell atlas.

It is intended to live as:

```text
crc-atlas-de-annotation/
├── adenoma/
└── crc/
```

## Structure

```text
crc/
├── README.md
├── config/
│   ├── CRC_unresolved_recluster_annotation_map.csv
│   └── CRC_parent9_12_subcluster_annotation_map.csv
└── scripts/
    ├── 01_CRC_original_CCA_pseudobulk_DESeq2.R
    ├── 02_CRC_unresolved_CCA_reclustering.R
    ├── 03_CRC_unresolved_recluster_DESeq2_meanfit.R
    ├── 04_CRC_recluster1_edgeR_glmLRT.R
    ├── 05_CRC_reclusters_2_5_6_8_10_limma_voom.R
    ├── 06_CRC_recluster16_Wilcoxon_annotation_only.R
    ├── 07_CRC_recluster18_Wilcoxon_annotation_only.R
    ├── 08_CRC_merge_recluster_annotations_into_full_object.R
    ├── 09_CRC_final_numbered_UMAP_and_export.R
    ├── 10_restore_GSE225857_donors_for_parent9_12.R
    ├── 11_edgeR_parent9_12_subclusters.R
    └── 12_parent9_12_subcluster_cell_counts.R
```

## Main workflow

### 1. Original CRC CCA cluster pseudobulk DE

`01_CRC_original_CCA_pseudobulk_DESeq2.R`

Uses raw RNA counts and original `cca_snn_res.0.5` assignments. Counts are aggregated by biological replicate and cluster. Each target CCA cluster is tested against the remaining CCA clusters using a paired design.

Processed CELLxGENE matrices are excluded from raw-count formal DE.

### 2. Reclustering the originally unresolved CRC CCA clusters

`02_CRC_unresolved_CCA_reclustering.R`

Original unresolved global CCA clusters:

```text
3, 5, 8, 9, 10, 15, 16, 18, 19, 20, 25, 26, 27, 28
```

The subset contains 270,276 cells. Existing RNA HVGs are reused; VST is not rerun on the split Seurat-v5 subset. PCA uses 20 dimensions and the unresolved cells are reintegrated with CCA. The resulting local unresolved IDs are 0-18 in `unresolved_cluster_res0.5`.

### 3. Formal DE for local unresolved reclusters

`03_CRC_unresolved_recluster_DESeq2_meanfit.R`

Primary paired pseudobulk DESeq2 analysis for local unresolved reclusters. Uses the validated raw pseudobulk cache.

Common thresholds:

```text
minimum target cells / paired unit = 20
minimum rest cells / paired unit   = 20
minimum paired units               = 3
adjusted P / FDR                   < 0.05
significant up                     log2FC >= 0.25
significant down                   log2FC <= -0.25
strong up                          log2FC >= 1
```

Corrected DESeq2 fitting uses `sfType="poscounts"` and `fitType="mean"`.

### 4. Numerical fallbacks

`04_CRC_recluster1_edgeR_glmLRT.R`

Recluster 1 fallback:

```text
edgeR DGEList
filterByExpr
TMM
estimateDisp(robust=FALSE)
glmFit
glmLRT
```

`05_CRC_reclusters_2_5_6_8_10_limma_voom.R`

Low-memory fallback for reclusters 2, 5, 6, 8, and 10:

```text
filtered raw counts
TMM
paired design
edgeR::voomLmFit
limma::eBayes
```

This fallback does not use `estimateDisp()` or `glmLRT()`.

### 5. Cell-level annotation support

`06_CRC_recluster16_Wilcoxon_annotation_only.R`

Used to refine local recluster 16 because formal pseudobulk DE provided too few informative markers.

`07_CRC_recluster18_Wilcoxon_annotation_only.R`

Recluster 18 had only one eligible paired biological replicate and therefore no formal replicate-aware DE. Wilcoxon is used only for descriptive annotation.

These cell-level P values must not be represented as biological-replicate DE evidence.

## Current local unresolved-recluster annotations

See:

```text
config/CRC_unresolved_recluster_annotation_map.csv
```

Parent-level labels 9 and 12 are currently:

```text
B-cell / plasma-cell population
```

Local recluster 16 was refined to:

```text
DKK1+/MUC5AC+ stress-associated epithelial cells
```

Local recluster 18 was annotated:

```text
LGR5+ absorptive progenitor-like epithelial cells
```

## 6. Merge refined annotations into the full CRC object

`08_CRC_merge_recluster_annotations_into_full_object.R`

Resolved original CCA cells retain their original annotation. Only cells belonging to the originally unresolved CCA clusters are replaced with their local recluster-derived annotation.

The main refined field is:

```text
CCA_annotation_refined
```

The script also stores:

```text
unresolved_recluster_id
unresolved_recluster_annotation
annotation_source
```

The two requested simplified original labels are:

```text
Dendritic cells
Epithelial cells
```

rather than donor/mixed-subset wording.

## 7. Final numbering, UMAP, and export

`09_CRC_final_numbered_UMAP_and_export.R`

Creates the numbered CRC UMAP with the cell-type legend and saves final outputs under:

```text
<atlas_root>/final_annotated_conditions/
```

including:

```text
annotated_CRC.rds
annotated_CRC_UMAP.png
annotated_CRC_UMAP.pdf
annotated_CRC_annotation_key.csv
```

Normal and TAN can later be added to the same final-condition directory.

## Parent 9/12 deeper resolution branch

The current parent-9/12 subcluster analysis is preserved separately because it is a later refinement branch.

### Donor recovery

`10_restore_GSE225857_donors_for_parent9_12.R`

Restores GSE225857 biological donors from the original GEO cell-level metadata. Donors come from the GEO `patients` field, not SampleTags or pooled immune/nonimmune GSM titles. The script preserves donor provenance and rejects missing or ambiguous mappings.

The original metadata location used on OSCER is:

```text
<atlas_root>/GSE225857_original_metadata/
```

### edgeR subcluster DE

`11_edgeR_parent9_12_subclusters.R`

Uses the saved parent-9 and parent-12 subcluster assignments. The final method is paired edgeR quasi-likelihood:

```text
target subcluster vs other subclusters within the same parent
design = ~ pair + group
TMM
filterByExpr
glmQLFit(legacy=FALSE, robust=TRUE)
glmQLFTest
```

The analysis excludes exactly the previously audited unresolved-donor cells:

```text
GSE221575 = 65 cells
GSE231559 = 52 cells
```

and retains the earlier CELLxGENE raw-count exclusions.

The uploaded final script records this explicitly. fileciteturn39file1L407-L418

### Subcluster cell counts

`12_parent9_12_subcluster_cell_counts.R`

Reports original cells, donor exclusions, cells available for raw-count DE, and study representation by parent/subcluster. fileciteturn39file2L14-L30

### Current subcluster annotation map

See:

```text
config/CRC_parent9_12_subcluster_annotation_map.csv
```

This map includes only subclusters already reviewed. It should be extended as additional subcluster marker tables are annotated.

## Reproducibility

Large Seurat objects, raw-count matrices, pseudobulk caches, fitted model objects, and generated result directories remain on OSCER and should not be committed to GitHub.

The repository should contain:

- analysis scripts;
- annotation maps/configuration;
- README/method documentation.

The GSE225857 donor restoration script explicitly uses the saved original GEO metadata files and maps the `patients` field to `verified_donor_id`, while rejecting ambiguous cell matches. fileciteturn39file0L55-L115

## Author

Olajumoke B. Oladapo  
Ph.D. Biomedical Engineering  
University of Oklahoma  
Moussa Lab
