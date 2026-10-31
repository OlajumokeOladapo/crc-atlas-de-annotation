# CRC: Full DE and Annotation Workflow

This folder documents the complete CRC analysis used in the colorectal cancer single-cell atlas, from the original CCA clusters through refinement of previously unresolved cells, additional resolution of heterogeneous B-cell/plasma-cell populations, and final integration into the annotated CRC atlas.

## Analysis flow

```text
Full CRC CCA object
~702,029 cells
cca_snn_res.0.5
        |
        v
Raw-count donor-level pseudobulk DESeq2
cluster vs all other CCA clusters
        |
        v
Original CCA annotation
resolved + unresolved populations
        |
        v
Extract only unresolved CCA cells
270,276 cells
        |
        v
CCA reclustering
new local reclusters 0-18
        |
        v
Replicate-aware pseudobulk DE
DESeq2 primary analysis
        |
        +-------------------------------+
        |                               |
        v                               v
edgeR fallback                  limma-voom fallback
recluster 1                     reclusters 2,5,6,8,10
        |
        v
Wilcoxon annotation support
reclusters 16 and 18
        |
        v
Final annotations for local reclusters 0-18
        |
        v
Additional resolution of heterogeneous
B-cell / plasma-cell populations
parents 9 and 12
        |
        v
Map refined cells back to full CRC object
by exact cell barcode
        |
        v
Final full CRC atlas
numbered UMAP + final annotated Seurat object
```

## Scripts

Run the scripts in order unless a script is explicitly described as an optional refinement branch.

### 01 — Original CRC CCA pseudobulk DESeq2

```text
scripts/01_CRC_original_CCA_pseudobulk_DESeq2.R
```

Starts from the condition-specific CRC CCA object:

```text
Integration_by_condition/CRC/CRC_CCA_integrated.rds
```

and tests each original `cca_snn_res.0.5` cluster against all remaining CRC CCA clusters using donor-level raw-count pseudobulks.

The paired model is:

```text
design = ~ replicate + group
```

with:

```text
target cells per donor >= 20
rest cells per donor   >= 20
eligible paired donors >= 3
gene total count       >= 10
padj cutoff            0.05
|log2FC| reporting     >= 0.25
```

Raw-count pseudobulk DE excludes the CELLxGENE sources whose available atlas matrices are processed/normalized rather than raw integer counts:

```text
CELLxGENE_Chen_2021
CELLxGENE_Moorman2025
```

### 02 — Reclustering of the original unresolved CRC cells

```text
scripts/02_CRC_unresolved_CCA_reclustering.R
```

The original CCA clusters retained as unresolved were:

```text
3, 5, 8, 9, 10, 15, 16, 18, 19, 20, 25, 26, 27, 28
```

Together these contain:

```text
270,276 cells
```

Important implementation choices:

- reuses the existing RNA variable features from the full CRC object,
- does **not** rerun VST/`FindVariableFeatures()` on the unresolved subset,
- scales the retained features,
- runs PCA with 20 PCs,
- uses `CCAIntegration`,
- stores the integrated reduction as `integrated.cca.unresolved`,
- stores the UMAP as `umap.cca.unresolved`,
- clusters at resolution 0.5,
- stores the local recluster identity in `unresolved_cluster_res0.5`.

The resulting local unresolved reclusters are:

```text
0-18
```

These are **local unresolved-reclustering IDs** and must not be confused with the original global CCA cluster IDs.

### 03 — Primary DE for the local unresolved reclusters

```text
scripts/03_CRC_unresolved_recluster_DESeq2_meanfit.R
```

Uses the validated raw-count pseudobulk cache:

```text
DESeq2_pseudobulk_CRC_unresolved_reclusters_res0.5/
validated_pseudobulk_cache/
CRC_unresolved_reclusters_res0.5_raw_pseudobulk.rds
```

For each local recluster:

```text
target = selected local recluster
rest   = all other local unresolved reclusters
```

within each eligible biological replicate.

The paired design is:

```text
~ replicate + group
```

The corrected DESeq2 fitting uses:

```text
sfType  = poscounts
fitType = mean
```

Eligibility and reporting thresholds:

```text
target cells per paired unit >= 20
rest cells per paired unit   >= 20
minimum paired units         >= 3
total gene count             >= 10
>=10 counts in >=3 pseudobulk libraries
adjusted P / FDR             < 0.05
significant up               log2FC >= 0.25
significant down             log2FC <= -0.25
strong up                    log2FC >= 1
```

### 04 — edgeR fallback for recluster 1

```text
scripts/04_CRC_recluster1_edgeR_glmLRT.R
```

DESeq2 was numerically unstable for local recluster 1. The successful fallback preserved the same paired target-versus-rest construction and used:

```text
edgeR DGEList
filterByExpr
TMM normalization
estimateDisp(robust = FALSE)
glmFit
glmLRT
```

Positive `logFC` indicates higher expression in the target recluster.

### 05 — limma-voom fallback for reclusters 2, 5, 6, 8, and 10

```text
scripts/05_CRC_reclusters_2_5_6_8_10_limma_voom.R
```

For these numerically unstable comparisons, the workflow retains:

- raw-count pseudobulk construction,
- biological-replicate pairing,
- expression filtering,
- `filterByExpr`,
- TMM normalization,
- the same paired target-versus-rest design.

Differential expression is then performed using:

```text
edgeR::voomLmFit
+
limma::eBayes
```

This fallback intentionally does **not** run `estimateDisp()` or `glmLRT()`.

### 06 — Wilcoxon annotation support for recluster 16

```text
scripts/06_CRC_recluster16_Wilcoxon_annotation_only.R
```

The formal pseudobulk result for recluster 16 yielded too few informative positive markers for confident biological annotation.

A cell-level Seurat Wilcoxon test was therefore used only for descriptive marker discovery.

The richer marker profile supported:

```text
DKK1+/MUC5AC+ stress-associated epithelial cells
```

This Wilcoxon analysis is annotation-supporting and is not treated as replicate-aware differential-expression inference.

### 07 — Wilcoxon annotation-only analysis for recluster 18

```text
scripts/07_CRC_recluster18_Wilcoxon_annotation_only.R
```

Recluster 18 had only one eligible paired biological replicate and therefore did not meet the minimum replication criterion for formal pseudobulk DE.

Cell-level Wilcoxon was used only for descriptive marker discovery and supported:

```text
LGR5+ absorptive progenitor-like epithelial cells
```

The Wilcoxon P values are not interpreted as biological-replicate DE evidence.

## Final annotations of CRC local unresolved reclusters

The final working map is stored in:

```text
config/CRC_unresolved_recluster_annotation_map.csv
```

| Local recluster | Annotation |
|---:|---|
| 0 | LGR5/OLFM4+ epithelial stem/progenitor-like cells |
| 1 | Plasma/plasmablast-like secretory B cells |
| 2 | Cycling epithelial cells |
| 3 | Activated antigen-presenting B cells |
| 4 | GZMM+ regulatory T cells |
| 5 | Plasma cells / plasmablast-like cells |
| 6 | Activated regulatory T cells |
| 7 | ACKR1+/MADCAM1+ venous endothelial cells |
| 8 | RGS5+/PDGFRB+ mural/pericyte cells |
| 9 | B-cell / plasma-cell population |
| 10 | Cycling gamma-delta T/NK-like cytotoxic lymphocytes |
| 11 | MHC-II-high B cells |
| 12 | B-cell / plasma-cell population |
| 13 | COL2A1/COL11A2+ mesenchymal stromal cells |
| 14 | Cycling B-lineage / plasmablast-like cells |
| 15 | Gamma-delta T / NK-like cytotoxic lymphocytes |
| 16 | DKK1+/MUC5AC+ stress-associated epithelial cells |
| 17 | INHBA+/BMP2+ vascular/perivascular stromal-like cells |
| 18 | LGR5+ absorptive progenitor-like epithelial cells |

## Additional refinement of parent reclusters 9 and 12

Local reclusters 9 and 12 were retained at the parent level as:

```text
B-cell / plasma-cell population
```

because their deeper reclustering revealed multiple B-cell/plasma-cell states together with smaller lineage-divergent components.

The deeper subcluster analysis is preserved as a separate refinement branch rather than replacing the parent labels in the main CRC atlas.

### 10 — Restore GSE225857 biological donors

```text
scripts/10_restore_GSE225857_donors_for_parent9_12.R
```

The GSE225857 immune/nonimmune GSM titles do not identify biological donors.

The donor-recovery workflow therefore uses the original GEO cell-level metadata files:

```text
GSM7058754_immune_meta.txt.gz
GSM7058755_non_immune_meta.txt.gz
```

stored under:

```text
<atlas_root>/GSE225857_original_metadata/
```

Cell identifiers are normalized and matched back to the original GEO metadata. Biological donor identity is taken from:

```text
patients
```

and written as:

```text
verified_donor_id
```

SampleTags and pooled immune/nonimmune GSM IDs are **not** treated as donor IDs.

The recovery script rejects ambiguous or unmatched cell mappings and preserves mapping provenance.

### 11 — edgeR DE for parent 9 and 12 subclusters

```text
scripts/11_edgeR_parent9_12_subclusters.R
```

The final subcluster workflow:

- uses the saved subcluster assignments,
- pools libraries at study + biological-donor level,
- compares each target subcluster only with other subclusters from the same parent,
- uses paired design `~ pair + group`,
- applies TMM normalization,
- applies `filterByExpr`,
- fits edgeR quasi-likelihood models with:

```text
glmQLFit(legacy = FALSE, robust = TRUE)
glmQLFTest
```

The thresholds are:

```text
target cells per donor >= 20
rest cells per donor   >= 20
paired donors          >= 3
FDR                    < 0.05
|logFC|                >= 0.25
strong up              logFC >= 1
```

The previous processed CELLxGENE exclusions are retained.

The final audited unresolved-donor exclusions are:

```text
GSE221575 = 65 cells
GSE231559 = 52 cells
total     = 117 cells
```

No additional tissue-site exclusions are applied in this subcluster DE step.

### 12 — Parent 9/12 subcluster cell-count report

```text
scripts/12_parent9_12_subcluster_cell_counts.R
```

Reports:

- original cell count,
- donor-excluded cells,
- cells remaining after donor exclusion,
- cells available for raw-count DE,
- number of contributing studies.

The output is:

```text
subcluster_cell_counts.csv
```

The counts precede the per-comparison 20-cells-per-group and minimum-three-paired-donors eligibility filters.

## Parent 9 subcluster annotations

The current reviewed subclusters are:

| Parent | Subcluster | Annotation |
|---:|---:|---|
| 9 | 0 | IgA+ plasma cells / plasmablasts |
| 9 | 1 | RGS5+/PDGFRB+ activated mural stromal cells |
| 9 | 2 | IGL-high B-lineage cells |
| 9 | 3 | Unresolved mixed T-cell/epithelial population |
| 9 | 4 | Vascular endothelial cells |
| 9 | 5 | Activated regulatory T cells |
| 9 | 6 | IgG+ plasma/plasmablast-like cells |
| 9 | 7 | Activated antigen-presenting B cells |
| 9 | 8 | C1Q+ inflammatory macrophages / TAM-like cells |

The remaining very small subclusters are not used to redefine the parent population.

The parent-level atlas label remains:

```text
B-cell / plasma-cell population
```

## Parent 12 subcluster annotations

The reviewed parent-12 subclusters currently include:

| Parent | Subcluster | Annotation |
|---:|---:|---|
| 12 | 0 | IgA+ plasma cells / plasmablasts |
| 12 | 1 | IgG+ plasma/plasmablast-like cells |
| 12 | 2 | IGL-high B-lineage cells |
| 12 | 3 | Activated antigen-presenting B cells |
| 12 | 4 | CCL4+ activated lymphocyte-like cells |

The parent-level atlas label remains:

```text
B-cell / plasma-cell population
```

The subcluster annotation map is stored in:

```text
config/CRC_parent9_12_subcluster_annotation_map.csv
```

and can be extended as additional reviewed subclusters are finalized.

### 08 — Merge refined unresolved-cell annotations into the full CRC atlas

```text
scripts/08_CRC_merge_recluster_annotations_into_full_object.R
```

This step does **not** rename an entire old unresolved CCA cluster.

Instead, it:

1. loads the full first-pass annotated CRC object,
2. loads the unresolved-cell reclustered object,
3. matches cells by exact barcode,
4. preserves original annotations for cells that were already resolved,
5. assigns each formerly unresolved cell its local recluster annotation,
6. preserves the original `cca_snn_res.0.5` cluster,
7. records the local recluster identity,
8. records annotation provenance,
9. saves the refined full CRC object.

Two original first-pass labels were simplified for the final atlas:

```text
Dendritic cells — mixed activated subsets
        ->
Dendritic cells

Donor-limited lineage-divergent epithelial
        ->
Epithelial cells
```

The main final biological annotation field is:

```text
CCA_annotation_refined
```

### 09 — Final numbered UMAP and export

```text
scripts/09_CRC_final_numbered_UMAP_and_export.R
```

This script:

- orders the refined annotations for display,
- assigns final numeric labels,
- creates the numbered CRC UMAP,
- displays the full cell-type legend,
- saves the final annotated object and plot.

Final outputs are written under:

```text
<atlas_root>/final_annotated_conditions/
```

with simple names:

```text
annotated_CRC.rds
annotated_CRC_UMAP.png
annotated_CRC_UMAP.pdf
annotated_CRC_annotation_key.csv
```

The same directory is intended to eventually contain:

```text
annotated_adenoma.rds
annotated_normal.rds
annotated_TAN.rds
```

and their corresponding final UMAPs.


## Final CRC cluster numbering

The final CRC atlas contains **33 numbered cell populations (0–32)**. The first group includes the refined populations recovered from unresolved-cell reclustering, followed by the original CCA populations that were already resolved.

The final number-to-cell-type map used in the CRC UMAP is:

| Final | Annotation |
|---:|---|
| 0 | LGR5/OLFM4+ epithelial stem/progenitor-like cells |
| 1 | LGR5+ absorptive progenitor-like epithelial cells |
| 2 | Cycling epithelial cells |
| 3 | DKK1+/MUC5AC+ stress-associated epithelial cells |
| 4 | Activated antigen-presenting B cells |
| 5 | MHC-II-high B cells |
| 6 | B-cell / plasma-cell population |
| 7 | Plasma/plasmablast-like secretory B cells |
| 8 | Plasma cells / plasmablast-like cells |
| 9 | Cycling B-lineage / plasmablast-like cells |
| 10 | GZMM+ regulatory T cells |
| 11 | Activated regulatory T cells |
| 12 | Gamma-delta T / NK-like cytotoxic lymphocytes |
| 13 | Cycling gamma-delta T/NK-like cytotoxic lymphocytes |
| 14 | ACKR1+/MADCAM1+ venous endothelial cells |
| 15 | RGS5+/PDGFRB+ mural/pericyte cells |
| 16 | COL2A1/COL11A2+ mesenchymal stromal cells |
| 17 | INHBA+/BMP2+ vascular/perivascular stromal-like cells |
| 18 | Activated CAF-like fibroblasts |
| 19 | Cytotoxic CD8 effector/effector-memory T cells |
| 20 | Dendritic cells |
| 21 | Epithelial cells |
| 22 | Enterocyte-like tumor epithelium |
| 23 | Goblet/secretory epithelium with Paneth-like features |
| 24 | Inflammatory squamous-like tumor epithelium |
| 25 | MAIT / Th17-like innate-like T cells |
| 26 | Mast cells |
| 27 | Naive/central-memory CD4 T cells |
| 28 | NK / gamma-delta T cytotoxic lymphocytes |
| 29 | Proliferating germinal-center-like B cells |
| 30 | Proliferating myeloid cells |
| 31 | S100A8/A9/A12+ inflammatory myeloid cells |
| 32 | TREM2+ C1Q+ tumor-associated macrophages |

A machine-readable version is saved with the final atlas as:

```text
final_annotated_conditions/annotated_CRC_annotation_key.csv
```

The numbered UMAP is saved as:

```text
final_annotated_conditions/annotated_CRC_UMAP.png
final_annotated_conditions/annotated_CRC_UMAP.pdf
```

The final annotated Seurat object is:

```text
final_annotated_conditions/annotated_CRC.rds
```

The final numbering is used for visualization and legend organization; the original CCA cluster identity remains preserved in the Seurat metadata.

## Final Seurat metadata

The refined CRC object preserves both the original and refined analysis.

```text
cca_snn_res.0.5
    Original global CRC CCA cluster

CCA_annotation_final
    Original first-pass CCA annotation

unresolved_recluster_id
    Local recluster 0-18 for cells that required refinement
    NA for cells that did not require unresolved-cell reclustering

unresolved_recluster_annotation
    Biological identity assigned from the local unresolved recluster

CCA_annotation_refined
    Final biological annotation used in the refined full CRC atlas

annotation_source
    "Original CCA annotation"
    or
    "Unresolved-cell reclustering"

final_cluster_number
    Numeric label used on the final CRC UMAP

final_cluster_label
    Number + biological annotation used in the UMAP legend
```

## Why the original unresolved cluster IDs are not final labels

The original unresolved CCA populations were heterogeneous.

A single original unresolved cluster could redistribute into multiple local reclusters after refinement. It would therefore be incorrect to simply rename every cell from an old unresolved CCA cluster with one new label.

The refinement is performed cell-by-cell:

```text
original unresolved cell barcode
        ->
local unresolved recluster assignment
        ->
refined biological annotation
        ->
mapped back into the full CRC object
```

This preserves the original CCA identity while allowing the final atlas annotation to reflect the finer biological structure recovered during reclustering.

## Reproducibility principle

The repository intentionally preserves the distinction between:

- original CCA integration and clustering,
- donor-level formal pseudobulk DE,
- numerical fallback models,
- unresolved-cell CCA refinement,
- annotation-oriented cell-level marker testing,
- deeper parent-9/12 subcluster analysis,
- donor-recovery provenance,
- and final metadata integration.

Large Seurat objects, raw-count matrices, pseudobulk caches, fitted model objects, and generated DE result directories remain on OSCER rather than being committed to GitHub.

The repository contains the reproducible scripts, configuration maps, and documentation required to trace the final CRC annotations back to the analysis stage that produced them.

---

## Author

**Olajumoke B. Oladapo**  
Ph.D. Biomedical Engineering  
University of Oklahoma  
Moussa Lab
