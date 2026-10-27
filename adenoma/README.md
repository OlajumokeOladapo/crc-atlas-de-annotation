# Adenoma/Polyp: Full DE and Annotation Workflow

This folder documents the complete Adenoma/Polyp analysis used in the CRC single-cell atlas, from the original CCA clusters through refinement of previously unresolved cells and final integration into a 29-population annotated atlas.

## Analysis flow

```text
Full Adenoma CCA object
306,396 cells
cca_snn_res.0.5 clusters 0-24
        |
        v
Raw-count donor-level pseudobulk DESeq2
cluster vs all other clusters
        |
        v
Original CCA annotation
14 resolved populations
11 unresolved CCA clusters
        |
        v
Extract only unresolved cells
        |
        v
CCA reclustering
new reclusters 0-14
        |
        v
Strict-Adenoma DE / marker refinement
        |
        v
Final annotations for reclusters 0-14
        |
        v
Map refined cells back to full object by cell barcode
        |
        v
Final full Adenoma atlas
clusters 1-29
```

## Scripts

Run the scripts in order.

### 01 — Original all-cluster pseudobulk DESeq2

```text
scripts/01_adenoma_CCA_DESeq2_pseudobulk.R
```

Starts from:

```text
Adenoma_CCA_integrated_WITH_CHEN.rds
```

and tests every original `cca_snn_res.0.5` cluster against all remaining CCA clusters using donor-level raw-count pseudobulks.

Original settings:

```text
target cells per donor >= 20
rest cells per donor   >= 20
eligible paired donors >= 3
gene total count       >= 10
design                 ~ replicate + group
DESeq2 test            Wald
size-factor method     poscounts
padj cutoff            0.05
|log2FC| reporting     >= 0.25
```

`CELLxGENE_Chen_2021` remains part of the integrated atlas but is excluded from raw-count DE because the available matrix is processed/normalized rather than raw integer counts.

### 02 — Apply original CCA annotations

```text
scripts/02_apply_original_adenoma_annotations.R
```

Creates:

```text
CCA_annotation_final
```

The original CCA clusters classified as unresolved were:

```text
0, 1, 2, 4, 5, 7, 15, 17, 18, 23, 24
```

Cluster 17 was initially called:

```text
Ribosome/OXPHOS-high unresolved
```

The complete original map is stored in:

```text
config/original_CCA_annotation_map.csv
```

### 03 — Reclustering of original unresolved cells

```text
scripts/03_recluster_original_unresolved_cells_CCA.R
```

This script extracts the cells from the 11 original unresolved CCA clusters and reclusters only those cells.

Important implementation choices:

- uses the retained normalized RNA expression,
- uses the existing condition-level variable features,
- does **not** rerun VST/`FindVariableFeatures()` on synthetic integration count layers,
- uses CCA integration for consistency with the parent Adenoma atlas,
- uses PCs 1–20,
- clusters at resolution 0.5,
- stores the new cluster in `unresolved_cluster_res0.5`,
- saves the reclustered UMAP as `umap.cca.unresolved`.

### 04 — DE and annotation of the new unresolved-cell reclusters

```text
scripts/04_unresolved_DE_and_annotation.R
```

This downstream stage uses the new reclusters 0–14.

For formal raw-count pseudobulk DE, the later strict-Adenoma rule is applied to GSE261388:

```text
P1   inflammatory polyp -> excluded
P12  inflammatory polyp -> excluded
P3   adenomatous polyp  -> retained
P4   adenomatous polyp  -> retained
```

The primary formal analysis uses DESeq2 LRT when the recluster has adequate biological replication.

Special cases:

- **Recluster 1**: the formal DE models did not provide a stable positive marker program. Annotation was supported using descriptive `log2((CPM + 1)/(rest CPM + 1))` enrichment and a paired biological-sample Wilcoxon test on `log2(CPM + 1)`.
- **Reclusters 9 and 14**: these were overwhelmingly contributed by `GSM6061707`; cell-level Seurat Wilcoxon was therefore used only for marker discovery/annotation, not as independent biological-replicate inference.

Final recluster annotations:

| Recluster | Annotation |
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

### 05 — Integrate refined cells into the full atlas and renumber

```text
scripts/05_merge_refined_annotations_final_numbering.R
```

This script does **not** rerun clustering or UMAP.

Instead, it:

1. loads the full annotated Adenoma object,
2. loads the refined unresolved-cell object,
3. matches cells by exact cell barcode,
4. preserves original annotations for previously resolved cells,
5. replaces only the former unresolved cells with their refined recluster annotation,
6. preserves the original CCA cluster in `cca_snn_res.0.5`,
7. stores the new recluster in `unresolved_recluster_id`,
8. assigns the final cluster numbers 1–29,
9. generates the final full-atlas numbered UMAP,
10. saves the final Seurat object.

## Final cluster numbering

The first 14 final numbers correspond to original CCA populations that were already resolved. The old unresolved CCA clusters do not retain final cluster numbers; their individual cells are redistributed according to the new reclustering. New reclusters 0–14 continue the numbering as final clusters 15–29.

| Final | Source | Annotation |
|---:|---|---|
| 1 | CCA 3 | WNT-high stem/progenitor epithelial |
| 2 | CCA 6 | Myofibroblast-like stroma |
| 3 | CCA 8 | GPNMB+ / MMP12+ macrophages |
| 4 | CCA 9 | Inflammatory/regenerative epithelial |
| 5 | CCA 10 | Vascular mural/endothelial-like cells |
| 6 | CCA 11 | Metabolic absorptive epithelial |
| 7 | CCA 12 | Activated antigen-presenting B cells |
| 8 | CCA 13 | Activated T / gamma-delta T–NK-like lymphocytes |
| 9 | CCA 14 | Enteroendocrine cells |
| 10 | CCA 16 | Enteric glia / neural-like cells |
| 11 | CCA 19 | Goblet-like secretory epithelial |
| 12 | CCA 20 | Mature absorptive enterocyte/colonocyte |
| 13 | CCA 21 | Endothelial cells |
| 14 | CCA 22 | CD163+ tissue-resident-like macrophages |
| 15 | Recl. 0 | LGR5+ stem/progenitor-like epithelial cells |
| 16 | Recl. 1 | CA1/CA2+ differentiated absorptive colonocyte-like cells |
| 17 | Recl. 2 | Goblet/secretory epithelial cells |
| 18 | Recl. 3 | Cycling/proliferating cells |
| 19 | Recl. 4 | Inflammatory/regenerative crypt epithelial cells |
| 20 | Recl. 5 | Mature absorptive colonocytes |
| 21 | Recl. 6 | Activated/plasmablast-like B cells |
| 22 | Recl. 7 | Smooth muscle/myofibroblast-like stromal cells |
| 23 | Recl. 8 | POU2F3+ tuft-like/tuft-precursor cells |
| 24 | Recl. 9 | OLFM4+ intestinal stem/progenitor epithelial cells |
| 25 | Recl. 10 | WNT-high intestinal stem/progenitor epithelial cells |
| 26 | Recl. 11 | Tuft cells |
| 27 | Recl. 12 | BEST4+ OTOP2+ absorptive epithelial cells |
| 28 | Recl. 13 | G2/M-high proliferating cells |
| 29 | Recl. 14 | Cycling OLFM4+ stem/progenitor epithelial cells |

A machine-readable version is provided in:

```text
config/final_cluster_number_map.csv
```

## Final Seurat metadata

The final object preserves both the original and refined analysis:

```text
cca_snn_res.0.5
    Original CCA cluster 0-24

CCA_annotation_final
    Original first-pass annotation

unresolved_recluster_id
    New 0-14 recluster ID for formerly unresolved cells
    NA for cells that did not require reclustering

CCA_annotation_refined
    Final biological annotation

annotation_source
    "Original CCA annotation"
    or
    "Unresolved-cell reclustering"

final_cluster_number
    Final atlas number 1-29

final_cluster_source_id
    e.g. CCA_22 or RECLUSTER_9

final_cluster_label
    e.g. "24: OLFM4+ intestinal stem/progenitor epithelial cells"
```

The final saved object is:

```text
Adenoma_polyp_CCA_FINAL_REFINED_RENUMBERED.rds
```

## Why the old unresolved cluster numbers disappear

The unresolved clusters were heterogeneous. A single old cluster could contribute cells to many different new reclusters. Therefore it would be incorrect to simply rename an entire old unresolved CCA cluster.

The refinement is done cell-by-cell:

```text
old unresolved cell barcode
        ->
new recluster assignment
        ->
new refined annotation
        ->
mapped back into the full Adenoma object
```

For example, old CCA cluster 17 (`Ribosome/OXPHOS-high unresolved`) was strongly dominated after refinement by the `OLFM4+ intestinal stem/progenitor epithelial` recluster, while smaller fractions entered other epithelial states. The final atlas therefore uses each cell's refined recluster identity rather than treating old CCA cluster 17 as one final population.

## Reproducibility principle

The repository intentionally preserves the distinction between:

- original atlas integration/clustering,
- donor-level formal DE,
- unresolved-cell refinement,
- annotation-oriented marker testing,
- and final metadata integration.

This makes it possible to trace every final Adenoma cell back to both its original CCA cluster and, when applicable, its refined recluster.
