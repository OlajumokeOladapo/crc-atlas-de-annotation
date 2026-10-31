# Normal: Full DE and Annotation Workflow

This folder documents the Normal-condition analysis in the colorectal cancer single-cell atlas, from the original CCA clusters through direct refinement of the remaining unresolved CCA populations and final export.

## Analysis flow

```text
Full Normal CCA object
~151,684 cells
        |
        v
Original CCA donor-level pseudobulk DESeq2
        |
        v
First-pass annotation
Unresolved CCA clusters: 0, 10, 21
        |
        v
Direct paired limma-voom refinement
        |
        v
Final annotations
        |
        v
annotated_normal.rds + numbered UMAP
```

## Scripts

1. `scripts/01_Normal_original_CCA_pseudobulk_DESeq2.R`
2. `scripts/02_Normal_unresolved_CCA_limma_voom.R`
3. `scripts/03_apply_final_Normal_annotations.R`
4. `scripts/04_Normal_final_numbered_UMAP_and_export.R`

## Unresolved-cluster refinement

The original first-pass Normal annotation left CCA clusters `0`, `10`, and `21` unresolved. These were directly tested against all other original Normal CCA clusters using paired donor-level pseudobulk limma-voom.

Final refined labels:

| CCA | Final annotation |
|---:|---|
| 0 | Paneth-like / LGR5+ crypt epithelial cells |
| 10 | Plasma cells / plasmablasts |
| 21 | Cycling secretory/crypt epithelial cells |

No second unresolved-cell reclustering is claimed for Normal.

## Complete final Normal annotation map

| CCA | Annotation |
|---:|---|
| 0 | Paneth-like / LGR5+ crypt epithelial cells |
| 1 | Differentiated surface absorptive epithelium, metallothionein-high |
| 2 | CD4 naive/memory T cells, ICOS+ / CTLA4+ |
| 3 | Mature distal-type colonocytes |
| 4 | Goblet cells |
| 5 | Naive/mature B cells |
| 6 | HOXB13+ OXPHOS-high colonic epithelium |
| 7 | Cycling transit-amplifying epithelial cells |
| 8 | Intestinal stem/progenitor cells |
| 9 | Mature absorptive enterocytes/colonocytes |
| 10 | Plasma cells / plasmablasts |
| 11 | Subepithelial/telocyte-like fibroblasts |
| 12 | Cytotoxic CD8 / gamma-delta T / NK-like lymphocytes |
| 13 | Smooth muscle cells |
| 14 | CFD+ ADH1B+ SCARA5+ fibroblasts |
| 15 | CD163+ tissue-resident macrophages |
| 16 | BEST4+ OTOP2+ absorptive epithelium |
| 17 | Blood vascular endothelial cells, mixed subtypes |
| 18 | Mixed enteric glia/neural and mast-cell cluster |
| 19 | Tuft cells |
| 20 | Enteroendocrine cells |
| 21 | Cycling secretory/crypt epithelial cells |
| 22 | Lymphatic endothelial cells |

The machine-readable map is `config/Normal_final_CCA_annotation_map.csv`.

Final outputs are written under `<atlas_root>/final_annotated_conditions/` as `annotated_normal.rds`, `annotated_normal_UMAP.png`, `annotated_normal_UMAP.pdf`, and `annotated_normal_annotation_key.csv`.
