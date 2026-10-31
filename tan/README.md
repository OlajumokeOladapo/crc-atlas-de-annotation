# TAN: Full DE and Annotation Workflow

This folder documents the tumor-adjacent-normal (TAN) analysis in the colorectal cancer single-cell atlas, from original CCA clusters through direct refinement of unresolved populations and explicit exclusion of an unannotatable CELLxGENE_Moorman2025-only cluster.

## Analysis flow

```text
Full TAN CCA object
~287,894 cells
        |
        v
Original CCA donor-level pseudobulk DESeq2
        |
        v
First-pass annotation
Unresolved CCA clusters: 17, 26, 34
        |
        v
Direct paired limma-voom refinement
17 and 26
        |
        +------------------------------+
        |                              |
        v                              v
17 resolved                       34 excluded
26 resolved                       CELLxGENE_Moorman2025-only
        |
        v
annotated_TAN.rds + numbered UMAP
```

## Scripts

1. `scripts/01_TAN_original_CCA_pseudobulk_DESeq2.R`
2. `scripts/02_TAN_unresolved_CCA_limma_voom.R`
3. `scripts/03_apply_final_TAN_annotations_and_exclude34.R`
4. `scripts/04_TAN_final_numbered_UMAP_and_export.R`

## Unresolved-cluster refinement

Final refined labels:

| CCA | Final annotation |
|---:|---|
| 17 | Cycling epithelial/B-lineage mixed population |
| 26 | POU2F3+/TRPM5+ tuft-like epithelial cells |
| 34 | Excluded - CELLxGENE_Moorman2025 |

Cluster 34 is excluded because it is restricted to `CELLxGENE_Moorman2025` and cannot be confidently annotated from the available evidence. It is retained in the full audit object with an exclusion reason but removed from the final downstream TAN object.

## Complete final TAN annotation map

The full machine-readable map is `config/TAN_final_CCA_annotation_map.csv`. It preserves all original resolved CCA annotations, the refined labels for clusters 17 and 26, and the explicit exclusion status for cluster 34.

Final outputs are written under `<atlas_root>/final_annotated_conditions/` as `annotated_TAN.rds`, `annotated_TAN_UMAP.png`, `annotated_TAN_UMAP.pdf`, `annotated_TAN_annotation_key.csv`, `TAN_annotation_summary.csv`, and `TAN_cluster34_excluded_cells.csv`.

---

## Author

**Olajumoke B. Oladapo**  
Ph.D. Biomedical Engineering  
University of Oklahoma  
Moussa Lab

