# CRC Atlas Differential Expression and Cell-Type Annotation

This repository contains the reproducible **differential-expression, unresolved-population refinement, and final cell-type annotation workflows** for the condition-specific colorectal cancer single-cell atlas.

The repository covers four biological conditions:

- **Adenoma/Polyp**
- **Colorectal cancer (CRC)**
- **Normal**
- **Tumor-adjacent normal (TAN)**

The condition-specific workflows begin from CCA-integrated Seurat objects and preserve biological replication during formal differential-expression testing. Each condition is documented separately because the unresolved populations, available biological replication, numerical behavior of the statistical models, and annotation-refinement strategy were not identical across conditions.

---

## Repository structure

```text
crc-atlas-de-annotation/
├── README.md
├── adenoma/
│   ├── README.md
│   ├── config/
│   └── scripts/
├── crc/
│   ├── README.md
│   ├── config/
│   └── scripts/
├── normal/
│   ├── README.md
│   ├── config/
│   └── scripts/
└── tan/
    ├── README.md
    ├── config/
    └── scripts/
```

Each condition directory contains its own detailed README, configuration/annotation maps, and analysis scripts.

---

## Atlas conditions

| Condition | Approximate cells | Refinement strategy |
|---|---:|---|
| Adenoma/Polyp | 306,396 | Unresolved-cell CCA reclustering followed by DE/marker-based annotation |
| CRC | 702,029 | Unresolved-cell CCA reclustering followed by replicate-aware DE, fallback models where required, and annotation refinement |
| Normal | 151,684 | Direct DE-based refinement of unresolved original CCA clusters |
| TAN | 287,894 | Direct DE-based refinement of unresolved original CCA clusters; one source-restricted cluster explicitly excluded |

The original condition-level clustering variable is:

```text
cca_snn_res.0.5
```

When unresolved-cell reclustering was performed, the new cluster IDs are stored separately and are **not treated as original CCA cluster IDs**.

---

# Analysis framework

## 1. Condition-specific CCA integration

The annotation workflows begin from condition-specific CCA-integrated Seurat objects.

The repository does not treat cells from different studies as independent biological replicates during formal differential-expression testing.

---

## 2. Raw-count pseudobulk differential expression

Where raw RNA counts and biological replicate identifiers are available, cells are aggregated into pseudobulk profiles at the biological-replicate level.

For a selected population:

```text
target = selected cluster
rest   = all remaining comparison clusters
```

Eligible biological replicate units generally require:

```text
target cells >= 20
rest cells   >= 20
paired biological replicates >= 3
```

The primary original-CCA model uses a paired design of the form:

```r
~ replicate + group
```

with positive fold change representing higher expression in the target population relative to the rest.

Typical reporting thresholds are:

```text
adjusted P value / FDR < 0.05
|log2 fold change| >= 0.25
strong upregulation: log2 fold change >= 1
```

---

## 3. Processed CELLxGENE matrices

The atlas contains datasets obtained through CELLxGENE in addition to GEO-derived datasets.

Processed/normalized CELLxGENE matrices are not treated as raw integer-count matrices for formal count-based pseudobulk DE.

In particular:

```text
CELLxGENE_Chen_2021
CELLxGENE_Moorman2025
```

are excluded from formal raw-count DE when only processed expression matrices are available.

These cells can still contribute to the integrated atlas and to biological interpretation where appropriate.

---

# Condition-specific workflows

## Adenoma/Polyp

Directory:

```text
adenoma/
```

The Adenoma/Polyp workflow begins with the original CCA clusters, identifies unresolved populations, reclusters the unresolved cells, performs differential-expression/marker refinement, and maps the refined annotations back to the full Adenoma atlas.

Original unresolved CCA clusters were separated from the already-resolved populations before secondary clustering.

The final workflow produces a full annotated Adenoma/Polyp object and numbered UMAP.

See:

```text
adenoma/README.md
```

for the complete analysis history, scripts, and final annotation map.

---

## CRC

Directory:

```text
crc/
```

CRC contains the largest condition-specific atlas and required the most extensive unresolved-population refinement.

The original unresolved CRC CCA clusters were:

```text
3, 5, 8, 9, 10, 15, 16, 18, 19, 20, 25, 26, 27, 28
```

These cells were extracted and reclustered using CCA, producing local unresolved-recluster IDs that are intentionally kept distinct from the original global CCA IDs.

Formal replicate-aware DE was performed using the most appropriate stable model for each population. The workflow includes:

- DESeq2 pseudobulk analysis,
- edgeR fallback analysis,
- low-memory limma-voom fallback analysis,
- annotation-oriented Wilcoxon analysis where replicate-level formal DE was not supported,
- additional refinement of selected heterogeneous B-cell/plasma-cell populations,
- final merge-back into the full CRC atlas.

The final CRC atlas contains **33 numbered populations (0–32)**.

See:

```text
crc/README.md
```

for the complete analysis history, statistical branches, annotation tables, and final population map.

---

## Normal

Directory:

```text
normal/
```

The first-pass Normal annotation left three original CCA clusters unresolved:

```text
0
10
21
```

Unlike CRC and Adenoma, these populations were refined directly at the original CCA-cluster level rather than undergoing a second unresolved-cell reclustering step.

Final refined annotations are:

| Original CCA cluster | Final annotation |
|---:|---|
| 0 | Paneth-like / LGR5+ crypt epithelial cells |
| 10 | Plasma cells / plasmablasts |
| 21 | Cycling secretory/crypt epithelial cells |

The final Normal object retains the original CCA cluster identity and adds the refined annotation.

See:

```text
normal/README.md
```

for the complete annotation map and scripts.

---

## Tumor-adjacent normal (TAN)

Directory:

```text
tan/
```

The first-pass TAN annotation left three original CCA clusters unresolved:

```text
17
26
34
```

Clusters 17 and 26 were refined using their differential-expression profiles:

| Original CCA cluster | Final annotation |
|---:|---|
| 17 | Cycling epithelial/B-lineage mixed population |
| 26 | POU2F3+/TRPM5+ tuft-like epithelial cells |

### TAN cluster 34

Cluster 34 is restricted to:

```text
CELLxGENE_Moorman2025
```

and could not be confidently assigned a biological identity from the available evidence.

Rather than forcing an unsupported annotation, the workflow records it as:

```text
Excluded - CELLxGENE_Moorman2025
```

with an explicit exclusion reason.

Cluster 34 is retained in the audit information but excluded from the final downstream TAN analysis object.

See:

```text
tan/README.md
```

for the complete TAN annotation map and exclusion workflow.

---

# Statistical methods used across the repository

The exact model used for a population depends on the available biological replication and numerical stability.

Methods represented in the repository include:

### DESeq2

Used for donor-aware raw-count pseudobulk testing where the model could be fit reliably.

Typical design:

```r
~ replicate + group
```

### edgeR

Used as a replicate-aware count-based fallback for selected comparisons when the original DESeq2 route was numerically unstable.

### limma-voom

Used for selected paired pseudobulk comparisons where a low-memory or numerically stable linear-model workflow was required.

The general workflow includes:

```r
edgeR::filterByExpr()
edgeR::calcNormFactors(method = "TMM")
edgeR::voomLmFit()
limma::eBayes()
```

### Wilcoxon marker analysis

Used only for annotation support where formal replicate-level differential expression was not justified or where the purpose was specifically marker discovery rather than population-level inference.

Wilcoxon results are not presented as a replacement for replicate-aware pseudobulk inference.

---

# Reproducibility principles

The repository follows several important rules:

1. **Biological replication is preserved.**  
   Individual cells are not treated as independent biological replicates for formal pseudobulk DE.

2. **Raw-count testing uses raw counts.**  
   Processed expression matrices are not silently substituted for raw integer-count input.

3. **Original and reclustered IDs remain distinct.**  
   Original `cca_snn_res.0.5` cluster IDs are not confused with local unresolved-recluster IDs.

4. **Unresolved populations are not forced into unsupported identities.**  
   Populations without sufficient biological evidence can remain unresolved or be explicitly excluded.

5. **Fallback statistical models are documented.**  
   A successful fallback analysis does not overwrite the history of the primary model attempt.

6. **Large data objects remain outside GitHub.**  
   Seurat objects, pseudobulk caches, large intermediate matrices, model objects, and full result directories remain on OSCER storage.

7. **Analysis provenance is retained.**  
   Scripts save summaries, annotation maps, diagnostic information, and `sessionInfo()` where applicable.

---

# Atlas storage

The analysis was developed using the atlas root:

```text
/ourdisk/hpc/moussalab/jumi/dont_archive/atlas
```

Many scripts allow the root to be overridden with:

```bash
export CRC_ATLAS_ROOT=/path/to/atlas
```

Large analysis outputs are intentionally not version-controlled.

---

# Final annotated-condition outputs

Final condition-level objects and UMAPs are stored under:

```text
<atlas_root>/final_annotated_conditions/
```

with outputs including:

```text
annotated_adenoma.rds
annotated_adenoma_UMAP.png
annotated_adenoma_UMAP.pdf

annotated_CRC.rds
annotated_CRC_UMAP.png
annotated_CRC_UMAP.pdf

annotated_normal.rds
annotated_normal_UMAP.png
annotated_normal_UMAP.pdf

annotated_TAN.rds
annotated_TAN_UMAP.png
annotated_TAN_UMAP.pdf
```

Condition-specific annotation keys and audit tables are also generated where applicable.

---

# Recommended use of this repository

This repository is intended to document and reproduce the **DE and annotation stage** of the colorectal cancer single-cell atlas.

For a particular condition, use the corresponding directory:

```text
adenoma/
crc/
normal/
tan/
```

and follow the execution order documented in that condition's README.

Because the refinement strategy differs by condition, the individual READMEs should be considered the authoritative description of the exact analysis history for each condition.

---

## Author

**Olajumoke B. Oladapo**  
Ph.D. Biomedical Engineering  
University of Oklahoma  
Moussa Lab
