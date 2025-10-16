## Single-Cell RNA-Seq Workshop
We are facing a Medical Mystery in this course. two patients get the same virus. David recovers in 5 days. Emily ends up in the ICU. Their standard blood tests look nearly identical. What's different?
Spoiler: Bulk RNA-seq can't tell us, but single-cell can.

### What You'll Learn
We'll walk through the complete single-cell analysis pipeline:

* Quality Control: Filter dead cells, empty droplets, and doublets
* Normalization: Use SCTransform to remove technical noise
* Dimensionality Reduction: PCA and UMAP to visualize thousands of genes in 2D
* Clustering: Group similar cells together
* Cell Type Annotation: Figure out what each cluster actually is (T cells? Monocytes?)
* Differential Expression: Compare David vs Emily the right way (pseudo-bulk, not cell-level!)

### Prerequisites
Basic R knowledge. Install these packages:
```r
install.packages(c("Seurat", "tidyverse", "patchwork"))
BiocManager::install(c("scDblFinder", "SingleR", "celldex"))
```
### Questions or concerns? 
Open an issue and we'll help you out.
