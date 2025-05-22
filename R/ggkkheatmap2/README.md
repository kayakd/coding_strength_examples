# ggkkheatmap2

`ggkkheatmap2` is a flexible and publication-ready heatmap visualization function built on `ggplot2`, developed by **Koray Doğan Kaya**. It enables seamless integration of hierarchical clustering, sample and gene group annotations, and intelligent display controls for large gene expression matrices.

---

## Features

- Clustered heatmaps with optional **row** and **column dendrograms**
- **Gene and sample grouping** with automatic cluster boundary lines
- Smart label display: hides gene names when too many rows
- Supports **z-score**, **logCPM**, and **viridis** color palettes
- Custom row and column order via `hclust` or `dendrogram` objects
- **Colored labels** for both row and column clusters
- Clean layout using `ggplot2`, `patchwork`, and `gridExtra`

---

## Installation

```r
# Ensure required packages are installed
install.packages(c("ggplot2", "ggdendro", "RColorBrewer", "scales", "reshape2", "patchwork", "gtable", "gridExtra", "viridis"))
```
# Simulate test data
```r
set.seed(42)
mat <- matrix(rnorm(20 * 10), nrow = 20, dimnames = list(paste0("Gene", 1:20), paste0("Sample", 1:10)))

# Create row and column clusters
row_dist <- dist(mat)
col_dist <- dist(t(mat))
row_clust <- hclust(row_dist)
col_clust <- hclust(col_dist)

# Sample grouping: 3 clusters by hand
col_groups <- cutree(col_clust, k = 3)
row_groups <- cutree(row_clust, k = 4)

# Visualize
ggkkheatmap2(mat,
             row_clust = row_clust,
             col_clust = col_clust,
             row_groups = row_groups,
             col_groups = col_groups,
             scale_values = TRUE,
             palette = "zscore",
             base_size = 11,
             border_linewidth = 0.01)
```
