#' Plot PCA from log-transformed expression matrix
#'
#' Created by Koray Doğan Kaya
#'
#' This function computes PCA from a log2-transformed expression matrix (genes × samples)
#' and generates a PCA plot using sample group metadata.
#'
#' @param x            A numeric matrix of log-transformed expression values (genes × samples).
#' @param sample_meta  A data frame containing sample metadata. Row names must match colnames(x).
#' @param group_col    The column name in sample_meta to color by (default: "Group").
#' @param title        Plot title (default: "").
#'
#' @return A ggplot2 object visualizing the first two principal components.
#'
plot_pca <- function(x, sample_meta, group_col = "Group", title = "") {
  # Dependencies
  if (!require(ggplot2, quietly = TRUE)) stop("Please install 'ggplot2'.")
  if (!require(ggrepel, quietly = TRUE)) stop("Please install 'ggrepel'.")

  # Check input consistency
  if (!is.matrix(x)) stop("Input 'x' must be a numeric matrix.")
  if (!all(colnames(x) %in% rownames(sample_meta))) stop("Column names of x must match rownames of sample_meta.")
  if (!group_col %in% colnames(sample_meta)) stop(paste("Column", group_col, "not found in sample_meta."))

  # PCA
  pca <- prcomp(t(x), scale. = TRUE)

  # Prepare data
  pca_data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    Group = sample_meta[colnames(x), group_col],
    Sample = colnames(x)
  )

  # Plot
  ggplot(pca_data, aes(x = PC1, y = PC2, color = Group, label = Sample)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 10) +
    labs(
      title = title,
      x = paste0("PC1 (", round(summary(pca)$importance[2, 1] * 100, 1), "% variance)"),
      y = paste0("PC2 (", round(summary(pca)$importance[2, 2] * 100, 1), "% variance)")
    ) +
    theme_minimal()
}

