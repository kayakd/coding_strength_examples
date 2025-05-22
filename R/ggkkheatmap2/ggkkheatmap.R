#' Enhanced Heatmap with Dendrograms and Group Annotations
#'
#' This function generates a ggplot-based heatmap with optional row and column dendrograms, gene/sample clustering,
#' group annotations, and scaling. It is especially useful for visualizing clustered gene expression matrices.
#'
#' @param mat A numeric matrix with genes as rows and samples as columns.
#' @param row_clust Optional hierarchical clustering (`hclust` or `dendrogram`) for rows (genes).
#' @param col_clust Optional hierarchical clustering (`hclust` or `dendrogram`) for columns (samples).
#' @param row_groups Optional vector indicating group assignment for each row. Used for color annotations.
#' @param col_groups Optional vector indicating group assignment for each column. Used for color annotations.
#' @param scale_values Logical; whether to z-score scale each row.
#' @param palette Name of the color palette to use. Supported: "zscore", "logcpm", or any viridis-compatible.
#' @param base_size Base font size for the plot.
#' @param border_linewidth Line width of border lines drawn between cluster groups.
#'
#' @return A ggplot2 heatmap grob object, optionally combined with dendrograms.
#'
#' @examples
#' mat <- matrix(rnorm(100), nrow = 10)
#' rownames(mat) <- paste0("Gene", 1:10)
#' colnames(mat) <- paste0("Sample", 1:10)
#' dend <- hclust(dist(mat))
#' ggkkheatmap2(mat, row_clust = dend, scale_values = TRUE)
#'
#' @author Created by Koray Doğan Kaya


 ggkkheatmap2 <- function(
    mat,
    row_clust = NULL,
    col_clust = NULL,
    row_groups = NULL,
    col_groups = NULL,
    scale_values = FALSE,
    palette = "zscore",
    base_size = 11,
    border_linewidth = 0.3
) {
  require(ggplot2)
  require(scales)
  require(ggdendro)
  require(reshape2)
  require(patchwork)
  require(RColorBrewer)
  require(grid)
  require(ggplotify)
  require(gtable)
  require(gridExtra)

  if (!is.matrix(mat)) stop("Input 'mat' must be a matrix")

  if (scale_values) {
    mat <- t(scale(t(mat)))
  }

  if (is.null(rownames(mat))) {
    rownames(mat) <- paste0("Gene_", seq_len(nrow(mat)))
  }

  if (!is.null(row_clust)) {
    if (inherits(row_clust, "hclust")) {
      dend <- as.dendrogram(row_clust)
    } else if (inherits(row_clust, "dendrogram")) {
      dend <- row_clust
    } else {
      stop("row_clust must be an 'hclust' or 'dendrogram' object")
    }

    row_order <- order.dendrogram(dend)

    if (is.character(row_order)) {
      if (!all(row_order %in% rownames(mat))) stop("Character-based dendrogram order does not match rownames of mat")
      mat <- mat[row_order, , drop = FALSE]
      if (!is.null(row_groups)) row_groups <- row_groups[match(row_order, rownames(mat))]
    } else {
      if (length(row_order) != nrow(mat)) stop("row_clust order does not match number of matrix rows")
      mat <- mat[row_order, , drop = FALSE]
      if (!is.null(row_groups)) row_groups <- row_groups[row_order]
    }
  }

  if (!is.null(col_clust)) {
    if (inherits(col_clust, "hclust")) {
      col_dend <- as.dendrogram(col_clust)
    } else if (inherits(col_clust, "dendrogram")) {
      col_dend <- col_clust
    } else {
      stop("col_clust must be an 'hclust' or 'dendrogram' object")
    }

    col_order <- order.dendrogram(col_dend)

    if (is.character(col_order)) {
      if (!all(col_order %in% colnames(mat))) stop("Character-based dendrogram order does not match colnames of mat")
      mat <- mat[, col_order, drop = FALSE]
      if (!is.null(col_groups)) col_groups <- col_groups[match(col_order, colnames(mat))]
    } else {
      if (length(col_order) != ncol(mat)) stop("col_clust order does not match number of matrix columns")
      mat <- mat[, col_order, drop = FALSE]
      if (!is.null(col_groups)) col_groups <- col_groups[col_order]
    }
  }

  mat <- mat[rowSums(is.na(mat)) < ncol(mat), , drop = FALSE]

  df_long <- reshape2::melt(mat)
  colnames(df_long) <- c("Gene", "Sample", "value")
  df_long$Gene <- factor(df_long$Gene, levels = rownames(mat), ordered = TRUE)
  df_long$Sample <- factor(df_long$Sample, levels = colnames(mat), ordered = TRUE)

  if (!is.null(row_groups)) df_long$row_group <- rep(row_groups, times = ncol(mat))
  if (!is.null(col_groups)) df_long$col_group <- rep(col_groups, each = nrow(mat))

  fill_scale <- if (palette == "logcpm") {
    cols1 <- rev(brewer.pal(n = 6, name = "BuPu"))
    cols2 <- brewer.pal(n = 9, name = "YlOrRd")
    cols <- colorRampPalette(c(cols1, cols2))(20)
    scale_fill_gradientn(colors = cols, name = palette)
  } else if (palette == "zscore") {
    scale_fill_gradientn(colors = c("blue", "gray", "red"), limits = c(-4, 4), name = palette)
  } else {
    scale_fill_gradientn(colors = viridis::viridis(100), name = palette)
  }

  show_ylab <- nrow(mat) <= 100
  gene_lab_df <- data.frame(Gene = rownames(mat), y = seq_len(nrow(mat)))
  gene_lab_df$color <- if (!is.null(row_groups)) {
    hue_pal()(length(unique(row_groups)))[as.numeric(factor(row_groups))]
  } else "black"

  sample_lab_df <- data.frame(Sample = colnames(mat))
  sample_lab_df$color <- if (!is.null(col_groups)) {
    hue_pal()(length(unique(col_groups)))[as.numeric(factor(col_groups))]
  } else "black"

  vlines <- if (!is.null(col_groups)) which(diff(as.numeric(factor(col_groups))) != 0) + 0.5 else NULL
  hlines <- if (!is.null(row_groups)) which(diff(as.numeric(factor(row_groups))) != 0) + 0.5 else NULL

  p_heat <- ggplot(df_long, aes(x = Sample, y = Gene, fill = value)) +
    scale_y_discrete(position = "right") +
    geom_tile() +
    fill_scale +
    theme_minimal(base_size = base_size) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, color = sample_lab_df$color),
      axis.text.y = element_text(hjust = 0, color = gene_lab_df$color, size = if (show_ylab) base_size else 0),
      axis.ticks.y.right = element_line(),
      axis.text.y.right = element_text(color = gene_lab_df$color),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.margin = margin(2, 0, 2, 0, unit = "pt")
    )

  if (!is.null(vlines)) p_heat <- p_heat + geom_vline(xintercept = vlines, color = "white", linewidth = border_linewidth)
  if (!is.null(hlines)) p_heat <- p_heat + geom_hline(yintercept = hlines, color = "white", linewidth = border_linewidth)

  if (!is.null(row_groups) && !show_ylab) {
    row_annot_df <- data.frame(
      Gene = factor(rownames(mat), levels = rownames(mat)),
      Group = factor(row_groups)
    )
    row_annot_df$y <- as.numeric(row_annot_df$Gene)
    row_annot_df$color <- hue_pal()(length(levels(row_annot_df$Group)))[as.numeric(row_annot_df$Group)]

    p_heat <- p_heat +
      geom_tile(data = row_annot_df, aes(x = 0.5, y = y, fill = Group), width = 0.2, inherit.aes = FALSE) +
      guides(fill = guide_legend(title = "row group")) +
      theme(legend.position = "right")
  }

  if (!is.null(row_clust)) {
    dendro_data <- ggdendro::dendro_data(dend, type = "rectangle")
    row_dendro <- ggplot() +
      geom_segment(data = dendro_data$segments,
                   aes(x = y, xend = yend, y = x, yend = xend),
                   linewidth = 0.2) +
      scale_y_continuous(expand = c(0, 0), breaks = NULL) +
      scale_x_reverse() +
      coord_cartesian(ylim = c(0.5, length(rownames(mat)) + 0.5)) +
      theme_void() +
      theme(plot.margin = margin(2, 0, 2, 0, unit = "pt"))
  }

  if (!is.null(col_clust)) {
    col_dendro_data <- ggdendro::dendro_data(col_dend, type = "rectangle")
    col_dendro <- ggplot() +
      geom_segment(data = col_dendro_data$segments,
                   aes(x = x, xend = xend, y = y, yend = yend),
                   linewidth = 0.2) +
      scale_x_continuous(expand = c(0, 0), breaks = NULL) +
      theme_void() +
      coord_cartesian(xlim = c(0.5, length(colnames(mat)) + 0.5)) +
      theme(plot.margin = margin(2, 0, 2, 0, unit = "pt"))
  }

  heatmap_grob <- ggplotGrob(p_heat)
  if (!is.null(row_clust)) row_dendro_grob <- ggplotGrob(row_dendro)
  if (!is.null(col_clust)) col_dendro_grob <- ggplotGrob(col_dendro)

  if (!is.null(row_clust)) row_dendro_grob$heights <- heatmap_grob$heights
  if (!is.null(col_clust)) col_dendro_grob$widths <- heatmap_grob$widths

  if (!is.null(row_clust) && !is.null(col_clust)) {
    empty <- grid::nullGrob()
    final_plot <- arrangeGrob(
      arrangeGrob(empty, col_dendro_grob, ncol = 2, widths = c(0.05, 1)),
      arrangeGrob(row_dendro_grob, heatmap_grob, ncol = 2, widths = c(0.05, 1)),
      nrow = 2, heights = c(0.05, 1)
    )
    grid::grid.newpage()
    grid::grid.draw(final_plot)
    return(invisible(final_plot))
  } else if (!is.null(row_clust)) {
    final_plot <- arrangeGrob(row_dendro_grob, heatmap_grob, ncol = 2, widths = c(0.05, 1))
    grid::grid.newpage()
    grid::grid.draw(final_plot)
    return(invisible(final_plot))
  } else {
    return(p_heat)
  }
}

set_panel_heights <- function(grob, ref_heights) {
  grob$heights <- ref_heights
  grob
}
