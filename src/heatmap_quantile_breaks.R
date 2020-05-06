#From https://slowkow.com/notes/pheatmap-tutorial/#quantile-breaks
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

#mat_breaks <- quantile_breaks(mat, n = 11)

#pheatmap(
#  mat               = mat,
#  color             = inferno(length(mat_breaks) - 1),
# OR
#  color              = colorRampPalette(rev(brewer.pal(n = 11,
                                        #name = "RdYlBu")))(length(mat_breaks) - 1),
#  breaks            = mat_breaks,
#  border_color      = NA,
#  show_colnames     = FALSE,
#  show_rownames     = FALSE,
#  annotation_col    = mat_col,
#  annotation_colors = mat_colors,
#  drop_levels       = TRUE,
#  fontsize          = 14,
#  main              = "Quantile Color Scale"
#)