#' Plot SNP distance cutoff and geometric model
#'
#' Visualizes observed and fitted distance distributions between multicopy SNPs.
#'
#' @param distances Vector of pairwise distances.
#' @param dist_cutoff Cutoff to highlight on plot.
#' @param p1 Fitted p1 from EM.
#' @param p2 Fitted p2 from EM.
#' @param outpath Output directory.
#' @param ID Optional run ID for filename.
#' @import ggplot2
#' @export
plot_distance_cutoff <- function(distances, dist_cutoff, p1, p2, outpath = ".", ID = "") {
  legend_position_key <- if (utils::packageVersion("ggplot2") >= "3.5.0") "legend.position.inside" else "legend.position"
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.55, 0.8)), legend_position_key))

  md <- max(distances)
  x <- 0:md
  y_emp <- tabulate(factor(distances, levels = x)) / length(distances)
  y1 <- dgeom(x, prob = p1)
  y2 <- dgeom(x, prob = p2)

  disttab <- data.frame(
    bin = rep(x, 3),
    freq = c(y_emp, y1, y2),
    dist = factor(rep(c("empirical", "within haplotype", "between haplotype"), each = length(x)))
  )

  mdf <- max(disttab$freq)
  mdf_y <- round(mdf * 1.1, 2)
  mindf_y <- -round((mdf_y - mdf), 2)

  dist_plot <- ggplot2::ggplot(disttab[disttab$bin < 1001, ], aes(x = bin, y = freq, color = dist)) +
    geom_point() + geom_path() +
    scale_x_continuous(limits = c(-50, 1050), expand = c(0, 0)) +
    scale_y_continuous(limits = c(mindf_y, mdf_y), expand = c(0, 0)) +
    labs(x = "distance", y = "relative frequency") +
    scale_color_manual(values = c("purple", "brown2", "darkorange"),
                       breaks = c("empirical", "within haplotype", "between haplotype")) +
    theme(panel.background = element_rect(fill = NA, color = "white"),
          axis.title = element_text(size = 22),
          axis.text = element_text(size = 18),
          legend.key = element_rect(color = "white", fill = NA),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 16)) +
    legend_theme +
    geom_segment(x = 0, xend = 1000, y = mindf_y, yend = mindf_y, linewidth = 2, color = "black") +
    geom_segment(x = -50, xend = -50, y = 0, yend = mdf_y, linewidth = 2, color = "black") +
    geom_vline(xintercept = dist_cutoff, linetype = "dashed", linewidth = 1) +
    guides(color = guide_legend(override.aes = list(size = 5))) +
    coord_cartesian()

  pdf(file = file.path(outpath, paste0(ID, "_dist.pdf")), width = 16, height = 8)
  print(dist_plot)
  dev.off()
}
