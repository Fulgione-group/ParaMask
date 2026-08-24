#' Plot read ratio deviation and allele ratio density
#'
#' Visualizes read ratio deviations and allele ratio distributions across EM classes.
#'
#' @param het A data.frame including EM_class, Het.allele.ratio, Het.allele.deviation, Heterozygous.geno.freq.
#' @param rs A vector of row indices (e.g. subsample) to plot.
#' @param outpath Directory to save output.
#' @param ID Label for filenames.
#' @param q1 Optional lower cutoff line (e.g. from CI).
#' @param q2 Optional upper cutoff line.
#' @import ggplot2
#' @import patchwork
#' @export
plot_read_ratio_deviation <- function(het, rs, outpath, ID = "", q1 = NULL, q2 = NULL) {
  legend_position_key <- if (utils::packageVersion("ggplot2") >= "3.5.0") "legend.position.inside" else "legend.position"
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.5, 0.9)), legend_position_key))

  plot_data <- het[rs, ]
  plot_data <- plot_data[!is.na(plot_data$Het.allele.ratio) & !is.na(plot_data$Het.allele.deviation), ]

  # Allele ratio plot
  dplot <- ggplot2::ggplot(plot_data) +
    geom_point(aes(x = Heterozygous.geno.freq, y = Het.allele.ratio,
                   color = EM_class, shape = EM_class, fill = EM_class),
               stroke = 1, alpha = 0.2, size = 0.7) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    labs(x = "Heterozygote frequency", y = "Read ratio deviation (D)") +
    scale_color_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_fill_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_shape_manual(values = c(0, 1, 2)) +
    theme(panel.background = element_rect(fill = NA, color = "white"),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = NA, color = "white"),
          legend.key = element_rect(fill = NA, color = "white")) +
    legend_theme +
    geom_segment(x = -0.05, xend = -0.05, y = 0, yend = 1, color = "black") +
    geom_segment(x = 0, xend = 1, y = -0.05, yend = -0.05, color = "black") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

  dens1 <- ggplot2::ggplot(plot_data, aes(x = Het.allele.ratio, color = EM_class, fill = EM_class)) +
    geom_density(alpha = 0.4) +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    scale_color_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_fill_manual(values = c("blue3", "darkgreen", "brown2")) +
    theme_void() +
    theme(legend.position = "none") +
    coord_flip()

  ARplot <- (dplot | dens1) + patchwork::plot_layout(ncol = 2, widths = c(4, 1))

  pdf(file = file.path(outpath, paste0(ID, "_AR.pdf")), width = 16, height = 8)
  print(ARplot)
  dev.off()

  # RRD deviation plot
  dplot2 <- ggplot2::ggplot(plot_data) +
    geom_point(aes(x = Heterozygous.geno.freq, y = Het.allele.deviation,
                   color = EM_class, shape = EM_class, fill = EM_class),
               stroke = 1, alpha = 0.2, size = 0.7) +
    geom_hline(yintercept = q1, linetype = "dashed") +
    geom_hline(yintercept = q2, linetype = "dashed") +
    scale_x_continuous(limits = c(-0.05, 1.05), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-85, 85), breaks = seq(-80, 80, by = 20), expand = c(0, 0)) +
    labs(x = "Heterozygote frequency", y = "Read ratio deviation (D)") +
    scale_color_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_fill_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_shape_manual(values = c(0, 1, 2)) +
    theme(panel.background = element_rect(fill = NA, color = "white"),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          legend.title = element_blank(),
          legend.direction = "horizontal",
          legend.background = element_rect(fill = NA, color = "white"),
          legend.key = element_rect(fill = NA, color = "white")) +
    legend_theme +
    geom_segment(x = -0.05, xend = -0.05, y = -80, yend = 80, color = "black") +
    geom_segment(x = 0, xend = 1, y = -85, yend = -85, color = "black") +
    guides(color = guide_legend(override.aes = list(size = 5, alpha = 1)))

  dens2 <- ggplot2::ggplot(plot_data, aes(x = Het.allele.deviation, color = EM_class, fill = EM_class)) +
    geom_density(alpha = 0.4) +
    scale_x_continuous(limits = c(-85, 85), breaks = seq(-80, 80, by = 20), expand = c(0, 0)) +
    scale_color_manual(values = c("blue3", "darkgreen", "brown2")) +
    scale_fill_manual(values = c("blue3", "darkgreen", "brown2")) +
    theme_void() +
    theme(legend.position = "none") +
    coord_flip()

  RRDplot <- (dplot2 | dens2) + patchwork::plot_layout(ncol = 2, widths = c(4, 1))

  pdf(file = file.path(outpath, paste0(ID, "_RRD.pdf")), width = 16, height = 8)
  print(RRDplot)
  dev.off()
}
