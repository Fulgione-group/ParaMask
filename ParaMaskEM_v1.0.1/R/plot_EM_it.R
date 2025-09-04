#' Plot EM iteration results
#'
#' Generates plots of fitted heterozygosity vs. MAF, with model predictions overlaid.
#'
#' @param het The original het data frame.
#' @param fit The VGAM fit object from the EM loop.
#' @param weight1 A vector of weights for class K.
#' @param rs Indices of SNPs used for plotting (subsample).
#' @param iteration Integer indicating the EM iteration number.
#' @param outpath Path to output folder.
#' @param ID Optional sample or run ID to use in the filename.
#'
#' @return A PDF file with two diagnostic plots.
#' @import ggplot2
#' @export
plot_EM_it <- function(het, fit, weight1, rs, rf, iteration, outpath, ID = "") {

  # Handle ggplot2 version compatibility for line width
  line_aes_key <- if (utils::packageVersion("ggplot2") >= "3.4.0") "linewidth" else "size"
  # Legend position key depending on ggplot2 version
  legend_position_key <- if (utils::packageVersion("ggplot2") >= "3.5.0") "legend.position.inside" else "legend.position"

  # Build the legend theme as a theme object using do.call
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.2, 0.6)), legend_position_key))


  # Prepare prediction data
  predMeans1 <- data.frame(MAF = seq((0.5 / 10000), 0.5, by = (0.5 / 10000)), Z = "K")
  predMeans1$MAF <- 1 - predMeans1$MAF
  predMeans2 <- data.frame(MAF = 0, Z = "D")
  predMeans2 <- predMeans2[rep(1, nrow(predMeans1)), ]

  predMeans1$Mean <- as.vector(VGAM::predictvglm(fit, newdata = predMeans1, type = "response"))
  predMeans2$Mean <- as.vector(VGAM::predictvglm(fit, newdata = predMeans2, type = "response"))

  predMeans1$MAFT <- 1 - predMeans1$MAF
  predMeans2$MAFT <- predMeans1$MAFT
  predMeans1$MeanUT <- predMeans1$Mean
  predMeans2$MeanUT <- predMeans2$Mean
  predMeans1$MeanT <- predMeans1$Mean * predMeans1$MAFT * 2
  predMeans2$MeanT <- predMeans2$Mean * predMeans2$MAFT * 2

  if(!is.null(rf)){
    plot_het <- het[rf,][rs,]
  }else{
    plot_het <- het[rs,]
  }
  # Plot 1
  rplot1 <- ggplot2::ggplot(data = plot_het) +
    geom_point(aes(x = Minor.allele.freq,
                   y = Heterozygous.geno.freq / (2 * Minor.allele.freq),
                   color = weight1[rs])) +
    do.call(geom_path, c(list(data = predMeans1,
                              mapping = aes(x = MAFT, y = MeanUT, linetype = "Paralogs"),
                              color = "black"),
                         setNames(list(1.2), line_aes_key))) +
    do.call(geom_path, c(list(data = predMeans2,
                              mapping = aes(x = MAFT, y = MeanUT, linetype = "SNPs"),
                              color = "black"),
                         setNames(list(1.2), line_aes_key))) +
    scale_x_continuous(limits = c(-0.025, 0.525), breaks = c(0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_linetype_manual(name = "SNP", labels = c("single-copy SNP", "multicopy SNP"),
                          values = c("solid", "twodash")) +
    scale_color_gradient(name = "SNP", limits = c(0, 1), breaks = c(1, 0),
                         labels = c("single-copy SNP", "multicopy SNP"),
                         high = "blue3", low = "brown2") +
    labs(x = "Minor allele frequency (MAF)", y = "Heterozygote frequency / (2×MAF)") +
    coord_fixed(ratio = 0.5) +
    do.call(geom_segment, c(list(aes(x = -0.025, y = -0.003, xend = -0.025, yend = 1.003)),
                            setNames(list(1.2), line_aes_key))) +
    do.call(geom_segment, c(list(aes(x = -0.0013, y = -0.05, xend = 0.5013, yend = -0.05)),
                            setNames(list(1.2), line_aes_key))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.ticks = do.call(element_line, setNames(list(1), line_aes_key)),
      axis.ticks.length = unit(0.25, "cm"),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.key = element_rect(color = "white", fill = NA),
      legend.key.size = unit(1.5, "cm"),
      legend.position = "right"
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)))

  # Plot 2
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.2, 0.6)), legend_position_key))
  rplot2 <- ggplot2::ggplot(data = plot_het) +
    geom_point(aes(x = Minor.allele.freq,
                   y = Heterozygous.geno.freq,
                   color = weight1[rs])) +
    do.call(geom_path, c(list(data = predMeans1,
                              mapping = aes(x = MAFT, y = MeanT, linetype = "Paralogs"),
                              color = "black"),
                         setNames(list(1.2), line_aes_key))) +
    do.call(geom_path, c(list(data = predMeans2,
                              mapping = aes(x = MAFT, y = MeanT, linetype = "SNPs"),
                              color = "black"),
                         setNames(list(1.2), line_aes_key))) +
    scale_x_continuous(limits = c(-0.025, 0.525), breaks = c(0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.05), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_linetype_manual(name = "SNP", labels = c("single-copy SNP", "multicopy SNP"),
                          values = c("solid", "twodash")) +
    scale_color_gradient(name = "SNP", limits = c(0, 1), breaks = c(1, 0),
                         labels = c("single-copy SNP", "multicopy SNP"),
                         high = "blue3", low = "brown2") +
    labs(x = "Minor allele frequency (MAF)", y = "Heterozygote frequency") +
    coord_fixed(ratio = 0.5) +
    do.call(geom_segment, c(list(aes(x = -0.025, y = -0.003, xend = -0.025, yend = 1.003)),
                            setNames(list(1.2), line_aes_key))) +
    do.call(geom_segment, c(list(aes(x = -0.0013, y = -0.05, xend = 0.5013, yend = -0.05)),
                            setNames(list(1.2), line_aes_key))) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.ticks = do.call(element_line, setNames(list(1), line_aes_key)),
      axis.ticks.length = unit(0.25, "cm"),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.key = element_rect(color = "white", fill = NA),
      legend.key.size = unit(1.5, "cm"),
    ) +
    legend_theme +
    guides(color = guide_legend(override.aes = list(size = 4)))

  # Save to PDF
  outfile <- file.path(outpath, paste0(ID, "_EM_iteration", iteration, ".pdf"))
  pdf(file = outfile, width = 7, height = 7)
  print(rplot1)
  print(rplot2)
  dev.off()
}
