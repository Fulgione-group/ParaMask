#' Plot SNP LLR classification
#'
#' Plots the log-likelihood ratio (LLR) for each SNP across MAF, colored by EM class.
#'
#' @param het A data.frame that includes columns `LLR`, `Minor.allele.freq`, and `EM_class`.
#' @param rs Vector of row indices (e.g., subsampled SNPs) to plot.
#' @param LLR_CO1 Numeric. Upper LLR cutoff (default: log(0.99/0.01)).
#' @param LLR_CO2 Numeric. Lower LLR cutoff (default: log(0.01/0.99)).
#' @param outpath Output directory.
#' @param ID Optional string to label output file.
#' @return A PDF file with the LLR classification plot.
#' @import ggplot2
#' @export

plot_llr_classification <- function(het, rs, cutoff = 0.99, outpath = ".", ID = "") {

  legend_position_key <- if (utils::packageVersion("ggplot2") >= "3.5.0") "legend.position.inside" else "legend.position"
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.15, 0.8)), legend_position_key))

  LLR_CO1 <- log(cutoff / (1 - cutoff))
  LLR_CO2 <- log((1 - cutoff) / cutoff)
  # Determine plotting bounds
  minLLR <- min(het$LLR) - abs(diff(range(het$LLR))) / 20
  maxLLR <- max(het$LLR) + abs(diff(range(het$LLR))) / 20
  minLLR <- signif(minLLR, 2)
  maxLLR <- signif(maxLLR, 2)

  breakstep <- floor((maxLLR - minLLR) / 10)
  minLLR <- minLLR - minLLR %% breakstep
  maxLLR <- maxLLR + breakstep - maxLLR %% breakstep

  # Build plot
  LLR_plot <- ggplot2::ggplot(data = het[rs, ], aes(x = Minor.allele.freq, y = LLR, color = EM_class)) +
    geom_point(size = 4) +
    scale_x_continuous(expand = c(0, 0), limits = c(-0.015, 0.515),
                       breaks = seq(0, 0.5, by = 0.25), labels = c(0, 0.25, 0.5)) +
    scale_y_continuous(limits = c(minLLR - (maxLLR - minLLR) / 20,
                                  maxLLR + (maxLLR - minLLR) / 20),
                       expand = c(0, 0),
                       breaks = seq(minLLR, maxLLR, by = breakstep)) +
    theme(panel.background = element_rect(fill = NA, color = "white"),
          axis.title = element_text(size = 22),
          axis.text = element_text(size = 18),
          axis.text.x = element_text(vjust = 0.5),
          legend.key = element_rect(color = "white", fill = NA),
          legend.title = element_blank(),
          legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 16)) +
    legend_theme +
    geom_segment(x = 0, y = minLLR - (maxLLR - minLLR) / 20,
                 xend = 0.5, yend = minLLR - (maxLLR - minLLR) / 20,
                 linewidth = 2.4, color = "black") +
    geom_segment(x = -0.015, y = minLLR, xend = -0.015, yend = maxLLR,
                 linewidth = 2.4, color = "black") +
    scale_color_manual(
      labels = c("high confidence single copy SNP", "Uncertain", "high confidence paralogous SNP"),
      values = c("blue3", "darkgreen", "brown2")
    ) +
    labs(y = "logLikelihood-ratio", x = "minor allele frequency (maf)") +
    geom_hline(yintercept = LLR_CO1, linewidth = 1.4) +
    annotate("text", x = 0.45, y = LLR_CO2 - 4, label = "LLR lower cutoff", size = 6) +
    geom_hline(yintercept = LLR_CO2, linewidth = 1.4) +
    annotate("text", x = 0.45, y = LLR_CO1 + 4, label = "LLR upper cutoff", size = 6)

  # Output file
  outfile <- file.path(outpath, paste0(ID, "_LLR.pdf"))
  pdf(file = outfile, width = 16, height = 8)
  print(LLR_plot)
  dev.off()
}
