#' Plot EM log-likelihood over iterations
#'
#' Plots the log-likelihood values from EM iterations and saves to a PDF.
#'
#' @param log_likelihood A matrix or data.frame with columns "iteration" and "logLikelihood".
#' @param outpath Output directory.
#' @param ID Optional identifier to include in the filename.
#' @return A saved PDF plot.
#' @import ggplot2
#' @export
plot_em_loglikelihood <- function(log_likelihood, outpath, ID = "") {
  if (!all(c("iteration", "logLikelihood") %in% colnames(log_likelihood))) {
    stop("Input must have 'iteration' and 'logLikelihood' columns.")
  }

  log_likelihood <- as.data.frame(log_likelihood)

  minLL <- min(log_likelihood$logLikelihood)
  maxLL <- max(log_likelihood$logLikelihood)

  # Adjust axis limits
  minLL <- signif(minLL - abs(minLL) / 10, digits = 2)
  maxLL <- signif(maxLL + abs(minLL) / 10, digits = 2)

  # Handle ggplot2 version compatibility
  line_aes_key <- if (utils::packageVersion("ggplot2") >= "3.4.0") "linewidth" else "size"
  legend_position_key <- if (utils::packageVersion("ggplot2") >= "3.5.0") "legend.position.inside" else "legend.position"
  legend_theme <- do.call(ggplot2::theme, setNames(list(c(0.2, 0.6)), legend_position_key))

  ll_plot <- ggplot2::ggplot(data = log_likelihood) +
    geom_point(aes(x = iteration, y = logLikelihood), size = 2) +
    scale_x_continuous(
      limits = c(0, max(log_likelihood$iteration) + 2),
      breaks = seq(1, max(log_likelihood$iteration) + 1, by = 2),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(minLL, maxLL + (maxLL - minLL) / 10),
      breaks = seq(minLL, maxLL, by = (maxLL - minLL) / 5),
      expand = c(0, 0)
    ) +
    do.call(geom_segment, c(list(x = 1, y = minLL, xend = max(log_likelihood$iteration) + 1, yend = minLL),
                            setNames(list(1.2), line_aes_key))) +
    do.call(geom_segment, c(list(x = 0, y = minLL, xend = 0, yend = maxLL),
                            setNames(list(1.2), line_aes_key))) +
    coord_fixed(ratio = 17 / (maxLL - minLL)) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(size = 14),
      axis.ticks = do.call(element_line, setNames(list(1), line_aes_key)),
      axis.ticks.length = unit(0.25, "cm"),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_blank(),
      legend.key = element_rect(color = "white", fill = NA),
      legend.key.size = unit(1.5, "cm"),
    )+
    legend_theme

  # Save to PDF
  outfile <- file.path(outpath, paste0(ID, "_LL.pdf"))
  pdf(file = outfile, width = 16, height = 8)
  print(ll_plot)
  dev.off()
}
