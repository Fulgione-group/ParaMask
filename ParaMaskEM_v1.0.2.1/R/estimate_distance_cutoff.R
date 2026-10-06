#' Estimate distance cutoff from multicopy SNPs
#'
#' Runs EM algorithm on distances between multicopy SNPs to identify cutoff separating clusters.
#'
#' @param het Data frame including Chromosome, Position, and EM_class.
#' @param dist_em_rep Number of EM replicates (e.g. 100).
#' @param min_dist Minimum allowed cutoff.
#' @param max_dist Maximum allowed cutoff.
#' @param cdist Logical. If TRUE, use conservative cutoff.
#' @param em_algorithm_geom Function performing EM on geometric mixture.
#' @param verbose Logical, print progress.
#'
#' @return A list with `dist_cutoff`, `p1`, `p2`, `dist_cutoff_samples`, and `distances`.
#' @export


estimate_distance_cutoff <- function(het, dist_em_rep = 100, min_dist = 5, max_dist = 1000,
                                     cdist = FALSE, verbose = FALSE) {
  if (verbose) message("Calculating distance cutoff...")

  chrlist <- unique(het$Chromosome)
  distances <- c()
  for (chr in chrlist) {
    het_par <- het[het$EM_class == 2 & het$Chromosome == chr, ]
    if (nrow(het_par) > 1) {
      dists <- diff(het_par$Position)
      distances <- c(distances, dists[dists > 0] - 1)
    }
  }

  dist_cutoff_samples <- data.frame()
  for (i in seq_len(dist_em_rep)) {
    result <- em_algorithm_geom(data = distances, max_iter = dist_em_rep)
    p1 <- result[[1]]; p2 <- result[[2]]; w <- result[[3]]
    if (p1 < p2) { tmp <- p1; p1 <- p2; p2 <- tmp }

    cutoff <- NA
    for (j in 0:1000) {
      d1 <- dgeom(j, prob = p1)
      d2 <- dgeom(j, prob = p2)
      if (d2 / d1 > 0.05) {
        conservative_cutoff <- j
        if (d2 / d1 > 1) {
          cutoff <- j
          break
        }
      }
    }

    dist_cutoff <- if (cdist) conservative_cutoff else cutoff
    dist_cutoff_samples <- rbind(dist_cutoff_samples, c(dist_cutoff, p1, p2))
  }

  dist_cutoff_samples <- as.data.frame(dist_cutoff_samples)
  colnames(dist_cutoff_samples) <- c("cutoff", "p1", "p2")

  final_cutoff <- round(median(dist_cutoff_samples$cutoff, na.rm = TRUE))
  final_cutoff <- max(final_cutoff, min_dist)
  final_cutoff <- min(final_cutoff, max_dist)

  closest <- which.min(abs(dist_cutoff_samples$cutoff - final_cutoff))
  p1 <- dist_cutoff_samples$p1[closest]
  p2 <- dist_cutoff_samples$p2[closest]

  return(list(
    dist_cutoff = final_cutoff,
    p1 = p1,
    p2 = p2,
    dist_cutoff_samples = dist_cutoff_samples,
    distances = distances
  ))
}
