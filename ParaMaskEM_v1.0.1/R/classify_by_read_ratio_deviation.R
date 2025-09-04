#' Classify SNPs by read ratio deviation
#'
#' Uses a normal distribution CI based on single-copy SNPs to reclassify ambiguous SNPs.
#'
#' @param het A data.frame including columns `Het.allele.deviation` and `EM_class`.
#' @param verbose Logical, whether to print info.
#'
#' @return A modified `het` with updated `EM_class` and `allele.deviation.seed`.
#' @export
classify_by_read_ratio_deviation <- function(het, verbose = FALSE) {
  ad_mu <- mean(het$Het.allele.deviation[het$EM_class == 0], na.rm = TRUE)
  ad_sd <- max(sd(het$Het.allele.deviation[het$EM_class == 0], na.rm = TRUE), 1)

  q1 <- qnorm(0.025, mean = ad_mu, sd = ad_sd)  # lower
  q2 <- qnorm(0.975, mean = ad_mu, sd = ad_sd)  # upper

  if (verbose) {
    message("Allele ratio CI mean = ", ad_mu)
    message("Allele ratio CI SD   = ", ad_sd)
    message("Lower CI cutoff      = ", q1)
    message("Upper CI cutoff      = ", q2)
  }

  het$allele.deviation.seed <- 0
  reclassify <- het$EM_class == 1 & (het$Het.allele.deviation < q1 | het$Het.allele.deviation > q2)
  het$allele.deviation.seed[reclassify] <- 1
  het$EM_class[reclassify] <- 2

  het$EM_class <- factor(het$EM_class, levels = c(0, 1, 2))

  return(list(
         het = het,
         q1 = q1,
         q2 = q2
  ))
}
