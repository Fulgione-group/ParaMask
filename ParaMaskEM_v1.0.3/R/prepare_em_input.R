#' Prepare EM input data
#'
#' Prepares regression data and weights for the ParaMaskEM algorithm.
#'
#' @param het A data.frame with heterozygosity statistics.
#' @param subsample Integer: number of variants to subsample for plotting (default = 10000).
#' @param verbose Whether to print debug output.
#' @param nSNPs Number of SNPs to use for EM fitting; 0 uses all eligible SNPs.
#' @param removeZeroHet Whether to exclude SNPs with heterozygote frequency = 0 from EM fitting.
#'
#' @return A list containing regData, regData2, regData3, weights, prediction data, and sampling index.
#' @export
prepare_em_input <- function(het, subsample = 10000, nSNPs = 0, pruneRRD = FALSE, pruneRRD_cutoff = 1.96, removeZeroHet = FALSE, verbose = FALSE) {
  if (verbose) message("Preparing input for EM...")

  # Build regData
  regData <- cbind(
    (1 - het$Minor.allele.freq),
    round(het$Heterozygous.geno.freq * het$Non.missing),
    round(het$Minor.allele.freq * 2 * het$Non.missing),
    sample(c(0, 1), size = nrow(het), replace = TRUE)
  )
  regData <- as.data.frame(regData)
  colnames(regData) <- c("MAF", "het", "N", "cluster")

  # Duplicate for Z (K and D labels)
  regData2 <- rbind(regData, regData)
  regData2$Z <- c(rep("K", nrow(regData)), rep("D", nrow(regData)))
  regData2$MAF[regData2$Z == "D"] <- 0
  colnames(regData2)[1:3] <- c("MAF", "het", "N")

  # Rows eligible for EM fitting
  rf <- seq_len(nrow(het))
  # Optional removal of zero-heterozygosity SNPs from EM fitting
  if (removeZeroHet) {
    rf <- rf[het$Heterozygous.geno.freq[rf] > 0]
    if (verbose) {
      message("Zero-heterozygosity filtering retained ", length(rf), " of ", nrow(het)," SNPs for EM fitting.")
    }
  }
  # Optional RRD pruning
  if (pruneRRD) {
    rf <- rf[is.na(het$Het.allele.deviation[rf]) |abs(het$Het.allele.deviation[rf]) <= pruneRRD_cutoff]

    if (verbose) {
      message(
        "RRD pruning retained ",
        length(rf), " of ", nrow(het),
        " SNPs for EM fitting."
      )
    }
  }

  # Optional subsampling for reduced fitting
  if (nSNPs > 0 && nSNPs < length(rf)) {
    rf <- sample(rf, size = nSNPs, replace = FALSE)
  }

  # If fitting uses only a subset, keep full regData2 for later classification
  if (length(rf) < nrow(het)) {
    regData3 <- regData2

    regData <- regData[rf, ]
    regData2 <- regData2[c(rf, rf + nrow(het)), ]

    # plotting indices relative to reduced fitting data
    rs <- sample(seq_len(length(rf)), size = min(subsample, length(rf)), replace = FALSE)

  } else {
    regData3 <- NULL
    rf <- NULL

    rs <- sample(
      seq_len(nrow(het)),
      size = min(subsample, nrow(het)),
      replace = FALSE
    )
  }
  # Compute binomial cutoffs
  regData$InCutoff <- apply(regData, 1, function(x) {
    qbinom(p = 0.999, size = x[3], prob = x[1], lower.tail = TRUE)
  })
  regData$InCutoff2 <- apply(regData, 1, function(x) {
    qbinom(p = 0.999, size = x[3], prob = x[1], lower.tail = FALSE)
  })

  # Initialize weights
  weight <- runif(n = nrow(regData), min = 0.01, max = 0.99)
  weight[regData$het > regData$InCutoff] <- 0.01
  weight[regData$het == regData$N] <- 0.01
  weight[regData$het < regData$InCutoff2] <- 0.99

  weight1 <- weight
  weight2 <- 1 - weight

  # Prepare prediction data
  predDF <- data.frame(MAF = regData2$MAF, Z = regData2$Z)

  return(list(
    regData = regData,
    regData2 = regData2,
    regData3 = regData3,
    predDF = predDF,
    weight1 = weight1,
    weight2 = weight2,
    rs = rs,
    rf = rf
  ))
}
