#' Classify SNPs based on log-likelihood ratio
#'
#' Classifies each SNP based on the LLR between the two EM models.
#'
#' @param het A data.frame containing `probs1` (K) and `probs2` (D) columns.
#' @param probs1 Vector of likelihoods under the single-copy model.
#' @param probs2 Vector of likelihoods under the multicopy model.
#' @param cutoff Numeric LLR threshold corresponding to 0.99/0.01 (default).
#' @param regData3 If a subset of nSNPs was used to fit model parameters, then other SNPs in regData3 are also classified.
#' @param fit fitted model parameters.
#'
#' @return A modified `het` data.frame with `L1`, `L2`, `LLR`, and `EM_class` columns.
#' @export
classify_snps_by_llr <- function(het,probs1, probs2, regData3 = NULL, fit = NULL, cutoff = 0.99) {

  # if nSNPs is active regData3 and fit will be set
  if( !is.null(regData3)  & !is.null(fit) ){
    predDF <- data.frame(MAF = regData3$MAF, Z = regData3$Z)

    pars3 <- data.frame(mu = rep(NA, nrow(predDF)), rho = rep(NA, nrow(predDF)))

    pars3 <- VGAM::predictvglm(fit, newdata = predDF, type = "link", untransform = TRUE)
    pars3 <- as.data.frame(pars3)
    names(pars3) <- c("mu", "rho")

    regData3$pred_fit_mu <- pars3$mu
    regData3$pred_fit_rho <- pars3$rho
    probs <- apply(regData3, 1, function(x) {
      VGAM::dbetabinom(x = as.numeric(x["het"]),
                       size = as.numeric(x["N"]),
                       prob = as.numeric(x["pred_fit_mu"]),
                       rho = as.numeric(x["pred_fit_rho"]))
    })

    probs1 <- probs[regData3$Z == "K"]
    probs2 <- probs[regData3$Z == "D"]

  }

  # Avoid log(0)
  probs1[probs1 == 0] <- 1e-40
  probs2[probs2 == 0] <- 1e-40

  # Assign log-likelihoods and ratio
  het$L1 <- probs1
  het$L2 <- probs2
  het$LLR <- log(het$L1 / het$L2)

  # Calculate symmetric cutoffs
  LLR_CO1 <- log(cutoff / (1 - cutoff))
  LLR_CO2 <- log((1 - cutoff) / cutoff)

  # Classify:
  # 0 = likely single-copy
  # 1 = ambiguous
  # 2 = likely multicopy
  het$EM_class <- 0
  het$EM_class[het$LLR < LLR_CO1 & het$LLR > LLR_CO2] <- 1
  het$EM_class[het$LLR <= LLR_CO2] <- 2

  het$EM_class <- factor(het$EM_class, levels = c(0, 1, 2))

  return(het)
}
