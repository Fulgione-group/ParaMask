#' Run EM loop for ParaMaskEM
#'
#' Iteratively refits the VGAM model using updated class weights until convergence.
#'
#' @param regData2 Data frame used for modeling (Z, MAF, het, N).
#' @param predDF Data frame used for prediction (contains MAF and Z).
#' @param fit Initial VGAM model fit.
#' @param coef_fit Initial coefficients from model fit.
#' @param offsetfit Logical: was the initial model fit using a fixed offset?
#' @param weight1 Vector of initial weights for class K.
#' @param weight2 Vector of initial weights for class D.
#' @param maxiter Maximum number of EM iterations.
#' @param boundfit Logical: whether to enforce slope boundaries.
#' @param lboundary Lower slope boundary.
#' @param uboundary Upper slope boundary.
#' @param tolerance Convergence threshold.
#' @param verbose Print progress?
#'
#' @return A list with the final fit, coefficients, log-likelihood, and offsetfit status.
#' @import VGAM
#' @export
run_em_loop <- function(regData2, predDF, fit, coef_fit, offsetfit,
                        weight1, weight2,
                        maxiter = 100,
                        boundfit = FALSE, lboundary = NULL, uboundary = NULL,
                        tolerance = 0.001,
                        verbose = FALSE) {

  if (verbose) message("Starting EM loop...")

  log_likelihood <- c()

  for (iteration in 1:maxiter) {
    if (verbose) message("EM Iteration: ", iteration)

    # --- E-step ---
    pars3 <- data.frame(mu = rep(NA, nrow(predDF)), rho = rep(NA, nrow(predDF)))

    if (boundfit && offsetfit) {
      pars3$mu <- apply(predDF, 1, function(x) {
        z <- ifelse(x["Z"] == "K", 1, 0)
        lp <- coef_fit[1] + z * coef_fit[3] + as.numeric(x["MAF"]) * z * coef_fit[4]
        VGAM::cloglog(lp, inverse = TRUE)
      })
      pars3$rho <- VGAM::logitlink(coef_fit[2], inverse = TRUE)
      offsetfit <- FALSE
    } else {
      pars3 <- VGAM::predictvglm(fit, newdata = predDF, type = "link", untransform = TRUE)
      pars3 <- as.data.frame(pars3)
      names(pars3) <- c("mu", "rho")
    }

    regData2$pred_fit_mu <- pars3$mu
    regData2$pred_fit_rho <- pars3$rho

    # Calculate betabinomial probabilities
    probs <- apply(regData2, 1, function(x) {
      VGAM::dbetabinom(x = as.numeric(x["het"]),
                       size = as.numeric(x["N"]),
                       prob = as.numeric(x["pred_fit_mu"]),
                       rho = as.numeric(x["pred_fit_rho"]))
    })

    probs1 <- probs[regData2$Z == "K"]
    probs2 <- probs[regData2$Z == "D"]

    weight1 <- (weight1 * probs1) / (weight1 * probs1 + weight2 * probs2)
    weight2 <- (weight2 * probs2) / (weight1 * probs1 + weight2 * probs2)

    weight1[weight1 < 1e-100] <- 1e-100
    weight2[weight2 < 1e-100] <- 1e-100

    # --- M-step ---
    if (verbose) message("M-step: refitting model...")

    fit <- tryCatch({
      VGAM::vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + MAF,
                 data = regData2,
                 family = VGAM::betabinomial(lmu = "clogloglink",
                                             lrho = "logitlink"),
                 trace = verbose,
                 subset = regData2$N > 1,
                 weights = c(weight1, weight2))
    }, error = function(e) {
      stop("ERROR in M-step during VGAM fitting: ", conditionMessage(e))
    })

    coef_fit_new <- stats::coefficients(fit)

    # Optional boundary refit
    if (boundfit && (coef_fit_new[4] > uboundary || coef_fit_new[4] < lboundary)) {
      cboundary <- ifelse(coef_fit_new[4] > uboundary, uboundary, lboundary)
      if (verbose) message("Refitting with slope boundary offset: ", cboundary)

      fit <- tryCatch({
        VGAM::vglm(cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF * cboundary),
                   data = regData2,
                   family = VGAM::betabinomial("clogloglink", "logitlink"),
                   start = c(2, 1, 0, 2),
                   trace = verbose,
                   subset = regData2$N > 1,
                   weights = c(weight1, weight2))
      }, error = function(e) {
        stop("ERROR in offset M-step during VGAM fitting: ", conditionMessage(e))
      })

      coef_fit_new <- stats::coefficients(fit)[1:3]
      coef_fit_new[4] <- cboundary
      offsetfit <- TRUE
    }

    # Convergence check
    if (verbose) {
      message("Old parameters: ", paste(round(coef_fit, 4), collapse = ", "))
      message("New parameters: ", paste(round(coef_fit_new, 4), collapse = ", "))
    }

    if (all(abs(coef_fit_new - coef_fit) < tolerance)) {
      if (verbose) message("EM converged at iteration ", iteration)
      break
    }

    # Update coefficients
    coef_fit <- coef_fit_new

    # Store likelihood
    log_likelihood <- rbind(log_likelihood, c(iteration, stats::logLik(fit)))
    colnames(log_likelihood) <- c("iteration", "logLikelihood")
  }

  return(list(
    fit = fit,
    coef_fit = coef_fit,
    offsetfit = offsetfit,
    log_likelihood = log_likelihood,
    iteration = iteration,
    weight1 = weight1,
    probs1 = probs1,
    probs2 = probs2
  ))
}
