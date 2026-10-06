#' Run EM loop for ParaMaskEM
#'
#' Iteratively refits the beta-binomial model using updated class weights
#' until convergence.
#'
#' @param regData2 Data frame used for modeling (Z, MAF, het, N).
#' @param predDF Data frame used for prediction (contains MAF and Z).
#' @param fit Initial ParaMask beta-binomial fit.
#' @param coef_fit Initial coefficients from model fit.
#' @param offsetfit Retained for backwards compatibility. No longer used.
#' @param weight1 Vector of initial weights for class K.
#' @param weight2 Vector of initial weights for class D.
#' @param maxiter Maximum number of EM iterations.
#' @param boundfit Logical: whether to enforce slope boundaries.
#' @param lboundary Lower slope boundary.
#' @param uboundary Upper slope boundary.
#' @param tolerance Convergence threshold.
#' @param verbose Print progress?
#'
#' @return A list with the final fit, coefficients, log-likelihood,
#'         class weights and probabilities.
#'
#' @export
run_em_loop <- function(regData2, predDF, fit, coef_fit, offsetfit,
                        weight1, weight2,
                        maxiter = 100,
                        boundfit = FALSE,
                        lboundary = NULL,
                        uboundary = NULL,
                        tolerance = 0.001,
                        verbose = FALSE) {

  if (verbose) {
    message("Starting EM loop...")
  }

  ## offset fitting is no longer necessary because the optimizer
  ## directly enforces coefficient boundaries
  offsetfit <- FALSE

  ## Set slope bounds
  if (boundfit) {

    if (is.null(lboundary)) {
      lboundary <- -Inf
    }

    if (is.null(uboundary)) {
      uboundary <- Inf
    }

  } else {

    lboundary <- -Inf
    uboundary <- Inf
  }

  log_likelihood <- matrix(
    numeric(0),
    ncol = 2,
    dimnames = list(NULL, c("iteration", "logLikelihood"))
  )

  probs1 <- NULL
  probs2 <- NULL

  for (iteration in seq_len(maxiter)) {

    if (verbose) {
      message("EM iteration: ", iteration)
    }

    # ------------------------------------------------------------
    # E-step
    # ------------------------------------------------------------

    pars3 <- predict_paramask_bb(
      fit = fit,
      newdata = predDF
    )

    regData2$pred_fit_mu <- pars3$mu
    regData2$pred_fit_rho <- pars3$rho

    ## Calculate beta-binomial probabilities
    probs <- VGAM::dbetabinom(
      x = regData2$het,
      size = regData2$N,
      prob = regData2$pred_fit_mu,
      rho = regData2$pred_fit_rho
    )

    probs1 <- probs[regData2$Z == "K"]
    probs2 <- probs[regData2$Z == "D"]

    ## ------------------------------------------------------------
    ## Update posterior class probabilities
    ##
    ## IMPORTANT:
    ## calculate denominator BEFORE changing either weight vector.
    ## ------------------------------------------------------------

    denom <- weight1 * probs1 +
      weight2 * probs2

    ## Protect against numerical underflow
    denom[!is.finite(denom) | denom <= 0] <- .Machine$double.xmin

    weight1_new <- (weight1 * probs1) / denom
    weight2_new <- (weight2 * probs2) / denom

    ## Numerical protection
    eps_weight <- 1e-100

    weight1_new <- pmax(weight1_new, eps_weight)
    weight2_new <- pmax(weight2_new, eps_weight)

    ## Renormalize after applying floor
    total_weight <- weight1_new + weight2_new

    weight1 <- weight1_new / total_weight
    weight2 <- weight2_new / total_weight

    # ------------------------------------------------------------
    # M-step
    # ------------------------------------------------------------

    if (verbose) {
      message("M-step: refitting beta-binomial model...")
    }

    fit_new <- tryCatch({

      fit_paramask_bb(
        regData2 = regData2,
        weight1 = weight1,
        weight2 = weight2,

        ## use previous EM estimates as starting values
        start = coef_fit,

        lower_slope = lboundary,
        upper_slope = uboundary,

        verbose = verbose
      )

    }, error = function(e) {

      stop(
        "ERROR in beta-binomial M-step: ",
        conditionMessage(e)
      )

    })

    coef_fit_new <- fit_new$coefficients

    # ------------------------------------------------------------
    # Store likelihood
    # ------------------------------------------------------------

    log_likelihood <- rbind(
      log_likelihood,
      c(
        iteration = iteration,
        logLikelihood = fit_new$logLik
      )
    )

    # ------------------------------------------------------------
    # Diagnostics
    # ------------------------------------------------------------

    if (verbose) {

      message(
        "Old parameters: ",
        paste(round(coef_fit, 4), collapse = ", ")
      )

      message(
        "New parameters: ",
        paste(round(coef_fit_new, 4), collapse = ", ")
      )

      if (boundfit) {

        if (abs(coef_fit_new[4] - lboundary) < 1e-8) {
          message("MAF slope is at lower boundary: ", lboundary)
        }

        if (abs(coef_fit_new[4] - uboundary) < 1e-8) {
          message("MAF slope is at upper boundary: ", uboundary)
        }
      }
    }

    # ------------------------------------------------------------
    # Convergence
    # ------------------------------------------------------------

    parameter_change <- abs(coef_fit_new - coef_fit)

    ## IMPORTANT:
    ## update fit/coefs BEFORE potentially breaking
    fit <- fit_new
    coef_fit <- coef_fit_new

    if (all(parameter_change < tolerance)) {

      if (verbose) {
        message("EM converged at iteration ", iteration)
      }

      break
    }
  }

  return(list(
    fit = fit,
    coef_fit = coef_fit,
    offsetfit = FALSE,
    log_likelihood = as.data.frame(log_likelihood),
    iteration = iteration,
    weight1 = weight1,
    weight2 = weight2,
    probs1 = probs1,
    probs2 = probs2
  ))
}
