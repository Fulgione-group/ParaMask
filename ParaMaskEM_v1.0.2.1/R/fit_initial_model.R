#' Fit initial VGAM model
#'
#' Fits a betabinomial model using VGAM on the initial EM input data.
#'
#' @param regData2 Data.frame used for modeling (Z, MAF, het, N).
#' @param weight1 Vector of weights for class K.
#' @param weight2 Vector of weights for class D.
#' @param boundfit Logical: whether to enforce boundary constraints.
#' @param lboundary Lower boundary for slope (if `boundfit = TRUE`).
#' @param uboundary Upper boundary for slope (if `boundfit = TRUE`).
#' @param verbose Whether to print debug messages.
#'
#' @return A list with the fitted model (`fit3`), coefficient vector (`coef_fit3`), and flag (`offsetfit`).
#' @import VGAM
#' @export
fit_initial_model <- function(regData2, weight1, weight2,
                              boundfit = FALSE, lboundary = NULL, uboundary = NULL,
                              verbose = FALSE) {
  if (verbose) message("Fitting initial model...")

  offsetfit <- FALSE
  fit3 <- NULL

  # Fit model
  tryCatch({
    fit3 <- VGAM::vglm(
      cbind(regData2$het, regData2$N - regData2$het) ~ Z + MAF,
      data = regData2,
      family = VGAM::betabinomial(lmu = "clogloglink", lrho = "logitlink"),
      trace = verbose,
      subset = (regData2$N > 1),
      weights = c(weight1, weight2)
    )
  }, error = function(e) {
    stop("Error in initial model fitting: ", conditionMessage(e))
  })

  coef_fit3 <- stats::coefficients(fit3)

  # Check boundaries and possibly refit with offset
  if (boundfit && !is.null(coef_fit3[4]) &&
      (coef_fit3[4] > uboundary || coef_fit3[4] < lboundary)) {

    if (verbose) {
      message("Initial slope outside boundaries, refitting with offset...")
    }

    cboundary <- ifelse(coef_fit3[4] > uboundary, uboundary, lboundary)

    tryCatch({
      fit3 <- VGAM::vglm(
        cbind(regData2$het, regData2$N - regData2$het) ~ Z + offset(MAF * cboundary),
        data = regData2,
        family = VGAM::betabinomial(lmu = "clogloglink", lrho = "logitlink"),
        start = c(2, 1, 0, 2),
        trace = verbose,
        subset = (regData2$N > 1),
        weights = c(weight1, weight2)
      )

      coef_fit3 <- stats::coefficients(fit3)[1:3]
      coef_fit3[4] <- cboundary
      offsetfit <- TRUE

    }, error = function(e) {
      stop("Error during boundary refit: ", conditionMessage(e))
    })
  }

  return(list(
    fit = fit3,
    coef_fit = coef_fit3,
    offsetfit = offsetfit
  ))
}
