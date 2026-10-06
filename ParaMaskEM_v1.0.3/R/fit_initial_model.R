#' Fit initial ParaMask beta-binomial model
#'
#' @param regData2 Data.frame used for modeling (Z, MAF, het, N).
#' @param weight1 Vector of weights for class K.
#' @param weight2 Vector of weights for class D.
#' @param boundfit Logical: whether to enforce slope constraints.
#' @param lboundary Lower boundary for slope.
#' @param uboundary Upper boundary for slope.
#' @param verbose Whether to print debug messages.
#'
#' @return List containing fitted model, coefficients and offset flag.
#' @export
fit_initial_model <- function(regData2, weight1, weight2,
                              boundfit = FALSE,
                              lboundary = NULL,
                              uboundary = NULL,
                              verbose = FALSE) {

  if (verbose) {
    message("Fitting initial beta-binomial model...")
  }

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

  fit3 <- fit_paramask_bb(
    regData2 = regData2,
    weight1 = weight1,
    weight2 = weight2,
    start = c(2, -2, 0, 1),
    lower_slope = lboundary,
    upper_slope = uboundary,
    verbose = verbose
  )

  coef_fit3 <- fit3$coefficients

  list(
    fit = fit3,
    coef_fit = coef_fit3,
    offsetfit = FALSE
  )
}
