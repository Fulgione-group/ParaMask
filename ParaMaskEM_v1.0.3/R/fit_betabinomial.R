#' Negative log-likelihood for ParaMask beta-binomial model
#'
#' Parameter order:
#'   par[1] = intercept for mu
#'   par[2] = intercept for rho
#'   par[3] = Z == "K" effect
#'   par[4] = MAF slope
paramask_bb_nll <- function(par, regData2, weights) {

  keep <- regData2$N > 1

  dat <- regData2[keep, , drop = FALSE]
  w   <- weights[keep]

  z <- as.numeric(dat$Z == "K")

  ## linear predictor for beta-binomial mean
  eta_mu <- par[1] +
            par[3] * z +
            par[4] * dat$MAF

  ## inverse complementary-log-log link
  mu <- 1 - exp(-exp(eta_mu))

  ## rho is common to both classes
  rho <- stats::plogis(par[2])

  ## avoid numerical 0/1 values
  eps <- 1e-12
  mu  <- pmin(pmax(mu,  eps), 1 - eps)
  rho <- pmin(pmax(rho, eps), 1 - eps)

  log_prob <- VGAM::dbetabinom(
    x    = dat$het,
    size = dat$N,
    prob = mu,
    rho  = rho,
    log  = TRUE
  )

  if (any(!is.finite(log_prob))) {
    return(1e100)
  }

  -sum(w * log_prob)
}



#' Fit ParaMask beta-binomial model
fit_paramask_bb <- function(regData2,
                            weight1,
                            weight2,
                            start = c(2, -2, 0, 1),
                            lower_slope = -Inf,
                            upper_slope = Inf,
                            verbose = FALSE) {

  weights <- c(weight1, weight2)

  lower <- c(
    -Inf,          # mu intercept
    -Inf,          # rho intercept
    -Inf,          # Z effect
    lower_slope    # MAF slope
  )

  upper <- c(
    Inf,
    Inf,
    Inf,
    upper_slope
  )

  opt <- stats::optim(
    par = start,
    fn = paramask_bb_nll,
    regData2 = regData2,
    weights = weights,
    method = "L-BFGS-B",
    lower = lower,
    upper = upper,
    control = list(
      trace = if (verbose) 1 else 0,
      maxit = 1000
    )
  )

  if (opt$convergence != 0) {
    warning(
      "Beta-binomial optimization did not converge: ",
      opt$message
    )
  }

  fit <- list(
    coefficients = opt$par,
    logLik = -opt$value,
    convergence = opt$convergence,
    message = opt$message
  )

  class(fit) <- "paramask_fit"

  fit
}


#' Predict parameters from a ParaMask beta-binomial fit
#'
#' @param fit A paramask_fit object.
#' @param newdata Data frame containing Z and MAF.
#'
#' @return Data frame with beta-binomial mean (mu) and rho.
predict_paramask_bb <- function(fit, newdata) {

  par <- fit$coefficients

  z <- as.numeric(newdata$Z == "K")

  eta_mu <- par[1] +
    par[3] * z +
    par[4] * newdata$MAF

  mu <- 1 - exp(-exp(eta_mu))
  rho <- stats::plogis(par[2])

  eps <- 1e-12
  mu <- pmin(pmax(mu, eps), 1 - eps)
  rho <- pmin(pmax(rho, eps), 1 - eps)

  data.frame(
    mu = mu,
    rho = rep(rho, nrow(newdata))
  )
}
