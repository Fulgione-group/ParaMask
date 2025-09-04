#' EM algorithm on geometric mixture model
#'
#' Fits a two-component geometric mixture model using the EM algorithm.
#'
#' @param data A numeric vector of positive integers (distances).
#' @param max_iter Maximum number of EM iterations. Default is 1000.
#' @param tol Convergence tolerance. Default is 1e-6.
#'
#' @return A list with `prob1`, `prob2`, and `weights`, the estimated geometric parameters and posterior weights.
#' @export
em_algorithm_geom <- function(data, max_iter = 1000, tol = 1e-6) {
  n <- length(data)
  k <- 2  # Number of components

  # Initialize parameters
  min_prob <- 1e-5
  prob1 <- runif(1, min = min_prob, max = 1 - min_prob)
  prob2 <- runif(1, min = min_prob, max = 1 - min_prob)
  weights <- runif(n, min = 0.1, max = 0.9)

  for (iter in 1:max_iter) {
    # E-step
    likelihood1 <- dgeom(data, prob = prob1)
    likelihood2 <- dgeom(data, prob = prob2)
    likelihood1[likelihood1 == 0] <- 1e-300
    likelihood2[likelihood2 == 0] <- 1e-300

    resp1 <- weights * likelihood1
    resp2 <- (1 - weights) * likelihood2
    total_resp <- resp1 + resp2
    weights <- resp1 / total_resp

    # M-step
    prob1_new <- min(0.9, sum(weights) / sum(data * weights))
    prob2_new <- min(0.9, sum(1 - weights) / sum(data * (1 - weights)))

    # Check convergence
    if (abs(prob1_new - prob1) < tol && abs(prob2_new - prob2) < tol) {
      break
    }

    # Update parameters
    prob1 <- prob1_new
    prob2 <- prob2_new
  }

  list(prob1 = prob1, prob2 = prob2, weights = weights)
}
