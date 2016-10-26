#' @importFrom stats acf var
NULL

#' Estimated sample size from MCMC chain
#' 
#' Computation taken from the `rstan` package. Included here so that the package
#' can be as self-contained as possible.
#' 
#' @param x A matrix or vector of posterior output for a variable. Rows are samples, columns are chains.
#' @return The expected sample size for the variable (i.e. estimated number of independent samples).
#' @export
ess <- function(x) {
  if (is.vector(x)) {
    dim(x) <- c(length(x), 1)
  }
  chains    <- ncol(x)
  n_samples <- nrow(x)
  acov <- lapply(1:chains, function(i) {
    c <- acf(x[,i], lag.max = n_samples - 1, plot=FALSE, type='covariance')
    c$acf[,,1]
  })
  acov <- do.call(cbind, acov)
  means <- apply(x, 2, mean)
  mean_var <- mean(acov[1,]) * n_samples / (n_samples - 1)
  var_plus <- mean_var * (n_samples - 1) / n_samples;
  if (chains > 1)
    var_plus <- var_plus + var(means)
  rho_hat_sum <- 0
  for (t in 2:nrow(acov)) {
    rho_hat <- 1 - (mean_var - mean(acov[t, ])) / var_plus
    if (is.nan(rho_hat)) rho_hat <- 0
    if (rho_hat < 0) break
    rho_hat_sum <- rho_hat_sum + rho_hat
  }
  ess <- chains * n_samples
  if (rho_hat_sum > 0) ess <- ess / (1 + 2 * rho_hat_sum)
  ess
}

#' Compute the split-R-hat (a measure of convergence) for a variable from posterior samples.
#' 
#' Computation taken from the `rstan` package. Included here so that the package
#' can be as self-contained as possible.
#' 
#' @param x A matrix or vector of posterior output for a variable. Rows are samples, columns are chains.
#' @return The split-R-hat for the variable in question. Should be close to 1 on convergence.
#' @export
rhat <- function(x) {
  if (is.vector(x)) {
    dim(x) <- c(length(x), 1)
  }
  chains <- ncol(x)
  n_samples <- nrow(x)
  half_n <- floor(n_samples / 2)
  idx_2nd <- n_samples - half_n + 1

  split_chain_mean <- numeric(chains * 2)
  split_chain_var <- numeric(chains * 2)

  for (i in 1:chains) {
    split_chain_mean[i] <- mean(x[1:half_n, i])
    split_chain_var[i] <- var(x[1:half_n, i])
    split_chain_mean[chains + i] <- mean(x[idx_2nd:n_samples, i])
    split_chain_var[chains + i] <- var(x[idx_2nd:n_samples, i])
  }
  var_between <- half_n * var(split_chain_mean)
  var_within <- mean(split_chain_var)
  sqrt((var_between/var_within + half_n -1)/half_n)
}
