#' reset_acceptance
#' @param state the current state of the Markov chain
#' @return the state, with the acceptance probabilities reset
#' @export
reset_acceptance <- function(state) {
  state$acceptR <- rep(0, length(state$acceptR))
  state$acceptU <- rep(0, length(state$acceptU))
  state$acceptX <- rep(0, length(state$acceptX))
  state$rejectR <- rep(0, length(state$rejectR))
  state$rejectU <- rep(0, length(state$rejectU))
  state$rejectX <- rep(0, length(state$rejectX))
  return(state)
}

#' fit_chain
#' 
#' Fits an MCMC chain for the epiclustR model for the given dataset
#' @param chain the MCMC chain to run
#' @param data the data for the model
#' @param state the current state of the Markov chain
#' @param prior details on the priors for the model
#' @param control details on the control for the model
#' @return the posterior for this chain
fit_chain <- function(chain, data, state, prior, control) {
  if (chain == 1)
    pb <- txtProgressBar(min=0, max=control$burnin+control$samples, width=56, style=3)
  for (i in seq_len(control$burnin)) {
    state <- update(data, state, prior, control)
    state <- reset_acceptance(state)
    if (chain==1) setTxtProgressBar(pb, i)
  }
  posterior <- list()
  for (i in seq_len(control$samples)) {
    state <- update(data, state, prior, control)
    posterior[[i]] <- state
    state <- reset_acceptance(state)
    if (chain==1) setTxtProgressBar(pb, i+control$burnin)
  }
  if (chain==1) close(pb)
  posterior
}

#' fit_model
#' 
#' Fits the epiclustR model for the given dataset
#' @param data the data for the model
#' @param state the current state of the Markov chain
#' @param prior details on the priors for the model
#' @param control details on the control for the model
#' @return the posterior for this chain
#' @export
fit_model <- function(data, state, prior, control) {
  if (control$parallel) {
    no_cores <- min(parallel::detectCores(), control$chains)
    cluster <- parallel::makeCluster(no_cores, outfile="")
    parallel::clusterSetRNGStream(cluster, 1)
    posterior <- parallel::clusterApply(cluster, 1:control$chains, fit_chain, data=data, state=state, prior=prior, control=control)
    parallel::stopCluster(cluster)
  } else {
    set.seed(1)
    posterior <- lapply(1:control$chains, fit_chain, data=data, state=state, prior=prior, control=control)
  }
  posterior
}
