#' Initialize priors
#'
#' @param aR, bR Gamma hyper-prior on precision of R's, defaults to 5, 5/500
#' @param aU, bU Gamma hyper-prior on precision of U's, defaults to 1, 0.5
#' @param aX, bX Beta hyper-prior on probability of an outbreak per year per region, defaults to 1,51
#' @param abetaX, bbetaX Gamma hyper-prior on the size of an outbreak, defaults to 1,1
#' @return a list of priors
#' @export
init_priors <- function(aR = 5, bR = 5/500,
                        aU = 1, bU = 0.5,
                        aX = 1, bX = 51,
                        abetaX = 1, bbetaX = 1) {
  list(aR=aR, bR=bR,
       aU=aU, bU=bU,
       aX=aX, bX=bX, abetaX=abetaX, bbetaX=bbetaX)
}

#' Initialize the MCMC control parameters for fitting the epiclustR model
#'
#' @param samples The number of posterior samples required. These will be distributed
#'                across the chains. Defaults to 1000.
#' @param chains  The number of posterior chains to run. Each chain will generate
#'                $samples/chains$ samples. Defaults to 4.
#' @param thinning The amount of samples to skip per posterior sample for better
#'                 independence between samples, defaults to 50.
#' @param burnin   The number of posterior samples to discard prior to sampling. Defaults to 20
#' @param parallel Whether to use parallel processing. Defaults to TRUE if chains > 1.
#' @param sigmaR   Standard deviation of Metropolis-Hastings random walk jump for R, defaults to 1.
#' @param sigmaU   Standard deviation of Metropolis-Hastings random walk jump for U, defaults to 1.
#' @param sigmaX   Standard deviation of Metropolis-Hastings random walk jump for betaX, defaults to 1.
#' @param blockR   Block sizes for updating R from the conditional prior distribution. Should be
#'                 mutually prime. Defaults to `c(4,5,9,11)`.
#' @return list containing the MCMC control parameters.
#' @export
init_control <- function(samples = 1000, chains = 4, thinning = 50, burnin = 20, parallel = chains > 1,
                         sigmaR = 1, sigmaU = 1, sigmaX = 1, blockR = c(4,5,9,11)) {
  # Generate the block matrices
  blocks <- lapply(blockR, construct_block_matrices)

  list(thinning = thinning,
       chains = chains,
       samples = samples / chains,
       burnin = burnin,
       parallel = parallel,
       sigmaR = sigmaR,
       sigmaU = sigmaU,
       sigmaX = sigmaX,
       blockR = blocks)
}

#' Initialize the model state
#'
#' @param ntimes The number of time points
#' @param nspace The number of spatial units
#' @param nregions The number of regions
#' @param fe The starting value for the fixed effect
#' @param pX The starting value for pX, the probability of an outbreak
#' @param kU,kR The starting values for the precision of spatial and temporal components.
#' @return A list containing the model state
init_state <- function(ntimes, nspace, nregions, fe = -10, pX = 0.1, kU = 1, kR = 1) {
  list(fe = fe,
       R  = rnorm(ntimes,0,1),
       U = rnorm(nspace,0,1),
       X = matrix(rbinom(ntimes*nregions,1,pX), ntimes, nregions),
       betaX = rep(0.2,nregions),
       pX = pX,
       kU = kU,
       kR = kR,
       acceptR = rep(0, 1+4), # TODO: Move these somewhere else...
       rejectR = rep(0, 1+4),
       acceptU = rep(0,2),
       rejectU = rep(0,2),
       acceptX = 0,
       rejectX = 0)
}

#' Reset the acceptance variables in the model state
#'
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

#' Fits an MCMC chain for the epiclustR model for the given dataset
#'
#' @param chain the MCMC chain to run
#' @param data the data for the model
#' @param prior details on the priors for the model
#' @param control details on the control for the model
#' @return the posterior for this chain
fit_chain <- function(chain, data, prior, control) {

  state <- init_state(nrow(data$cases), ncol(data$cases), length(data$rgmb))

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
    posterior[[i]]$ecases <- cases_per_time(data, state, smoothed=FALSE)
    posterior[[i]]$scases <- cases_per_time(data, state, smoothed=TRUE)
    state <- reset_acceptance(state)
    if (chain==1) setTxtProgressBar(pb, i+control$burnin)
  }
  if (chain==1) close(pb)
  posterior
}

#' Checks the data and updates it, e.g. for no population where we have cases
#' 
#' @param data the data for the model
#' @return the checked (and possibly updated) data for the model
#' @export
check_data <- function(data) {
  # initialize the region lookup table
  data$rgmb <- init_region_lut(data$mbrg)

  # check we have population where we have cases
  check_popn <- function(n, cases) {
    wch <- which(n == 0 & apply(cases, 2, sum) > 0)
    if (length(wch > 0)) {
      cat("Setting population to 1 in spatial locations ", wch, " as we have cases there and popn=0\n")
    }
    n[wch] = 1
    n
  }
  data$popn <- check_popn(data$popn, data$cases)

  data
}

#' Fits the epiclustR model for the given dataset
#'
#' @param data the data for the model
#' @param prior details on the priors for the model
#' @param control details on the control for the model
#' @param the random number seed to set. Defaults to NULL (a random seed)
#' @return the posterior for this chain
#' @export
fit_model <- function(data, prior, control, seed = NULL) {
  # make sure the data has been checked
  data <- check_data(data)

  # run our chains
  if (control$parallel) {
    no_cores <- min(parallel::detectCores(), control$chains)
    cluster <- parallel::makeCluster(no_cores, outfile="")
    parallel::clusterSetRNGStream(cluster, seed)
    posterior <- parallel::clusterApply(cluster, 1:control$chains, fit_chain, data=data, prior=prior, control=control)
    parallel::stopCluster(cluster)
  } else {
    set.seed(seed)
    posterior <- lapply(1:control$chains, fit_chain, data=data, prior=prior, control=control)
  }
  posterior
}
