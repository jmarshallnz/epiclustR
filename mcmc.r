set.seed(1)
source("Functions/MCMCFunctions.r")

params     <- defaults()
params$incR <- TRUE
params$incU <- TRUE
params$incX <- TRUE
params$Xmode <- 2
params$regionchoice <- "2"

# TODO: create the data files

params$tps   <- 312
params$mbs   <- 1834
params$maxne <- 17
params$baseDeviance <- 21000 # TODO: How is this defined?

params$burnin     <- 10
params$iters      <- 100
params$samplefreq <- 1
params$datapath   <- "."
params$outpath    <- "MidCentral/current_version"
sigmaR<-1

# TODO: Initialize currently loads the data. This needs to be changed
state <- Initialise("MidCentral", setpriors = 2)

RLikelihood<-RLikelihoodRUX2
ULikelihood<-ULikelihoodRUX2
XLikelihood<-XLikelihoodRUX2
betaXLikelihood<-betaXLikelihoodRUX2

XFullOutput <- TRUE

# burnin and samples are in terms of posterior samples
control <- list(thinning = 5, samples = 100, burnin = 2)

print(system.time({
for (i in seq_len(control$burnin)) {
  state <- Update(state, control, burnin=TRUE)
  state <- InitAcceptance(state)
  cat("burnin sample:", i, "of", control$burnin, "\n")
}
for (i in seq_len(control$samples)) {
  state <- Update(state, control)
  Sample(state)
  state <- InitAcceptance(state)
  cat("posterior sample:", i, "of", control$samples, "\n")
#  cat(file=outputfile, "iteration:", i, "\n")
}
}))
#Convergence(state)

compare_var <- function(variable) {
  ref <- scan(file.path('MidCentral/RUX2_region2', paste0(variable, '.txt')))
  new <- scan(file.path('MidCentral/current_version', paste0(variable, '.txt')))
  print(all.equal(ref, new))
}

# compare versions
compare_var('fixedEffects')
compare_var('R')
compare_var('U')
compare_var('betaX')
compare_var('cumulativeX')
compare_var('expectedCases')
compare_var('smoothedCases')
