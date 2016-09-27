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

# the priors
prior <- list(aR=aR, bR=bR,
              aU=aU, bU=bU,
              aX=aX, bX=bX, abetaX=abetaX, bbetaX=bbetaX)

# burnin and samples are in terms of posterior samples
control <- list(thinning = 50, samples = 1000, burnin = 20)

# data
data <- list(cases=cases, popn=n, mbrg=mbrg, nb=weight, rgmb=wch)

# control of proposals
proposal <- list(sigmaR=sigmaR, Rbefore=Rbefore, Rafter=Rafter, Rsigma=Rsigma_eigen,
                 sigmaU=sigmaU,
                 sigmaX=sigmaX)

# TODO:
# 1. Move the initacceptance code into the package
# 2. Move the MCMC iteration code into the package
# 3. Get parallisation working
# 4. Move model initialisation into the package??
print(system.time({
for (i in seq_len(control$burnin)) {
  state <- update(data, state, prior, c(control, proposal))
  state <- InitAcceptance(state)
  cat("burnin sample:", i, "of", control$burnin, "\n")
}
posterior <- list()
for (i in seq_len(control$samples)) {
  state <- update(data, state, prior, c(control, proposal))
  posterior[[i]] <- state
  state <- InitAcceptance(state)
  cat("posterior sample:", i, "of", control$samples, "\n")
#  cat(file=outputfile, "iteration:", i, "\n")
}
}))

extract_variable <- function(x, variable) {
  x[[variable]]
}
# now dump out our posterior
R = simplify2array(lapply(posterior, extract_variable, 'R'))
U = simplify2array(lapply(posterior, extract_variable, 'U'))
betaX = simplify2array(lapply(posterior, extract_variable, 'betaX'))
fe = simplify2array(lapply(posterior, extract_variable, 'fe'))
kR = simplify2array(lapply(posterior, extract_variable, 'kR'))
kU = simplify2array(lapply(posterior, extract_variable, 'kU'))
X = simplify2array(lapply(posterior, extract_variable, 'X'))
ecases = simplify2array(lapply(posterior, ExpectedCases))
scases = simplify2array(lapply(posterior, ExpectedCases, smoothed=TRUE))

write.table(t(R), 'MidCentral/current_version/R.txt', row.names=FALSE, col.names=FALSE)
write.table(t(U), 'MidCentral/current_version/U.txt', row.names=FALSE, col.names=FALSE)
write.table(t(betaX), 'MidCentral/current_version/betaX.txt', row.names=FALSE, col.names=FALSE)
write.table(t(fe), 'MidCentral/current_version/fixedEffects.txt', row.names=FALSE, col.names=FALSE)
write.table(t(ecases), 'MidCentral/current_version/expectedCases.txt', row.names=FALSE, col.names=FALSE)
write.table(t(scases), 'MidCentral/current_version/smoothedCases.txt', row.names=FALSE, col.names=FALSE)

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
compare_var('expectedCases')
compare_var('smoothedCases')
