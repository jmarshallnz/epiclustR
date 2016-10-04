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
control <- list(thinning = 50, chains = 4, samples = 1000 / 4, burnin = 20, parallel = TRUE)

# data
data <- list(cases=cases, popn=n, mbrg=mbrg, nb=weight, rgmb=wch)

# control of proposals
proposal <- list(sigmaR=sigmaR, Rbefore=Rbefore, Rafter=Rafter, Rsigma=Rsigma_eigen,
                 sigmaU=sigmaU,
                 sigmaX=sigmaX)

# fit the model
print(system.time({
posterior <- fit_model(data, state, prior, control=c(control, proposal))
}))

# do analysis
ssapply <- function(x, fun, ...) {
  simplify2array(lapply(x, function(xx, fun, ...) { simplify2array(lapply(xx, fun, ...)) }, fun = fun, ...))
}

extract_variable <- function(x, variable) {
  x[[variable]]
}

# now dump out our posterior chains
R = ssapply(posterior, extract_variable, 'R')
U = ssapply(posterior, extract_variable, 'U')
betaX = ssapply(posterior, extract_variable, 'betaX')
fe = ssapply(posterior, extract_variable, 'fe')
kR = ssapply(posterior, extract_variable, 'kR')
kU = ssapply(posterior, extract_variable, 'kU')
X = ssapply(posterior, extract_variable, 'X')
pX = ssapply(posterior, extract_variable, 'pX')
ecases = ssapply(posterior, cases_per_time, data=data, smoothed=FALSE)
scases = ssapply(posterior, cases_per_time, data=data, smoothed=TRUE)
plot(apply(ecases, 1, median), type='l')
lines(apply(scases, 1, median), col='red')
#hist(apply(R, 1, ess))

#write.table(t(R), 'MidCentral/current_version/R.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(U), 'MidCentral/current_version/U.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(betaX), 'MidCentral/current_version/betaX.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(fe), 'MidCentral/current_version/fixedEffects.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(ecases), 'MidCentral/current_version/expectedCases.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(scases), 'MidCentral/current_version/smoothedCases.txt', row.names=FALSE, col.names=FALSE)

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
