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
state <- Initialise("MidCentral")

RLikelihood<-RLikelihoodRUX2
ULikelihood<-ULikelihoodRUX2
XLikelihood<-XLikelihoodRUX2
betaXLikelihood<-betaXLikelihoodRUX2

cat("Algorithm starts",date(),"\n")
for (i in 1:params$iters) {
  state <- RUpdate(i, state)
#  cat("done rupdate, running uupdate\n")
  state <- UUpdate(i, state)
  R  <- state$R
  fe <- state$fe
  U  <- state$U
#  cat("done uupdate, running xupdate\n")
  XUpdate(i)
#  cat("done xupdate\n")
  if (i%%params$samplefreq==0) {
    Sample(i, state)
    state <- InitAcceptance(state)
  }
  cat("iteration:", i, "\n")
#  cat(file=outputfile, "iteration:", i, "\n")
}
cat("Algorithm ends", date(),"\n")
Convergence(state)
