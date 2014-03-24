mcmc_file <- paste(zhome,"Functions/MCMCFunctions.r",sep="")
source(mcmc_file)
wd <- paste(zhome,"$REGION$/RUX2_region$SUBREGION$",sep="")
outputfile <- paste(zhome,"output.txt",sep="")
setwd(wd)
incR<-TRUE
incU<-TRUE
incX<-TRUE
Xmode<-2
regionchoice<-"$SUBREGION$"
tps<-$NUMTIMES$
burnin<-$BURNIN$
iters<-$ITERS$
samplefreq<-$SAMPLEFREQ$
tidyup<-FALSE
Initialise("$REGION$")
RLikelihood<-RLikelihoodRUX2
ULikelihood<-ULikelihoodRUX2
XLikelihood<-XLikelihoodRUX2
betaXLikelihood<-betaXLikelihoodRUX2
if (Continue()) {
  cat("Algorithm continues from iteration",st,"at",date(),"\n")
  for (i in st:iters) {
    RUpdate(i)
    UUpdate(i)
    XUpdate(i)
    if (i%%samplefreq==0)
    {
      Sample(i)
    }
    cat("iteration:", i, "\n")
    cat(file=outputfile, "iteration:", i, "\n")
  }
  cat("Algorithm ends",date(),"\n")
  Convergence()
}
