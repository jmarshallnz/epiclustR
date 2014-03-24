# Defaults
lags<-3
kV<-0.1 # prior precision
sigmaV<-1
V<-matrix(0,lags,1)
acceptV<-0
rejectV<-0
if (tidyup) {
  file.remove("V.txt")
  file.remove("acceptanceV.txt")
}
VInitialise3 <- function() {
  VLikelihood<<-VLikelihoodRUV3
  clV<<-list()
  mblV<<-list()
  tplV<<-list()
  for (j in 1:lags) {
    clV[[j]]<<-which(
  }   
}
if (lags==3) {VInitialise3()}
VUpdate <- function(i=0) {
  for (j in 1:lags) {
    # Metropolis Hastings Update
    proposal<-rnorm(V[j],sigmaV)
    ap<-VLikelihood(j,proposal)*exp(-kV*(proposal^2-V[j]^2)/2)
    un<-runif(1)
    if (un<=ap) {
      V[j]<<-proposal
      acceptV<<-acceptV+1
    } else {
      rejectV<<-rejectV+1
    }
  }
}
VLikelihoodRUV3 <- function(j,proposal) {
  if (j==1) {
    prod(dpois(cases[clV[[j]]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+proposal+V[2]*bin12+V[3]*bin13))/dpois(cases[clV[[j]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+V[1]+V[2]*bin12+V[3]*bin13))
  } else if (j==2) {
    prod(dpois(cases[clV[[j]]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+proposal+V[1]*bin21+V[3]*bin23))/dpois(cases[clV[[j]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+V[2]+V[1]*bin21+V[3]*bin23))
  } else if (j==3) {
    prod(dpois(cases[clV[[j]]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+proposal+V[1]*bin31+V[2]*bin32))/dpois(cases[clV[[j]],n[mblV[[j]]]*exp(fe+R[tplV[[j]]]+U[mblV[[j]]]+V[3]+V[1]*bin31+V[2]*bin32))
  }
}
    