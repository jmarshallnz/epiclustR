# Defaults
sigmaU<-1
kU<-1
U<-rnorm(params$mbs,0,1)

UInitialize <- function() {
  acceptU<<-matrix(0,2,1)
  rejectU<<-matrix(0,2,1)
  if (params$tidyup) {file.remove(file.path(params$outpath, "U.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "kU.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "acceptanceU.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "sumU.txt"))}
}

# computes the variance of the U's for the Gibbs gamma update
USumFunction <- function(U) {
  s<-0
  for (i in seq_along(U)) {
    s <- s + sum((U[i]-U[weight[i,2:(1+weight[i,1])]])^2)
  }
  s/2 # TODO: Why the divide by two here? Shouldn't that be done in the rate in the Gamma fnc?
}

USetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aU<<-1
    bU<<-10^-2
  } else if (setpriors==2) {
    aU<<-1
    bU<<-0.5
  } else if (setpriors==3) {
    aU<<-10^-4
    bU<<-10^-4
  }     
}

UUpdate <- function(i=0) {
  lenU <- length(U)
  # Gibb's step to update kU
  kU<<-rgamma(1,aU+(lenU-1)/2,rate=(bU+USumFunction(U)/2))
  for (j in 1:lenU) {
    if (i%%2==0) {
      proposal<-rnorm(1,mean(U[weight[j,2:(1+weight[j,1])]]),(1/kU/weight[j,1])^(0.5))
      ap<-ULikelihood(j,proposal)
    } else {
      proposal<-rnorm(1,U[j],sd=sigmaU)
      ap<-exp(-kU*sum((U[weight[j,2:(1+weight[j,1])]]-proposal)^2-(U[weight[j,2:(1+weight[j,1])]]-U[j])^2)/2)*ULikelihood(j,proposal)
    }
    un<-runif(1)
    if (un<=ap) {
      U[j]<<-proposal
      acceptU[1+i%%2]<<-acceptU[1+i%%2]+1
    } else {
      rejectU[1+i%%2]<<-rejectU[1+i%%2]+1
    }
  }
}

USample <- function(state) {
  state$fe <- state$fe + mean(U)
  U<<-U-mean(U) 
  cat(c(t(U),"\n"),file=file.path(params$outpath, "U.txt"),append=TRUE,sep=" ")
  cat(kU,file=file.path(params$outpath, "kU.txt"),append=TRUE,sep="\n")
  cat(c(t(acceptU[]/(acceptU[]+rejectU[])),"\n"),file=file.path(params$outpath, "acceptanceU.txt"),append=TRUE,sep=" ")
  acceptU<<-rep(0,2)
  rejectU<<-rep(0,2)   
  cat(USumFunction(U),"\n",file=file.path(params$outpath, "sumU.txt"),append=TRUE)
  return(state)
}

URisk <- function() {
  rep(U,each=params$tps)
}
ULikelihoodU <- function(j,proposal) {
  prod(dpois(cases[,j],n[j]*exp(fe+proposal))/dpois(cases[,j],n[j]*exp(fe+U[j])))
}
ULikelihoodRUX <- function(j,proposal) {
  prod(dpois(cases[,j],n[j]*exp(fe+R+proposal+betaX*X[,mbrg[j]]))/dpois(cases[,j],n[j]*exp(fe+R+U[j]+betaX*X[,mbrg[j]])))
}
ULikelihoodRUX2 <- function(j,proposal) {
  tps <- params$tps
  prod(dpois(cases[,j],n[j]*exp(fe+R+proposal+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]]))/dpois(cases[,j],n[j]*exp(fe+R+U[j]+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])))
}
ULikelihoodRUX3 <- function(j,proposal) {
  tps <- params$tps
  prod(dpois(cases[1,j],n[j]*exp(fe+R[1]+proposal+betaX[mbrg[j]]*X[1,mbrg[j]]))/dpois(cases[1,j],n[j]*exp(fe+R[1]+U[j]+betaX[mbrg[j]]*X[1,mbrg[j]])),dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+proposal+rep(betaX[mbrg[j]],tps-1)*(X[1:(tps-1),mbrg[j]]+X[2:tps,mbrg[j]])))/dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+U[j]+rep(betaX[mbrg[j]],tps-1)*(X[1:(tps-1),mbrg[j]]+X[2:tps,mbrg[j]]))))
}
ULikelihoodRUX4 <- function(j,proposal) {
  tps <- params$tps
  prod(dpois(cases[1,j],n[j]*exp(fe+R[1]+proposal+betaX[mbrg[j]]*X[1,mbrg[j]]))/dpois(cases[1,j],n[j]*exp(fe+R[1]+U[j]+betaX[mbrg[j]]*X[1,mbrg[j]])),dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+proposal+rep(betaX[mbrg[j]],tps-1)*pmax(X[1:(tps-1),mbrg[j]],X[2:tps,mbrg[j]])))/dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+U[j]+rep(betaX[mbrg[j]],tps-1)*pmax(X[1:(tps-1),mbrg[j]],X[2:tps,mbrg[j]]))))
}
ULikelihoodUX <- function(j,proposal) {
  prod(dpois(cases[,j],n[j]*exp(fe+proposal+betaX*X[,mbrg[j]]*logcases[,j]))/dpois(cases[,j],n[j]*exp(fe+U[j]+betaX*X[,mbrg[j]]*logcases[,j])))
}
ULikelihoodRU <- function(j,proposal) {
  prod(dpois(cases[,j],n[j]*exp(fe+R+proposal))/dpois(cases[,j],n[j]*exp(fe+R+U[j])))
}
ULikelihoodRUW <- function(j,proposal) {
  tps <- params$tps
  prod(dpois(cases[,j],n[j]*exp(fe+R+proposal+W[wthr[3:(tps+2),j]]+W[ws+wthr[2:(tps+1),j]]+W[2*ws+wthr[1:tps,j]]))/dpois(cases[,j],n[j]*exp(fe+R[]+U[j]+W[wthr[3:(tps+2),j]]+W[ws+wthr[2:(tps+1),j]]+W[2*ws+wthr[1:tps,j]])))
}

UConvergence <- function() {
  plotPairs("kU")
  plotPairs("sumU",p=F)
  plotPairs("acceptanceU",length(acceptU),F,ha=T)
}

UTraces <- function() {
  U<<-plotPairs("U",params$mbs,z=T)
}
