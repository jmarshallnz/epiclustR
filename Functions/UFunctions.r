# Defaults
sigmaU<-1

UInitialize <- function() {
  if (params$tidyup) {file.remove(file.path(params$outpath, "U.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "kU.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "acceptanceU.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "sumU.txt"))}

  state <- list(acceptU = matrix(0,2,1),
                rejectU = matrix(0,2,1),
                kU = 1,
                U = rnorm(params$mbs,0,1))
  return(state)
}

# computes the variance of the U's for the Gibbs gamma update
USumFunction <- function(U) {
  s<-0
  for (i in seq_along(U)) {
    s <- s + sum((U[i]-U[weight[i,2:(1+weight[i,1])]])^2)
  }
  s/2 # We divide by 2 here as the above sum will count all squared distances twice.
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

UUpdate <- function(i=0, state) {
  # save current state
  U <- state$U
  kU <- state$kU
  acceptU <- state$acceptU
  rejectU <- state$rejectU

  lenU <- length(U)
  # Gibb's step to update kU
  kU<-rgamma(1,aU+(lenU-1)/2,rate=(bU+USumFunction(U)/2))
  for (j in 1:lenU) {
    if (i%%2==0) {
      proposal<-rnorm(1,mean(U[weight[j,2:(1+weight[j,1])]]),(1/kU/weight[j,1])^(0.5))
      ap<-ULikelihood(j,U[j],proposal,state)
    } else {
      proposal<-rnorm(1,U[j],sd=sigmaU)
      prior_ratio <- -kU*sum((U[weight[j,2:(1+weight[j,1])]]-proposal)^2-(U[weight[j,2:(1+weight[j,1])]]-U[j])^2)/2
      ap<-ULikelihood(j,U[j],proposal,state)+prior_ratio
    }
    un<-runif(1)
    if (ap >= 0 || un<=exp(ap)) {
      U[j] <- proposal
      acceptU[1+i%%2] <- acceptU[1+i%%2]+1
    } else {
      rejectU[1+i%%2] <- rejectU[1+i%%2]+1
    }
  }
  state$kU <- kU
  state$fe <- state$fe + mean(U)
  state$U  <- U - mean(U)

  state$acceptU <- acceptU
  state$rejectU <- rejectU
  return(state)
}

USample <- function(state) {
  cat(c(t(state$U),"\n"),file=file.path(params$outpath, "U.txt"),append=TRUE,sep=" ")
  cat(state$kU,file=file.path(params$outpath, "kU.txt"),append=TRUE,sep="\n")
  cat(c(t(state$acceptU[]/(state$acceptU[]+state$rejectU[])),"\n"),file=file.path(params$outpath, "acceptanceU.txt"),append=TRUE,sep=" ")
  cat(USumFunction(state$U),"\n",file=file.path(params$outpath, "sumU.txt"),append=TRUE)
}

UInitAcceptance <- function(state) {
  state$acceptU <- rep(0,2)
  state$rejectU <- rep(0,2)
  return(state)
}

URisk <- function(state) {
  rep(state$U,each=length(state$R))
}

ULikelihoodU <- function(j,curr,prop,state) {
  fe <- state$fe
  sum(dpois(cases[,j],n[j]*exp(fe+prop), log=TRUE)-
      dpois(cases[,j],n[j]*exp(fe+curr), log=TRUE))
}
ULikelihoodRUX <- function(j,curr,prop,state) {
  fe <- state$fe
  R  <- state$R
  sum(dpois(cases[,j],n[j]*exp(fe+R+prop+betaX*X[,mbrg[j]]), log=TRUE)-
      dpois(cases[,j],n[j]*exp(fe+R+curr+betaX*X[,mbrg[j]]), log=TRUE))
}

u_likelihood_rux2 <- function(cases, n, fe, R, X, mbrg, betaX, curr, prop, j) {
  tps <- length(R)
  lambda_curr <- n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]])
  lambda_prop <- lambda_curr * exp(prop-curr)
  sum(cases[,j] * (prop - curr) - lambda_prop + lambda_curr)
}

ULikelihoodRUX2 <- function(j,curr,prop,state) {
  #  sum(dpois(cases[,j],n[j]*exp(fe+R+prop+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]]), log=TRUE)-
  #      dpois(cases[,j],n[j]*exp(fe+R+curr+rep(betaX[mbrg[j]],tps)*X[,mbrg[j]]), log=TRUE))
  u_likelihood_rux2(cases, n, state$fe, state$R, state$X, mbrg, state$betaX,
                    curr, prop, j)
}

ULikelihoodRUX3 <- function(j,curr,prop,state) {
  fe <- state$fe
  R  <- state$R
  tps <- length(R)
  sum(dpois(cases[1,j],n[j]*exp(fe+R[1]+prop+betaX[mbrg[j]]*X[1,mbrg[j]]), log=TRUE)-
      dpois(cases[1,j],n[j]*exp(fe+R[1]+curr+betaX[mbrg[j]]*X[1,mbrg[j]]), log=TRUE),
      dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+prop+rep(betaX[mbrg[j]],tps-1)*(X[1:(tps-1),mbrg[j]]+X[2:tps,mbrg[j]])), log=TRUE)-
      dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+curr+rep(betaX[mbrg[j]],tps-1)*(X[1:(tps-1),mbrg[j]]+X[2:tps,mbrg[j]])), log=TRUE))
}
ULikelihoodRUX4 <- function(j,curr,prop,state) {
  fe <- state$fe
  R  <- state$R
  tps <- length(R)
  sum(dpois(cases[1,j],n[j]*exp(fe+R[1]+prop+betaX[mbrg[j]]*X[1,mbrg[j]]), log=TRUE)-
      dpois(cases[1,j],n[j]*exp(fe+R[1]+U[j]+betaX[mbrg[j]]*X[1,mbrg[j]]), log=TRUE),
      dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+prop+rep(betaX[mbrg[j]],tps-1)*pmax(X[1:(tps-1),mbrg[j]],X[2:tps,mbrg[j]])), log=TRUE)-
      dpois(cases[2:tps,j],n[j]*exp(fe+R[2:tps]+curr+rep(betaX[mbrg[j]],tps-1)*pmax(X[1:(tps-1),mbrg[j]],X[2:tps,mbrg[j]])), log=TRUE))
}
ULikelihoodUX <- function(j,curr,prop,state) {
  fe <- state$fe
  # NOTE: This will fail, as logcases doesn't exist
  sum(dpois(cases[,j],n[j]*exp(fe+prop+betaX*X[,mbrg[j]]*logcases[,j]), log=TRUE)-
      dpois(cases[,j],n[j]*exp(fe+curr+betaX*X[,mbrg[j]]*logcases[,j]), log=TRUE))
}
ULikelihoodRU <- function(j,curr,prop,state) {
  fe <- state$fe
  R  <- state$R
  sum(dpois(cases[,j],n[j]*exp(fe+R+prop), log=TRUE)-
      dpois(cases[,j],n[j]*exp(fe+R+curr), log=TRUE))
}
ULikelihoodRUW <- function(j,curr,prop,state) {
  tps <- params$tps
  fe <- state$fe
  tps <- length(R)
  sum(dpois(cases[,j],n[j]*exp(fe+R+prop+W[wthr[3:(tps+2),j]]+W[ws+wthr[2:(tps+1),j]]+W[2*ws+wthr[1:tps,j]]), log=TRUE)-
      dpois(cases[,j],n[j]*exp(fe+R+curr+W[wthr[3:(tps+2),j]]+W[ws+wthr[2:(tps+1),j]]+W[2*ws+wthr[1:tps,j]]), log=TRUE))
}

UConvergence <- function(state) {
  plotPairs("kU")
  plotPairs("sumU",p=F)
  plotPairs("acceptanceU",length(state$acceptU),F,ha=T)
}

UTraces <- function() {
  U<<-plotPairs("U",params$mbs,z=T)
}
