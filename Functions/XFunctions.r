#
# Independent priors for X 
#
# Problems: 1) DIC should not be used.
#
# NB logcases has been removed from this version
# Xmode 0 = one betaX for everywhere
# Xmode 1 = one betaX, pX fixed
# Xmode 2 = betaX for each region
# Xmode 3 = betaX for each region, X covers two timepoints
# Xmode 4 = betaX for each region, X covers two timepoints but max(X_{t-1},X_t) is used.
# 
# Defaults
XOutput<-TRUE
XFullOutput<-FALSE
aX<-1
bX<-51 # On average one outbreak per year per region

XSetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aX<<-1
    bX<<-51
    abetaX<<-1
    bbetaX<<-1
  } else if (setpriors==2) {
    aX<<-1    # was 2
    bX<<-51   # was 102
    abetaX<<-1
    bbetaX<<-1
  } else if (setpriors==3) {
    aX<<-1
    bX<<-1
    abetaX<<-1
    bbetaX<<-1
  }     
}
XInitialise <- function() {
  tps <- params$tps
  mbs <- params$mbs
  Xmode <- params$Xmode

  mbrg<-scan(file.path(params$datapath,params$region,paste0("Regions",params$regionchoice,".txt")))
  mbrg<<-matrix(mbrg,1,mbs)
  rgs<<-max(mbrg)
  wch<<-list()
  lwch<<-matrix(0,1,rgs)
  for (j in 1:rgs) {
    wch[[j]]<<-which(mbrg==j)
    lwch[j]<<-length(wch[[j]])
  }
  pX <- 0.1
  X <- matrix(rbinom(tps*rgs,1,pX),tps,rgs)  
  if (Xmode==0) {
    betaX <- 0.2
    sigmaX<<-0.2
    XUpdate<<-XUpdate0
    XRisk<<-XRisk0
  } else if (Xmode==1) {
    betaX <- 0.2
    sigmaX<<-0.2
    
    XRisk<<-XRisk0
    XUpdate<<-XUpdate1
  } else if (Xmode==2) {
    if (params$tidyup) {file.remove(file.path(params$outpath, "betaXconditional.txt"))}
    betaX <- rep(0.2,rgs)
    sigmaX<<-1
    XUpdate<<-XUpdate2
    XRisk<<-XRisk2
  } else if (Xmode==3 | Xmode==4) {
    if (params$tidyup) {file.remove(file.path(params$outpath, "betaXconditional.txt"))}
    betaX <- rep(0.2,rgs)
    sigmaX<<-1
    abetaX<<-1
    bbetaX<<-1
    XUpdate<<-XUpdate3
    XRisk<<-XRisk3
    if (Xmode==4) {XRisk<-XRisk4}
  }

  acceptX <- 0
  rejectX <- 0
  if (params$tidyup) {
    file.remove(file.path(params$outpath, "betaX.txt"))
    file.remove(file.path(params$outpath, "pX.txt"))
    file.remove(file.path(params$outpath, "acceptanceX.txt"))
    file.remove(file.path(params$outpath, "smoothedCases.txt"))
    file.remove(file.path(params$outpath, "fullX.txt"))
  }
  state <- list(X = X,
                pX = pX,
                betaX = betaX,
                acceptX = acceptX,
                rejectX = rejectX)
  return(state)
}

XUpdate0 <- function(i=0, state) {
  pX <- state$pX
  betaX <- state$betaX
  acceptX <- state$acceptX
  rejectX <- state$rejectX

  tps <- nrow(state$X)
  rps <- ncol(state$X)

  # Update X
  state$X <- XLikelihood(state)
  #Update pX
  state$pX <- rbeta(1,aX+sum(state$X),bX+tps*rgs-sum(state$X)) 
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(betaX,proposal,state)
  un<-runif(1)
  if ((ap >= 0 || un<=exp(ap)) && proposal>0) {# This causes the prior for betaX to be the flat prior on the positive half-line.
    betaX<-proposal
    acceptX<-acceptX+1
  } else {
    rejectX<-rejectX+1
  }
  state$betaX <- betaX
  state$acceptX <- acceptX
  state$rejectX <- rejectX
  return(state)
}

XUpdate1 <- function(i=0, state) {
  pX <- state$pX
  betaX <- state$betaX
  acceptX <- state$acceptX
  rejectX <- state$rejectX

  tps <- nrow(state$X)
  rps <- ncol(state$X)

  # Update X
  state$X <- XLikelihood(state)
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(betaX,proposal,state)
  un<-runif(1)
  if ((ap >= 0 || un<=exp(ap)) && proposal>0) {# This causes the prior for betaX to be the flat prior on the positive half-line.
    betaX<-proposal
    acceptX<-acceptX+1
  } else {
    rejectX<-rejectX+1
  }
  state$pX <- pX
  state$betaX <- betaX
  state$acceptX <- acceptX
  state$rejectX <- rejectX
  return(state)
}
XUpdate2 <- function(i=0, state) {

  # Call straight into c-land
  state <- update_x(cases, n, wch, state, list(aX=aX, bX=bX, sigmaX=sigmaX, abetaX=abetaX, bbetaX=bbetaX))

  return(state)
}

XUpdate3 <- function(i=0, state) {
  pX <- state$pX
  betaX <- state$betaX
  acceptX <- state$acceptX
  rejectX <- state$rejectX

  tps <- nrow(state$X)
  rps <- ncol(state$X)

  # Update X
  for (j in 1:tps) {
    X$state[j,] <- XLikelihood(j,state)
  }
  #Update pX
  state$pX<-rbeta(1,aX+sum(state$X),bX+tps*rgs-sum(state$X))
  # Update betaX
  for (j in 1:rgs) {
    proposal<-rnorm(1,betaX[j],sigmaX)
    if (proposal < 0) {
      rejectX <- rejectX + 1
    } else {
      prior_ratio <- (abetaX - 1) * (log(proposal) - log(betaX[j])) - (proposal - betaX[j])*bbetaX
      ap <- betaXLikelihood(j,betaX[j],proposal,state)+prior_ratio
      un<-runif(1)
      if (ap >= 0 || un<=exp(ap)) {
        betaX[j]<-proposal
        acceptX<-acceptX+1
      } else {
        rejectX<-rejectX+1
      }
    }
  }
  state$betaX <- betaX
  state$acceptX <- acceptX
  state$rejectX <- rejectX
  return(state)
}

squashProd <- function(input) {# This may be faster when applied to a 3d ragged array tps*mbs*max(mb in a rg)
  tps <- params$tps # If X is already created, these are known?
  mbs <- params$mbs

  input<-matrix(input,tps,mbs)
  output<-matrix(0,tps,rgs)
  for (k in 1:rgs)
  {
    if (length(wch[[k]]) > 1)
    {
      output[,k]<-apply(input[,wch[[k]]],1,sum)
    }
    else
    {
      output[,k]<-input[,wch[[k]]]
    }
  }
  return(exp(output))
}

squashProd3 <- function(input) {

  output<-matrix(0,1,rgs)
  for (k in 1:rgs) {
    output[k]<-sum(input[wch[[k]]])
  }
  return(exp(output))
}

XRisk0 <- function(state) {
  state$X[,mbrg]*state$betaX
}
XRisk2 <- function(state) {
  state$X[,mbrg]*rep(state$betaX[mbrg],each=nrow(state$X))
}
XRisk3 <- function(state) {
  X <- state$X
  betaX <- state$betaX
  rbind(X[1,mbrg],X[1:(nrow(X)-1),mbrg]+X[2:nrow(X),mbrg])*rep(betaX[mbrg],each=nrow(X))
}
XRisk4 <- function(state) {
  X <- state$X
  betaX <- state$betaX
  rbind(X[1,mbrg],pmax(X[1:(nrow(X)-1),mbrg],X[2:nrow(X),mbrg]))*rep(betaX[mbrg],each=nrow(X))
}
XLikelihoodRUX <- function(state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X
  Xr <- squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+0*betaX), log=TRUE)-
                   dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+1*betaX), log=TRUE))
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}

XLikelihoodRUX2 <- function(state) {
  x_sample_rux2(cases, n, state$fe, state$R, state$U, state$betaX, state$pX, wch)
}

XLikelihoodRUX3 <- function(j,state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X

  # TODO: Figure out what each case is doing
  if (j>1 && j<nrow(X)) {
    Xr <- squashProd3((dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+0)*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(0+X[j+1,mbrg])*betaX[mbrg]), log=TRUE))-
                      (dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+1)*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(1+X[j+1,mbrg])*betaX[mbrg]), log=TRUE)))
  } else if (j==1) {
    Xr <- squashProd3((dpois(cases[j,],n*exp(fe+R[j]+U+0*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(0+X[j+1,mbrg])*betaX[mbrg]), log=TRUE))-
                      (dpois(cases[j,],n*exp(fe+R[j]+U+1*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(1+X[j+1,mbrg])*betaX[mbrg]), log=TRUE)))
  } else { # if j == nrow(X)
    Xr <- squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+0)*betaX[mbrg]), log=TRUE)-
                      dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+1)*betaX[mbrg]), log=TRUE))
  }
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}
XLikelihoodRUX4 <- function(j,state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X
  if (j>1 && j<nrow(X)) {
    Xr <- squashProd3((dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],0)*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(0,X[j+1,mbrg])*betaX[mbrg]), log=TRUE))-
                      (dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],1)*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(1,X[j+1,mbrg])*betaX[mbrg]), log=TRUE)))
  } else if (j==1) {
    Xr <- squashProd3((dpois(cases[j,],n*exp(fe+R[j]+U+0*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(0,X[j+1,mbrg])*betaX[mbrg]), log=TRUE))-
                      (dpois(cases[j,],n*exp(fe+R[j]+U+1*betaX[mbrg]), log=TRUE)
                      +dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(1,X[j+1,mbrg])*betaX[mbrg]), log=TRUE)))
  } else { # if j == nrow(X)
    Xr <- squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],0)*betaX[mbrg]), log=TRUE)-
                      dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],1)*betaX[mbrg]), log=TRUE))
  }
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}
XLikelihoodRX <- function(state) {
  fe  <- state$fe
  R   <- state$R
  betaX <- state$betaX
  Xr <- squashProd(dpois(cases,rep(n,each=length(R))*exp(fe+rep(R,ncol(n))+0*betaX), log=TRUE)-
                   dpois(cases,rep(n,each=length(R))*exp(fe+rep(R,ncol(n))+1*betaX), log=TRUE))
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}
XLikelihoodUX <- function(state) {
  tps <- params$tps
  fe  <- state$fe
  U   <- state$U
  betaX <- state$betaX
  Xr <- squashProd(dpois(cases,rep(n,each=tps)*exp(fe+rep(U,each=tps)+0*betaX), log=TRUE)-
                   dpois(cases,rep(n,each=tps)*exp(fe+rep(U,each=tps)+1*betaX), log=TRUE))
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}
XLikelihoodX <- function(state) {
  tps <- params$tps
  fe  <- state$fe
  betaX <- state$betaX
  Xr <- squashProd(dpois(cases,rep(n,each=tps)*exp(fe+0*betaX), log=TRUE)-
                   dpois(cases,rep(n,each=tps)*exp(fe+1*betaX), log=TRUE))
  Xfcd <- state$pX / (Xr*(1-state$pX) + state$pX)
  matrix(rbinom(length(Xfcd),1,Xfcd),nrow(Xfcd),ncol(Xfcd))
}
betaXLikelihoodRUX <- function(curr,prop,state) {
  tps <- length(R)
  mbs <- ncol(n)
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  X   <- state$X
  sum(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*prop), log=TRUE)-
      dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*curr), log=TRUE))
}

betaXLikelihoodRUX2 <- function(j,curr,prop,state) {
  betax_likelihood_rux2(cases, n, state$fe, state$R, state$U, state$X, wch[[j]],
                        curr, prop, j)
}

betaXLikelihoodRUX3 <- function(j,curr,prop,state) {
  tps <- length(R)
  mbs <- ncol(n)
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  X   <- state$X
  sum(dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*prop), log=TRUE)-
      dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*curr), log=TRUE),
      dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*prop), log=TRUE)-
      dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*curr), log=TRUE))
}
betaXLikelihoodRUX4 <- function(j,curr,prop,state) {
  tps <- length(R)
  mbs <- ncol(n)
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  X   <- state$X
  sum(dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*prop), log=TRUE)-
      dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*curr), log=TRUE),
      dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+pmax(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])],X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*prop), log=TRUE)-
      dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+pmax(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])],X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*curr), log=TRUE))
}

XSample <- function(state) {
  cat(state$pX,"\n",file=file.path(params$outpath, "pX.txt"),append=TRUE)
  cat(state$betaX,"\n",file=file.path(params$outpath, "betaX.txt"),append=TRUE)
  if (params$Xmode>1) {
    cat((1-apply(1-state$X,2,prod))*state$betaX,"\n",file=file.path(params$outpath, "betaXconditional.txt"),append=TRUE)
  }
  cat(state$acceptX/(state$acceptX+state$rejectX),"\n",file=file.path(params$outpath, "acceptanceX.txt"),append=TRUE)
  if (XOutput) {
    cat(state$X,"\n",file=file.path(params$outpath, "X.txt"))
  }
  if (XFullOutput) {
    out <- character(ncol(state$X))
    for (i in 1:ncol(state$X)) {
      out[i] <- paste(state$X[,i], collapse="")
    }
    cat(out,"\n",file=file.path(params$outpath, "fullX.txt"),append=TRUE)
  }
}

XInitAcceptance <- function(state) {
  state$acceptX<-0
  state$rejectX<-0
  return(state)
}

XConvergence <- function(state) {
  rgs <- ncol(state$X)
  pX<<-plotPairs("pX")
  if (params$Xmode>1) {
    betaX<<-plotPairs("betaX",rgs,use=T)
  } else {
    betaX<<-plotPairs("betaX",use=T)
    # trick for returning betaX to a vector
    betaX<<-sum(betaX)
  }
  plotPairs("acceptanceX",length(state$acceptX),F,half=T)
  if (params$Xmode>1) {
    try(input<-scan(file.path(params$outpath, "betaXconditional.txt")),T)
    input<-matrix(input,rgs,length(input)/rgs)
    rr<-matrix(0,rgs,1)
    for (k in 1:rgs) {
      w<-which(input[k,]>0)
      w<-w[which(w>=1+params$burnin/params$samplefreq)]
      if (length(w)==0) {
        betaX[k]<<-0
      } else {
        betaX[k]<<-mean(input[k,w])
        rr[k]<-mean(exp(input[k,w]))
      }
    }
    write.table(t(betaX),file.path(params$outpath, "posteriorbetaXconditional.txt"),row.names=F,col.names=F)
    write.table(t(rr),file.path(params$outpath, "relativeriskbetaXconditional.txt"),row.names=F,col.names=F)
  }
}
XContinue <- function() {
  if (params$Xmode<2) {
    betaX<<-Upload("betaX",1)
  } else {
    betaX<<-Upload("betaX",rgs)
  }
  pX<<-Upload("pX",1)
  X<<-scan("X.txt")
  X<<-matrix(X,params$tps,params$rgs)
}
