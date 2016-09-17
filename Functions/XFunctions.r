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
    aX<<-2
    bX<<-102
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
  cumX <- matrix(0,tps,rgs)
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
    betaX <- matrix(0.2,1,rgs)
    sigmaX<<-1
    XUpdate<<-XUpdate2
    XRisk<<-XRisk2
  } else if (Xmode==3 | Xmode==4) {
    if (params$tidyup) {file.remove(file.path(params$outpath, "betaXconditional.txt"))}
    betaX <- matrix(0.2,1,rgs)
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
                cumX = cumX,
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
  Xfcd<-XLikelihood(1, state)*pX/(XLikelihood(0, state)*(1-pX)+XLikelihood(1, state)*pX)
  state$X <- matrix(rbinom(tps*rgs,1,Xfcd),tps,rgs)
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
  # TODO: Ideally we'd remove this. Problem is full X is large to save
  #       in terms of the posterior, and if we want it right we have
  #       to get rid of the burnin period.
  #       I guess saving X to a binary file might be way more efficient
  #       though? About 16 times smaller.
  if (i>params$burnin) {
    state$cumX<-state$cumX+X
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
  Xfcd<-XLikelihood(1,state)*pX/(XLikelihood(0,state)*(1-pX)+XLikelihood(1,state)*pX)
  state$X<-matrix(rbinom(tps*rgs,1,Xfcd),tps,rgs)
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
  # TODO: Ideally we'd remove this. Problem is full X is large to save
  #       in terms of the posterior, and if we want it right we have
  #       to get rid of the burnin period.
  #       I guess saving X to a binary file might be way more efficient
  #       though? About 16 times smaller.
  if (i>params$burnin) {
    state$cumX<-state$cumX+X
  }
  state$pX <- pX
  state$betaX <- betaX
  state$acceptX <- acceptX
  state$rejectX <- rejectX
  return(state)
}
XUpdate2 <- function(i=0, state) {
  pX <- state$pX
  betaX <- state$betaX
  acceptX <- state$acceptX
  rejectX <- state$rejectX

  tps <- nrow(state$X)
  rps <- ncol(state$X)

  # Update X, pX
  X1 <- XLikelihood(1,state)*pX
  X0 <- XLikelihood(0,state)*(1-pX)
  state$X<-matrix(rbinom(tps*rgs,1,X1/(X0 + X1)),tps,rgs)
  state$pX<-rbeta(1,aX+sum(state$X),bX+tps*rgs-sum(state$X))
  # Update betaX
  for (j in 1:rgs) {
    proposal <- rnorm(1,betaX[j],sigmaX)
    un<-runif(1) # always take the same number of random numbers...
    if (proposal <= 0) {
      rejectX <- rejectX + 1
    } else {
      prior_ratio <- (abetaX - 1) * (log(proposal) - log(betaX[j])) - (proposal - betaX[j])*bbetaX
      ap <- betaXLikelihood(j,betaX[j],proposal,state)+prior_ratio
      if (ap >= 0 || un<=exp(ap)) {
        betaX[j]<-proposal
        acceptX<-acceptX+1
      } else {
        rejectX<-rejectX+1
      }
    }
  }
  # TODO: Ideally we'd remove this. Problem is full X is large to save
  #       in terms of the posterior, and if we want it right we have
  #       to get rid of the burnin period.
  #       I guess saving X to a binary file might be way more efficient
  #       though? About 16 times smaller.
  if (i>params$burnin) {
    state$cumX<-state$cumX+state$X
  }
  state$betaX <- betaX
  state$acceptX <- acceptX
  state$rejectX <- rejectX
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
    X1<-XLikelihood(j,1,state)*pX
    Xfcd<-X1/(XLikelihood(j,0,state)*(1-pX)+X1)
    X$state[j,]<-rbinom(rgs,1,Xfcd)
  }
  #Update pX
  state$pX<-rbeta(1,aX+sum(state$X),bX+tps*rgs-sum(state$X))
  # Update betaX
  for (j in 1:rgs) {
    proposal<-rnorm(1,betaX[j],sigmaX)
    un<-runif(1) # always take the same number of random numbers...
    if (proposal < 0) {
      rejectX <- rejectX + 1
    } else {
      prior_ratio <- (abetaX - 1) * (log(proposal) - log(betaX[j])) - (proposal - betaX[j])*bbetaX
      ap <- betaXLikelihood(j,betaX[j],proposal,state)+prior_ratio
      if (ap >= 0 || un<=exp(ap)) {
        betaX[j]<-proposal
        acceptX<-acceptX+1
      } else {
        rejectX<-rejectX+1
      }
    }
  }
  # TODO: Ideally we'd remove this. Problem is full X is large to save
  #       in terms of the posterior, and if we want it right we have
  #       to get rid of the burnin period.
  #       I guess saving X to a binary file might be way more efficient
  #       though? About 16 times smaller.
  if (i>params$burnin) {
    state$cumX<-state$cumX+X
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
      output[,k]<-apply(input[,wch[[k]]],1,prod)
    }
    else
    {
      output[,k]<-input[,wch[[k]]]
    }
  }
  return(output)
}

squashProd3 <- function(input) {

  output<-matrix(0,1,rgs)
  for (k in 1:rgs) {
    output[k]<-prod(input[wch[[k]]])
  }
  return(output)
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
XLikelihoodRUX <- function(x, state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X
  squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+x*betaX)))
}
XLikelihoodRUX2 <- function(x, state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X
  squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+x*rep(betaX[mbrg],each=nrow(X)))))
}
XLikelihoodRUX3 <- function(j,x,state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X

  # TODO: Figure out what each case is doing
  if (j>1 && j<nrow(X)) {
    return(squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+x)*betaX[mbrg]))*dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(x+X[j+1,mbrg])*betaX[mbrg]))))
  } else if (j==1) {
    return(squashProd3(dpois(cases[1,],n*exp(fe+R[1]+U+x*betaX[mbrg]))*dpois(cases[2,],n*exp(fe+R[2]+U+(x+X[2,mbrg])*betaX[mbrg]))))
  } else { # if j == nrow(X)
    return(squashProd3(dpois(cases[nrow(X),],n*exp(fe+R[nrow(X)]+U+(X[nrow(X)-1,mbrg]+x)*betaX[mbrg]))))
  }
}
XLikelihoodRUX4 <- function(j,x,state) {
  fe  <- state$fe
  R   <- state$R
  U   <- state$U
  betaX <- state$betaX
  X     <- state$X
  if (j>1 && j<nrow(X)) {
    return(squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],x)*betaX[mbrg]))*dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(x,X[j+1,mbrg])*betaX[mbrg]))))
  } else if (j==1) {
    return(squashProd3(dpois(cases[1,],n*exp(fe+R[1]+U+x*betaX[mbrg]))*dpois(cases[2,],n*exp(fe+R[2]+U+pmax(x,X[2,mbrg])*betaX[mbrg]))))
  } else { # if j == nrow(X)
    return(squashProd3(dpois(cases[nrow(X),],n*exp(fe+R[nrow(X)]+U+pmax(X[nrow(X)-1,mbrg],x)*betaX[mbrg]))))
  }
}
XLikelihoodRX <- function(x,state) {
  fe  <- state$fe
  R   <- state$R
  betaX <- state$betaX
  squashProd(dpois(cases,rep(n,each=length(R))*exp(fe+rep(R,ncol(n))+x*betaX)))
}
XLikelihoodUX <- function(x,state) {
  tps <- params$tps
  fe  <- state$fe
  U   <- state$U
  betaX <- state$betaX
  squashProd(dpois(cases,rep(n,each=tps)*exp(fe+rep(U,each=tps)+x*betaX)))
}
XLikelihoodX <- function(x,state) {
  tps <- params$tps
  fe  <- state$fe
  betaX <- state$betaX
  squashProd(dpois(cases,rep(n,each=tps)*exp(fe+x*betaX)))
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

betax_likelihood_rux2 <- function(cases, n, fe, R, U, X, rgmb_j, curr, prop, j) {
  tps <- length(R)
  Xj <- X[rep(j-1,each=tps)*tps+rep(1:tps,length(rgmb_j))]
  lambda_curr <- rep(n[rgmb_j],each=tps)*exp(fe+rep(R,length(rgmb_j))+rep(U[rgmb_j],each=tps)+Xj*curr)
  lambda_prop <- lambda_curr * exp(Xj*(prop-curr))
  #  sum(dpois(cases[,wch[[j]]],lambda_prop, log=TRUE) -
  #      dpois(cases[,wch[[j]]],lambda_curr, log=TRUE))
  sum(cases[,rgmb_j] * Xj * (prop - curr) - lambda_prop + lambda_curr)
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
    cat(state$cumX,"\n",file=file.path(params$outpath, "cumulativeX.txt"))
  }
  if (XFullOutput) {
    cat(state$X,"\n",file=file.path(params$outpath, "fullX.txt"),append=TRUE)
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
XAnalysis <- function() {
  tps <- nrow(X)
  rgs <- ncol(X)

  plot(1:tps,t="n",xlab="week",ylab="X",main="Epidemic Analysis",xlim=c(1,tps),ylim=c(0,1))
  for (i in 1:rgs) {
    lines(1:tps,cumX[,i]/(iters-burnin)*samplefreq)
  }
  for (i in 1:rgs) {
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]==1)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)==1)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="darkgreen",pch=20,cex=0.75)
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]==2)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)==2)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="orange",pch=20,cex=0.75)
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]>2)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)>2)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="darkred",pch=20,cex=0.75)
  }
  X<<-cumX/(iters-burnin)*samplefreq
  write.table(t(X),file.path(params$outpath, "posteriorX.txt"),row.names=F,col.names=F)
}
XAnalysis2 <- function() {
  tps <- nrow(X)
  rgs <- ncol(X)

  par(mfrow=c(3,4))
  for (i in 1:rgs) {
    plot(1:tps,cumX[,i]/(iters-burnin)*samplefreq,t="l",xlab="week",ylab="X",main=paste("X",i),xlim=c(1,tps),ylim=c(0,1))
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]==1)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)==1)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="darkgreen",pch=20,cex=0.75)
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]==2)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)==2)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="orange",pch=20,cex=0.75)
    if (lwch[i]==1) {
      casealert<-which(cases[,wch[[i]]]>2)
    } else {
      casealert<-which(apply(cases[,wch[[i]]],1,sum)>2)
    }
    points(casealert,cumX[casealert,i]/(iters-burnin)*samplefreq,col="darkred",pch=20,cex=0.75)
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
  cumX<<-scan("cumulativeX.txt")
  cumX<<-matrix(cumX,params$tps,params$rgs)
}
