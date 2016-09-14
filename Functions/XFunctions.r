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
pX<-0.1
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
  X<<-matrix(rbinom(tps*rgs,1,pX),tps,rgs)  
  cumX<<-matrix(0,tps,rgs)
  if (Xmode==0) {
    betaX<<-0.2
    sigmaX<<-0.2
    XUpdate<<-XUpdate0
  } else if (Xmode==1) {
    betaX<<-0.2
    sigmaX<<-0.2
    
    XUpdate<<-XUpdate1
  } else if (Xmode==2) {
    if (params$tidyup) {file.remove(file.path(params$outpath, "betaXconditional.txt"))}
    betaX<<-matrix(0.2,1,rgs)
    sigmaX<<-1
    XUpdate<<-XUpdate2
    XRisk<<-XRisk2
  } else if (Xmode==3 | Xmode==4) {
    if (params$tidyup) {file.remove(file.path(params$outpath, "betaXconditional.txt"))}
    betaX<<-matrix(0.2,1,rgs)
    sigmaX<<-1
    abetaX<<-1
    bbetaX<<-1
    XUpdate<<-XUpdate3
    XRisk<<-XRisk3
    if (Xmode==4) {XRisk<-XRisk4}
  }

  acceptX<<-0
  rejectX<<-0
  if (params$tidyup) {
    file.remove(file.path(params$outpath, "betaX.txt"))
    file.remove(file.path(params$outpath, "pX.txt"))
    file.remove(file.path(params$outpath, "acceptanceX.txt"))
    file.remove(file.path(params$outpath, "smoothedCases.txt"))
    file.remove(file.path(params$outpath, "fullX.txt"))
  }
}
XUpdate0 <- function(i=0) {
  tps <- nrow(X)
  rps <- ncol(X)

  # Update X
  Xfcd<-XLikelihood(1)*pX/(XLikelihood(0)*(1-pX)+XLikelihood(1)*pX)
  X<<-matrix(rbinom(tps*rgs,1,Xfcd),tps,rgs)
  #Update pX
  pX<<-rbeta(1,aX+sum(X),bX+tps*rgs-sum(X)) 
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(proposal)
  un<-runif(1)
  if (un<=ap && proposal>0) {# This causes the prior for betaX to be the flat prior on the positive half-line.
    betaX<<-proposal
    acceptX<<-acceptX+1
  } else {
    rejectX<<-rejectX+1
  }   
}
XUpdate1 <- function(i=0) {
  tps <- nrow(X)
  rps <- ncol(X)

  # Update X
  Xfcd<-XLikelihood(1)*pX/(XLikelihood(0)*(1-pX)+XLikelihood(1)*pX)
  X<<-matrix(rbinom(tps*rgs,1,Xfcd),tps,rgs)
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(proposal)
  un<-runif(1)
  if (un<=ap && proposal>0) {# This causes the prior for betaX to be the flat prior on the positive half-line.
    betaX<<-proposal
    acceptX<<-acceptX+1
  } else {
    rejectX<<-rejectX+1
  }   
}
XUpdate2 <- function(i=0) {
  tps <- nrow(X)
  rps <- ncol(X)

  # Update X
  Xfcd<-XLikelihood(1)*pX/(XLikelihood(0)*(1-pX)+XLikelihood(1)*pX)
  X<<-matrix(rbinom(tps*rgs,1,Xfcd),tps,rgs)
  #Update pX
  pX<<-rbeta(1,aX+sum(X),bX+tps*rgs-sum(X)) 
  # Update betaX
  for (j in 1:rgs) {
    proposal<-rnorm(1,betaX[j],sigmaX)
    ap<-betaXLikelihood(proposal,j)*(dgamma(proposal,abetaX,bbetaX)/dgamma(betaX[j],abetaX,bbetaX))
    un<-runif(1)
    if (un<=ap) {
      betaX[j]<<-proposal
      acceptX<<-acceptX+1
    } else {
      rejectX<<-rejectX+1
    }
  }   
}
XUpdate3 <- function(i=0) {
  tps <- nrow(X)
  rps <- ncol(X)

  # Update X
  for (j in 1:tps) {
    X1<-XLikelihood(j,1)*pX
    Xfcd<-X1/(XLikelihood(j,0)*(1-pX)+X1)
    X[j,]<<-rbinom(rgs,1,Xfcd)
  }
  #Update pX
  pX<<-rbeta(1,aX+sum(X),bX+tps*rgs-sum(X))
  # Update betaX
  for (j in 1:rgs) {
    proposal<-rnorm(1,betaX[j],sigmaX)
    ap<-betaXLikelihood(proposal,j)*(dgamma(proposal,abetaX,bbetaX)/dgamma(betaX[j],abetaX,bbetaX))
    un<-runif(1)
    if (un<=ap) {
      betaX[j]<<-proposal
      acceptX<<-acceptX+1
    } else {
      rejectX<<-rejectX+1
    }
  }   
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
XRisk <- function() {
  X[,mbrg]*betaX
}
XRisk2 <- function() {
  X[,mbrg]*rep(betaX[mbrg],each=nrow(X))
}
XRisk3 <- function() {
  rbind(X[1,mbrg],X[1:(nrow(X)-1),mbrg]+X[2:nrow(X),mbrg])*rep(betaX[mbrg],each=nrow(X))
}
XRisk4 <- function() {
  rbind(X[1,mbrg],pmax(X[1:(nrow(X)-1),mbrg],X[2:nrow(X),mbrg]))*rep(betaX[mbrg],each=nrow(X))
}
XLikelihoodRUX <- function(x) {
  mbs <- params$mbs
  squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,mbs)+rep(U,each=nrow(X))+x*betaX)))
}
XLikelihoodRUX2 <- function(x) {
  squashProd(dpois(cases,rep(n,each=nrow(X))*exp(fe+rep(R,ncol(n))+rep(U,each=nrow(X))+x*rep(betaX[mbrg],each=nrow(X)))))
}
XLikelihoodRUX3 <- function(j,x) {
  # TODO: Figure out what each case is doing
  if (j>1 && j<nrow(X)) {
    return(squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+(X[j-1,mbrg]+x)*betaX[mbrg]))*dpois(cases[j+1,],n*exp(fe+R[j+1]+U+(x+X[j+1,mbrg])*betaX[mbrg]))))
  } else if (j==1) {
    return(squashProd3(dpois(cases[1,],n*exp(fe+R[1]+U+x*betaX[mbrg]))*dpois(cases[2,],n*exp(fe+R[2]+U+(x+X[2,mbrg])*betaX[mbrg]))))
  } else { # if j == nrow(X)
    return(squashProd3(dpois(cases[nrow(X),],n*exp(fe+R[nrow(X)]+U+(X[nrow(X)-1,mbrg]+x)*betaX[mbrg]))))
  }
}
XLikelihoodRUX4 <- function(j,x) {
  if (j>1 && j<nrow(X)) {
    return(squashProd3(dpois(cases[j,],n*exp(fe+R[j]+U+pmax(X[j-1,mbrg],x)*betaX[mbrg]))*dpois(cases[j+1,],n*exp(fe+R[j+1]+U+pmax(x,X[j+1,mbrg])*betaX[mbrg]))))
  } else if (j==1) {
    return(squashProd3(dpois(cases[1,],n*exp(fe+R[1]+U+x*betaX[mbrg]))*dpois(cases[2,],n*exp(fe+R[2]+U+pmax(x,X[2,mbrg])*betaX[mbrg]))))
  } else { # if j == nrow(X)
    return(squashProd3(dpois(cases[nrow(X),],n*exp(fe+R[nrow(X)]+U+pmax(X[nrow(X)-1,mbrg],x)*betaX[mbrg]))))
  }
}
XLikelihoodRX <- function(x) {
  squashProd(dpois(cases,rep(n,each=length(R))*exp(fe+rep(R,ncol(n))+x*betaX)))
}
XLikelihoodUX <- function(x) {
  tps <- params$tps
  squashProd(dpois(cases,rep(n,each=tps)*exp(fe+rep(U,each=tps)+x*betaX)))
}
XLikelihoodX <- function(x) {
  tps <- params$tps
  squashProd(dpois(cases,rep(n,each=tps)*exp(fe+x*betaX)))
}
betaXLikelihoodRUX <- function(proposal) {
  tps <- length(R)
  mbs <- ncol(n)
  prod(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*proposal))/dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*betaX)))
}
betaXLikelihoodRUX2 <- function(proposal,j) {
  tps <- length(R)
  mbs <- ncol(n)
  prod(dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*proposal))/dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*betaX[j])))
}
betaXLikelihoodRUX3 <- function(proposal,j) {
  tps <- length(R)
  mbs <- ncol(n)
  prod(dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*proposal))/dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*betaX[j])),dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*proposal))/dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])]+X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*betaX[j])))
}
betaXLikelihoodRUX4 <- function(proposal,j) {
  tps <- length(R)
  mbs <- ncol(n)
  prod(dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*proposal))/dpois(cases[1,wch[[j]]],n[wch[[j]]]*exp(fe+rep(R[1],lwch[j])+U[wch[[j]]]+X[1,j]*betaX[j])),dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+pmax(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])],X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*proposal))/dpois(cases[2:tps,wch[[j]]],rep(n[wch[[j]]],each=tps-1)*exp(fe+rep(R[2:tps],lwch[j])+rep(U[wch[[j]]],each=tps-1)+pmax(X[rep((j-1)*tps,each=tps-1)+rep(2:tps,lwch[j])],X[rep((j-1)*tps,each=tps-1)+rep(1:(tps-1),lwch[j])])*betaX[j])))
}
XSample <- function(i) {
  # TODO: Remove this - it is completely ludicrous...
  if (i>params$burnin) {
    cumX<<-cumX+X
  }
  cat(pX,"\n",file=file.path(params$outpath, "pX.txt"),append=TRUE)
  cat(betaX,"\n",file=file.path(params$outpath, "betaX.txt"),append=TRUE)
  if (params$Xmode>1) {
    cat((1-apply(1-X,2,prod))*betaX,"\n",file=file.path(params$outpath, "betaXconditional.txt"),append=TRUE)
  }
  cat(acceptX/(acceptX+rejectX),"\n",file=file.path(params$outpath, "acceptanceX.txt"),append=TRUE)
  if (XOutput) {
    cat(X,"\n",file=file.path(params$outpath, "X.txt"))
    cat(cumX,"\n",file=file.path(params$outpath, "cumulativeX.txt"))
  }
  if (XFullOutput) {
    cat(X,"\n",file=file.path(params$outpath, "fullX.txt"),append=TRUE)
  } 
}

XInitAcceptance <- function(state) {
  acceptX<<-0
  rejectX<<-0
  return(state)
}

XConvergence <- function() {
  rgs <- ncol(X)
  pX<<-plotPairs("pX")
  if (params$Xmode>1) {
    betaX<<-plotPairs("betaX",rgs,use=T)
  } else {
    betaX<<-plotPairs("betaX",use=T)
    # trick for returning betaX to a vector
    betaX<<-sum(betaX)
  }
  plotPairs("acceptanceX",length(acceptX),F,half=T)
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
