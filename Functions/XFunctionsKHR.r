# Temporally correlated epidemic indicators, with forward filtering backward sampling updates.
# 
# Problems: 1) kX does not converge with logcases.
#           2) It is quite slow.
#           3) Sometimes all the X's go to one and fe -> -infty.
#
# Xmode 0 = one betaX for everywhere with gamma prior
# Xmode 1 = one betaX for everywhere with flat prior
# Xmode 2 = betaX for each region
#
# Defaults
XOutput<-TRUE
XFullOutput<-FALSE
kX<-c(0.05,0.5)
acceptX<-0
rejectX<-0
if (tidyup) {
  file.remove("betaX.txt")
  file.remove("kX.txt")
  file.remove("sumX.txt")
  file.remove("acceptanceX.txt")
  file.remove("smoothedCases.txt")
  file.remove("fullX.txt")
}
XSetPriors <- function(setpriors) {
  if (setpriors==0) {
    aX<<-c(1,1)
    bX<<-c(1,1)
    abetaX<<-1
    bbetaX<<-1
  } else if (setpriors==1) {
    aX<<-c(0.5,0.5)
    bX<<-c(0.5,0.5)
  } else if (setpriors==2) {
    aX<<-c(1,2)
    bX<<-c(51,2)
    abetaX<<-1
    bbetaX<<-1
  } else if (setpriors==3) {
    aX<<-c(1,1)
    bX<<-c(1,1)
    abetaX<<-1
    bbetaX<<-1
  }     
}
XInitialise <- function() {
  mbrg<-scan(paste(datapath,region,"\\Regions",regionchoice,".txt",sep=""))
  mbrg<<-matrix(mbrg,1,mbs)
  rgs<<-max(mbrg)
  wch<<-list()
  lwch<<-matrix(0,1,rgs)
  for (j in 1:rgs) {
    wch[[j]]<<-which(mbrg==j)
    lwch[j]<<-length(wch[[j]])
  }
  X<<-matrix(rbinom(tps*rgs,1,1/52),tps,rgs)  
  cumX<<-matrix(0,tps,rgs)
  if (Xmode==0) {
    betaX<<-0.2
    sigmaX<<-0.1
    XUpdate<<-XUpdate0
  } else if (Xmode==1) {
    betaX<<-0.2
    sigmaX<<-0.1
    XUpdate<<-XUpdate1
  } else if (Xmode==2) {
    if (tidyup) {file.remove("betaXconditional.txt")}
    betaX<<-matrix(0.2,1,rgs)
    sigmaX<<-1
    XUpdate<<-XUpdate2
    XRisk<<-XRisk2
  }
}
XJumpFunction<- function(X) {
  temp<-X[1:(tps-1),]+2*X[2:tps,]
  c(length(which(temp==0)),length(which(temp==1)),length(which(temp==2)),length(which(temp==3)))
}
XUpdate0 <- function(i=0) {
  # Update X
  pOne<-matrix(0,tps,rgs)#p(X[i,j]=1|Y_1,...,Y_{t-1})
  p<-matrix(0,tps,rgs)   #p(X[i,j]=1|Y_1,...,Y_t)
  pOne[1,]<-rep(kX[1]/sum(kX),rgs)
  for (j in 1:tps) {
    if (j>1) {pOne[j,]<-kX[1]*(1-p[j-1,])+(1-kX[2])*p[j-1,]}
    pStarZero<-XLikelihood(j,0)*(1-pOne[j,])
    pStarOne <-XLikelihood(j,1)*pOne[j,]
    p[j,]<-pStarOne/(pStarOne+pStarZero)
  }
  for (j in tps:1) {
    if (j==tps) {
      ap<-p[j,]
    } else {
      # ap = p(X[t+1]|X[t]=1)*p(X[t]=1|Y1..t)/p(X[t+1]|Y1..t)
      ap<-(X[j+1,]+(-1)^X[j+1,]*kX[2])*p[j,]/(1-X[j+1,]+(-1)^(1-X[j+1,])*pOne[j+1,])
    }
    un<-runif(rgs)
    X[j,which(un<=ap)]<<-1
    X[j,which(un>ap)]<<-0
  }
  # Update kX
  tm<-XJumpFunction(X)
  proposal<-c(rbeta(1,aX[1]+tm[3],bX[1]+tm[1]),rbeta(1,aX[2]+tm[2],bX[2]+tm[4]))# (p01,p10)
  ap<-prod(dbinom(X[1,],1,proposal[1]/sum(proposal))/dbinom(X[1,],1,kX[1]/sum(kX)))
  un<-runif(1)
  if (un<=ap) {kX<<-proposal}
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(proposal)*(dgamma(proposal,abetaX,bbetaX)/dgamma(betaX[j],abetaX,bbetaX))
  un<-runif(1)
  if (un<=ap) {
    betaX<<-proposal
    acceptX<<-acceptX+1
  } else {
    rejectX<<-rejectX+1
  }   
}
XUpdate1 <- function(i=0) {
  # Update X
  pOne<-matrix(0,tps,rgs)#p(X[i,j]=1|Y_1,...,Y_{t-1})
  p<-matrix(0,tps,rgs)   #p(X[i,j]=1|Y_1,...,Y_t)
  pOne[1,]<-rep(kX[1]/sum(kX),rgs)
  for (j in 1:tps) {
    if (j>1) {pOne[j,]<-kX[1]*(1-p[j-1,])+(1-kX[2])*p[j-1,]}
    pStarZero<-XLikelihood(j,0)*(1-pOne[j,])
    pStarOne <-XLikelihood(j,1)*pOne[j,]
    p[j,]<-pStarOne/(pStarOne+pStarZero)
  }
  for (j in tps:1) {
    if (j==tps) {
      ap<-p[j,]
    } else {
      # ap = p(X[t+1]|X[t]=1)*p(X[t]=1|Y1..t)/p(X[t+1]|Y1..t)
      ap<-(X[j+1,]+(-1)^X[j+1,]*kX[2])*p[j,]/(1-X[j+1,]+(-1)^(1-X[j+1,])*pOne[j+1,])
    }
    un<-runif(rgs)
    X[j,which(un<=ap)]<<-1
    X[j,which(un>ap)]<<-0
  }
  # Update kX
  tm<-XJumpFunction(X)
  proposal<-c(rbeta(1,aX[1]+tm[3],bX[1]+tm[1]),rbeta(1,aX[2]+tm[2],bX[2]+tm[4]))# (p01,p10)
  ap<-prod(dbinom(X[1,],1,proposal[1]/sum(proposal))/dbinom(X[1,],1,kX[1]/sum(kX)))
  un<-runif(1)
  if (un<=ap) {kX<<-proposal}
  # Update betaX
  proposal<-rnorm(1,betaX,sigmaX)
  ap<-betaXLikelihood(proposal)
  un<-runif(1)
  if (un<=ap & proposal>0) {# flat prior on (0,infty) symmetrical avoids double modes in posterior
    betaX<<-proposal
    acceptX<<-acceptX+1
  } else {
    rejectX<<-rejectX+1
  }   
}
XUpdate2 <- function(i=0) {
  # Update kX
  tm<-XJumpFunction(X)
  proposal<-c(rbeta(1,aX[1]+tm[3],bX[1]+tm[1]),rbeta(1,aX[2]+tm[2],bX[2]+tm[4]))# (p01,p10)
  ap<-prod(dbinom(X[1,],1,proposal[1]/sum(proposal))/dbinom(X[1,],1,kX[1]/sum(kX)))
  un<-runif(1)
  if (un<=ap) {kX<<-proposal}
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
  # Update X
  pOne<-matrix(0,tps,rgs)#p(X[i,j]=1|Y_1,...,Y_{t-1})
  p<-matrix(0,tps,rgs)   #p(X[i,j]=1|Y_1,...,Y_t)
  pOne[1,]<-rep(kX[1]/sum(kX),rgs)
  for (j in 1:tps) {
    if (j>1) {pOne[j,]<-kX[1]*(1-p[j-1,])+(1-kX[2])*p[j-1,]}
    pStarZero<-XLikelihood(j,0)*(1-pOne[j,])
    pStarOne <-XLikelihood(j,1)*pOne[j,]
    p[j,]<-pStarOne/(pStarOne+pStarZero)
  }
  for (j in tps:1) {
    if (j==tps) {
      ap<-p[j,]
    } else {
      # ap = p(X[t+1]|X[t]=1)*p(X[t]=1|Y1..t)/p(X[t+1]|Y1..t)
      ap<-(X[j+1,]+(-1)^X[j+1,]*kX[2])*p[j,]/(1-X[j+1,]+(-1)^(1-X[j+1,])*pOne[j+1,])
    }
    un<-runif(rgs)
    X[j,which(un<=ap)]<<-1
    X[j,which(un>ap)]<<-0
  }  
}
squashProd <- function(input) {
  output<-matrix(0,1,rgs)
  for (k in 1:rgs) {
    output[k]<-prod(input[wch[[k]]])    
  }
  return(output)
}
XRisk <- function() {X[,mbrg]*betaX}
XRisk2 <- function() {X[,mbrg]*rep(betaX[mbrg],each=tps)}
XLikelihoodRUX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+U+x*betaX)))}
XLikelihoodRUX2 <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+U+x*betaX[mbrg])))}
XLikelihoodRX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+x*betaX)))}
XLikelihoodUX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+U+x*betaX)))}
betaXLikelihoodRUX <- function(proposal) {prod(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*proposal))/dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*betaX)))}
betaXLikelihoodRUX2 <- function(proposal,j) {prod(dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*proposal))/dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*betaX[j])))}
betaXLikelihoodRX <- function(proposal) {prod(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*proposal))/dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*betaX)))}
XSample <- function(i) {
  if (i>burnin) {
    cumX<<-cumX+X
  }
  cat(c(t(kX),"\n"),file="kX.txt",append=TRUE)
  cat(t(XJumpFunction(X)),"\n",file="sumX.txt",append=TRUE)
  cat(betaX,"\n",file="betaX.txt",append=TRUE)
  if (Xmode>1) {
    cat((1-apply(1-X,2,prod))*betaX,"\n",file="betaXconditional.txt",append=TRUE)
  }
  cat(acceptX/(acceptX+rejectX),"\n",file="acceptanceX.txt",append=TRUE)
  acceptX<<-0
  rejectX<<-0
  if (XOutput) {
    cat(X,"\n",file="X.txt")
    cat(cumX,"\n",file="cumulativeX.txt")
  }
  if (XFullOutput) {
    cat(X,"\n",file="fullX.txt",append=TRUE)
  }  
}
XConvergence <- function() {
  plotPairs("kX",2)
  plotPairs("sumX",4,p=F)
  if (Xmode>1) {
    betaX<<-plotPairs("betaX",rgs,use=T)
  } else {
    betaX<<-plotPairs("betaX",use=T)
    # trick for returning betaX to a vector
    betaX<<-sum(betaX)
  }
  plotPairs("acceptanceX",length(acceptX),F,half=T)
  if (Xmode>1) {
    try(input<-scan("betaXconditional.txt"),T)
    input<-matrix(input,rgs,length(input)/rgs)
    rr<-matrix(0,rgs,1)
    for (k in 1:rgs) {
      w<-which(input[k,]>0)
      w<-w[which(w>=1+burnin/samplefreq)]
      if (length(w)==0) {
        betaX[k]<<-0
      } else {
        betaX[k]<<-mean(input[k,w])
        rr[k]<-mean(exp(input[k,w]))
      }
    }
    write.table(t(betaX),"posteriorbetaXconditional.txt",row.names=F,col.names=F)
    write.table(t(rr),"relativeriskbetaXconditional.txt",row.names=F,col.names=F)
  }
}
XAnalysis <- function() {
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
  write.table(t(X),"posteriorX.txt",row.names=F,col.names=F)
}
XAnalysis2 <- function() {
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
  if (Xmode<2) {
    betaX<<-Upload("betaX",1)
  } else {
    betaX<<-Upload("betaX",rgs)
  }
  kX<<-Upload("kX",2)
  X<<-scan("X.txt")
  X<<-matrix(X,tps,rgs)
  cumX<<-scan("cumulativeX.txt")
  cumX<<-matrix(cumX,tps,rgs)
}