#Defaults
library(MASS)
sigmaR<-1
kR<-1
Rblock<-c(4,5,9,11)
R<-matrix(rnorm(params$tps,0,1),1,params$tps)
fe <- -10

# This returns the sum of R squared, i.e. the variance of R's for the kR update
RSumFunction <- function(R) {
  lenR <- length(R)
  sum((R[1:(lenR-2)] - 2*R[2:(lenR-1)] + R[3:lenR])^2)
}

RInitialize <- function() {
  Rsigma<<-list()
  Rmu<<-list()
  for (i in 1:length(Rblock)) {
    K<-6*diag(Rblock[i])#Structure matrix for R
    for (j in 1:Rblock[i]) {
      if (j>2) {K[j,j-2]<-1}
      if (j>1) {K[j,j-1]<--4}
      if (j<Rblock[i]) {K[j,j+1]<--4}
      if (j<Rblock[i]-1) {K[j,j+2]<-1}
    }
    Rsigma[[i]]<<-solve(K)
    t1<-matrix(0,Rblock[i],1)
    t2<-matrix(0,Rblock[i],1)
    t3<-matrix(0,Rblock[i],1)
    t4<-matrix(0,Rblock[i],1)
    t1[1]<-1
    t2[1:2]<-c(-4,1)
    t3[(Rblock[i]-1):Rblock[i]]<-c(1,-4)
    t4[Rblock[i]]<-1
    Rmu[[i]]<<-cbind(-Rsigma[[i]]%*%t1,-Rsigma[[i]]%*%t2,-Rsigma[[i]]%*%t3,-Rsigma[[i]]%*%t4)
  }
  acceptR<<-matrix(0,1+length(Rblock),1)
  rejectR<<-matrix(0,1+length(Rblock),1)
  if (params$tidyup) {file.remove(file.path(params$outpath, "R.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "kR.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "acceptanceR.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "sumR.txt"))}
}

RSetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aR<<-1
    bR<<-10^-4
  } else if (setpriors==2) {
    aR<<-5
    bR<<-5/500
  } else if (setpriors==3) {
    aR<<-10^-4
    bR<<-10^-4
  }
}

RUpdate <- function(i=0) {
  lenR <- length(R)
  # Gibb's step to update kR
  kR<<-rgamma(1,aR+(lenR-2)/2,rate=bR+RSumFunction(R)/2)
  method<-1+i%%(1+length(Rblock))
  endmethod<-rbinom(1,1,0.5)
  j<-1 #start of update block
  while (j<=lenR) {
    k<-j# end of update block
    if (method>length(Rblock) || (endmethod==0 && (j<3 || j>lenR-2))) {
      # Metropolis Hastings proposal step to update R.
      proposal<-rnorm(1,R[j],sd=sigmaR)
      ap<-RLikelihood(j,k,proposal)
      # full conditional component of ap
      if (j>2) {ap<-ap*exp(-kR*((R[j-2]-2*R[j-1]+proposal)^2-(R[j-2]-2*R[j-1]+R[j])^2)/2)}
      if (j>1 && j<lenR) {ap<-ap*exp(-kR*((R[j-1]-2*proposal+R[j+1])^2-(R[j-1]-2*R[j]+R[j+1])^2)/2)}
      if (j<lenR-1) {ap<-ap*exp(-kR*((proposal-2*R[j+1]+R[j+2])^2-(R[j]-2*R[j+1]+R[j+2])^2)/2)}     
    } else {
      # Conditional Prior Proposal step to update R
      if (j==1) {
        proposal<-rnorm(1,2*R[2]-R[3],(1/kR)^(0.5))
      } else if (j==2) {
        proposal<-rnorm(1,0.4*R[1]+0.8*R[3]-0.2*R[4],(0.2/kR)^(0.5))
      } else if (j==lenR-1) {
        proposal<-rnorm(1,-0.2*R[lenR-3]+0.8*R[lenR-2]+0.2*R[lenR],(0.2/kR)^(0.5))
      } else if (j==lenR) {
        proposal<-rnorm(1,-R[lenR-2]+2*R[lenR-1],(1/kR)^(0.5))
      } else {
        #proposal<-rnorm(1,-R[j-2]/6+R[j-1]*4/6+R[j+1]*4/6-R[j+2]/6,(1/6/kR)^(0.5))
        k<-j+Rblock[method]-1
        if (k>=lenR-1) {   # TODO: This breaks in a whole bunch of cases if Rblock[method] is large enough
          j<-j-(k-lenR+2)
          k<-lenR-2
        }
        proposal<-mvrnorm(1,Rmu[[method]][,1]*R[j-2]+Rmu[[method]][,2]*R[j-1]+Rmu[[method]][,3]*R[k+1]+Rmu[[method]][,4]*R[k+2],Rsigma[[method]]/kR)
      }
      ap<-RLikelihood(j,k,proposal)
    }
    un<-runif(1)
#    cat("ap=", ap, "un=", un, "\n")
    if (un<=ap) {
      R[j:k]<<-proposal
      acceptR[method]<<-acceptR[method]+1
    } else {
      rejectR[method]<<-rejectR[method]+1
    }
    j<-k+1
  }
}

RRisk <- function() {
  rep(R, params$mbs)
}

RLikelihoodR <- function(j,k,proposal) {
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(R[j:k],mbs))))
}

RLikelihoodRS <- function(j,k,proposal) {
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal+S[j:k],mbs)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(R[j:k]+S[j:k],mbs))))
}

RLikelihoodRU <- function(j,k,proposal) {# calculate the likelihood ratio
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(R[j:k],mbs)+rep(U,each=k-j+1))))
}

RLikelihoodRUW <- function(j,k,proposal) {# calculate the likelihood ratio
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+W[wthr[j:k+2,]]+W[ws+wthr[j:k+1,]]+W[2*ws+wthr[j:k,]]))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+W[wthr[j:k+2,]]+W[ws+wthr[j:k+1,]]+W[2*ws+wthr[j:k,]])))
}

RLikelihoodRX <- function(j,k,proposal) {
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+X[j:k,mbrg]*betaX))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+X[j:k,mbrg]*betaX)))
}

RLikelihoodRUX <- function(j,k,proposal) {
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX)))
}

RLikelihoodRUX2 <- function(j,k,proposal) {
  mbs <- ncol(n)
  num <- dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1)))
  den <- dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1)))
  prod(num / den)
}

RLikelihoodRUX3 <- function(j,k,proposal) {
  mbs <- params$mbs
  if (j==1) {
    prod(dpois(cases[1,],n*exp(fe+rep(proposal,mbs)+U+X[1,mbrg]*betaX[mbrg]))/dpois(cases[1,],n*exp(fe+rep(R[1],mbs)+U+X[1,mbrg]*betaX[mbrg])))
  } else {
    prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+(X[(j-1):(k-1),mbrg]+X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+(X[(j-1):(k-1),mbrg]+X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1))))
  }
}

RLikelihoodRUX4 <- function(j,k,proposal) {
  mbs <- params$mbs
  if (j==1) {
    prod(dpois(cases[1,],n*exp(fe+rep(proposal,mbs)+U+X[1,mbrg]*betaX[mbrg]))/dpois(cases[1,],n*exp(fe+rep(R[1],mbs)+U+X[1,mbrg]*betaX[mbrg])))
  } else {
    prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+pmax(X[(j-1):(k-1),mbrg],X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+pmax(X[(j-1):(k-1),mbrg],X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1))))
  }
}

RLikelihoodBR <- function(j,k,proposal) {
  mbs <- params$mbs
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(B[Blink],each=k-j+1)))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(R[j:k],mbs)+rep(B[Blink],each=k-j+1))))
}

RSample <- function() {
  fe<<-fe+mean(R)
  R<<-R-mean(R)
  cat(c(t(R),"\n"),file=file.path(params$outpath, "R.txt"),append=TRUE,sep=" ")
  cat(kR,file=file.path(params$outpath, "kR.txt"),append=TRUE,sep="\n")
  cat(c(t(acceptR[]/(acceptR[]+rejectR[])),"\n"),file=file.path(params$outpath, "acceptanceR.txt"),append=TRUE,sep=" ")
  acceptR[]<<-rep(0,1+length(Rblock))
  rejectR[]<<-rep(0,1+length(Rblock))
  cat(RSumFunction(R),"\n",file=file.path(params$outpath, "sumR.txt"),append=TRUE)
}

RConvergence <- function() {
  plotPairs("kR")
  plotPairs("sumR",p=F)
  plotPairs("acceptanceR",length(Rblock)+1,F,ha=T)
}

RTraces <- function() {
  R<<-plotPairs("R",params$tps,z=T)
}

RAnalysis <- function() {
  tps <- params$tps

  input<-scan(file.path(params$outpath, "expectedCases.txt"))
  st<-1+burnin/samplefreq
  en<-length(input)/tps
  input<-matrix(input,tps,en)
  out1<-apply(cases,1,sum)
  out2<-apply(input[,st:en],1,mean)
  if (incX | incS) {
    input3<-scan(file.path(params$outpath, "smoothedCases.txt"))
    input3<-matrix(input3,tps,en)
    out3<-apply(input3[,st:en],1,mean)
  }
  plot(1:tps,t="n",xlab="week",ylab="",main="Full Temporal Analysis",xlim=c(1,tps),ylim=c(0,max(out1)))
  lines(1:tps,out1,col=3)
  lines(1:tps,out2,col=2)
  if (params$incX | params$incS) {lines(1:tps,out3,col=4)}
  if (tps==312) {
    for (i in 0:5) {
      plot(52*i+1:52,t="n",xlab="week",ylab="",main="Yearly Temporal Analysis",xlim=c(52*i+1,52*i+52),ylim=c(0,max(out1)))
      lines(52*i+1:52,out1[52*i+1:52],col=3)
      lines(52*i+1:52,out2[52*i+1:52],col=2)
      if (params$incX | params$incS) {lines(52*i+1:52,out3[52*i+1:52],col=4)}
    }
  }
}
