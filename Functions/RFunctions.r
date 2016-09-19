#Defaults
sigmaR<-1
Rblock<-c(4,5,9,11)

# load in likelihood
Rcpp::sourceCpp('src/likelihood.cpp')
Rcpp::sourceCpp('src/rmvnorm.cpp')
Rcpp::sourceCpp('src/update_r.cpp')
Rcpp::sourceCpp('src/rbernoulli.cpp')

# This returns the sum of R squared, i.e. the variance of R's for the kR update
RSumFunction <- function(R) {
  lenR <- length(R)
  sum((R[1:(lenR-2)] - 2*R[2:(lenR-1)] + R[3:lenR])^2)
}

# computes the scaled eigen-vectors of A
eigen_scale <- function(A) {
  eig <- eigen(A, symmetric=TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values,0)))
}

RInitialize <- function() {
  Rsigma<<-list()
  Rsigma_eigen<<-list()
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
    Rsigma_eigen[[i]]<<-eigen_scale(Rsigma[[i]])
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

  if (params$tidyup) {file.remove(file.path(params$outpath, "R.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "kR.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "acceptanceR.txt"))}
  if (params$tidyup) {file.remove(file.path(params$outpath, "sumR.txt"))}

  state <- list(R  = matrix(rnorm(params$tps,0,1),1,params$tps),
                kR = 1,
                fe = -10,
                acceptR = matrix(0,1+length(Rblock),1),
                rejectR = matrix(0,1+length(Rblock),1))
  return(state)
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

RUpdate <- function(i=0, state) {
  # save current state
  R  <- state$R
  acceptR <- state$acceptR
  rejectR <- state$rejectR

  state <- update_r(cases, n, i, state, list(aR=aR, bR=bR))
  kR <- state$kR

  lenR <- length(R)
  method<-1+i%%(1+length(Rblock))
  endmethod<-rbernoulli(0.5)
  j<-1 #start of update block
  while (j<=lenR) {
    k<-j# end of update block
    if (method>length(Rblock) || (endmethod==0 && (j<3 || j>lenR-2))) {
      # Metropolis Hastings proposal step to update R.
      state$R <- R
      out <- update_r_mh(cases, n, mbrg, state, list(sigmaR=sigmaR), j-1)
      proposal <- out$proposal
      ap <- out$ap
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
        mu <- Rmu[[method]][,1]*R[j-2]+Rmu[[method]][,2]*R[j-1]+Rmu[[method]][,3]*R[k+1]+Rmu[[method]][,4]*R[k+2];
        proposal <- rmvnorm(mu, Rsigma_eigen[[method]]/sqrt(kR))
      }
      ap<-RLikelihood(j,k,R[j:k],proposal,state)
    }
    un<-runif(1)
#    cat("ap=", ap, "un=", un, "\n")
    if (ap >= 0 | un<=exp(ap)) {
      R[j:k]<-proposal
      acceptR[method]<-acceptR[method]+1
    } else {
      rejectR[method]<-rejectR[method]+1
    }
    j<-k+1
  }
  state$kR <- kR
  state$fe <- state$fe + mean(R)
  state$R  <- R - mean(R)

#  cat("Ending RUpdate, dim(R)=", dim(state$R), "\n")
  state$acceptR <- acceptR
  state$rejectR <- rejectR
  return(state)
}

RRisk <- function(state) {
  rep(state$R, params$mbs)
}

RLikelihoodR <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)), log=TRUE))
}

RLikelihoodRS <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop+S[j:k],mbs)), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr+S[j:k],mbs)), log=TRUE))
}

RLikelihoodRU <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(U,each=k-j+1)), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)), log=TRUE))
}

RLikelihoodRUW <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(U,each=k-j+1)+W[wthr[j:k+2,]]+W[ws+wthr[j:k+1,]]+W[2*ws+wthr[j:k,]]), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+W[wthr[j:k+2,]]+W[ws+wthr[j:k+1,]]+W[2*ws+wthr[j:k,]]), log=TRUE))
}

RLikelihoodRX <- function(j,k,curr,propl) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+X[j:k,mbrg]*betaX), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+X[j:k,mbrg]*betaX), log=TRUE))
}

RLikelihoodRUX <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX), log=TRUE))
}

RLikelihoodRUX2 <- function(j,k,curr,prop,state) {
  # lambda here is the current rate
  #lambda_curr <- rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1))
  #lambda_prop <- lambda_curr * rep(exp(prop-curr),mbs)
  # sum(dpois(cases[j:k,],lambda_prop, log=TRUE)-
  #     dpois(cases[j:k,],lambda_curr, log=TRUE))
  r_likelihood_rux2(cases, n, state$fe, state$U, state$X, mbrg, state$betaX,
                    curr, prop, j, k)
}

RLikelihoodRUX3 <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  if (j==1) { # only called if k==j
    sum(dpois(cases[1,],n*exp(fe+rep(prop,mbs)+U+X[1,mbrg]*betaX[mbrg]), log=TRUE)-
        dpois(cases[1,],n*exp(fe+rep(curr,mbs)+U+X[1,mbrg]*betaX[mbrg]), log=TRUE))
  } else {
    sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(U,each=k-j+1)+(X[(j-1):(k-1),mbrg]+X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)), log=TRUE)-
        dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+(X[(j-1):(k-1),mbrg]+X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)), log=TRUE))
  }
}

RLikelihoodRUX4 <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  if (j==1) { # only called if k==j
    sum(dpois(cases[1,],n*exp(fe+rep(prop,mbs)+U+X[1,mbrg]*betaX[mbrg]), log=TRUE)-
        dpois(cases[1,],n*exp(fe+rep(curr,mbs)+U+X[1,mbrg]*betaX[mbrg]), log=TRUE))
  } else {
    sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(U,each=k-j+1)+pmax(X[(j-1):(k-1),mbrg],X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)), log=TRUE)-
        dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(U,each=k-j+1)+pmax(X[(j-1):(k-1),mbrg],X[j:k,mbrg])*rep(betaX[mbrg],each=k-j+1)), log=TRUE))
  }
}

RLikelihoodBR <- function(j,k,curr,prop) {
  mbs <- ncol(n)
  sum(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(prop,mbs)+rep(B[Blink],each=k-j+1)), log=TRUE)-
      dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(curr,mbs)+rep(B[Blink],each=k-j+1)), log=TRUE))
}

RSample <- function(state) {
  cat(c(t(state$R),"\n"),file=file.path(params$outpath, "R.txt"),append=TRUE,sep=" ")
  cat(state$kR,file=file.path(params$outpath, "kR.txt"),append=TRUE,sep="\n")
  cat(c(t(state$acceptR[]/(state$acceptR[]+state$rejectR[])),"\n"),file=file.path(params$outpath, "acceptanceR.txt"),append=TRUE,sep=" ")
  cat(RSumFunction(state$R),"\n",file=file.path(params$outpath, "sumR.txt"),append=TRUE)
}

RInitAcceptance <- function(state) {
  state$acceptR <- matrix(0,1+length(Rblock),1)
  state$rejectR <- matrix(0,1+length(Rblock),1)
  return(state)
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
