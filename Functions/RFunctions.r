#Defaults
sigmaR<-1
Rblock<-c(4,5,9,11)

library(epiclustR)

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

R_covariance_matrix <- function(n) {
  # Generates the structure matrix for a random walk of order 2
  K<-6*diag(2 + n + 2) # need 2 on either side for ends
  for (j in 1:nrow(K)) {
    if (j>2) {K[j,j-2]<-1}
    if (j>1) {K[j,j-1]<--4}
    if (j+1<=ncol(K)) {K[j,j+1]<--4}
    if (j+2<=ncol(K)) {K[j,j+2]<-1}
  }
  K[1,1:2] <- K[n+4,n+4:3] <- c(1, -2)
  K[2,1:2] <- K[n+3,n+4:3] <- c(-2, 5)
  K
}

RInitialize <- function() {
  Rsigma<<-list()
  Rsigma_eigen<<-list()
  Rbefore<<-list()
  Rafter<<-list()
  for (i in 1:length(Rblock)) {
    o <- 2 # offset
    K <- R_covariance_matrix(Rblock[i])
    invK <- solve(K[1:Rblock[i]+o, 1:Rblock[i]+o])
    Rsigma[[i]]<<- invK
    Rsigma_eigen[[i]]<<-eigen_scale(Rsigma[[i]])
    Kbef <- K[1:Rblock[i]+o,seq_len(min(o,2))+max(o,2)-2, drop=FALSE]
    Kaft <- K[1:Rblock[i]+o,o + Rblock[i] + seq_len(min(4-o,2)), drop=FALSE]
    Rbefore[[i]] <<- -invK %*% Kbef
    Rafter[[i]] <<-  -invK %*% Kaft
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

Update <- function(state, control, burnin = FALSE) {
  # transfer the priors
  prior <- list(aR=aR, bR=bR,
                aU=aU, bU=bU,
                aX=aX, bX=bX, abetaX=abetaX, bbetaX=bbetaX)

  # transfer the control of proposals
  proposal <- list(sigmaR=sigmaR, Rbefore=Rbefore, Rafter=Rafter, Rsigma=Rsigma_eigen,
                   sigmaU=sigmaU,
                   sigmaX=sigmaX)

  # data
  data <- list(cases=cases, popn=n, mbrg=mbrg, nb=weight, rgmb=wch)

  # do the updates
  state <- update(data, state, prior, c(control, proposal))

  # TODO: Ideally we'd remove this. Problem is full X is large to save
  #       in terms of the posterior, and if we want it right we have
  #       to get rid of the burnin period.
  #       I guess saving X to a binary file might be way more efficient
  #       though? About 16 times smaller.
  if (!burnin) {
    state$cumX<-state$cumX+state$X
  }
  return(state)
}

RUpdate <- function(i=0, state) {
  # call into C++ land
  update_r(cases, n, mbrg, i, state, list(aR=aR, bR=bR, sigmaR=sigmaR), Rmu, Rsigma_eigen)
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
