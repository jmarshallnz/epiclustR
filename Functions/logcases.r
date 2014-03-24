# logcase = 0: all ones - this should never be needed
# logcase = 1: one only if the previous week contains a case
# logcase = 2: log(Y_{t-1,j})
if (logcase==0) {
  logcases<-matrix(1,tps,mbs)
} else {
  logcases<-matrix(0,tps,mbs)#This is mbs - see likelihoods.
  for (j in 1:rgs) {
    cs<-apply(cases[1:(tps-1),wch[[j]]],1,sum)
    if (logcase==1) {logcases[1+which(cs>0),wch[[j]]]<-1}
    if (logcase==2) {logcases[2:tps,wch[[j]]]<-rep(log(1+cs),lch[j])}
  }
}
rm(cs)
RLikelihoodRX <- function(j,k,proposal) {
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+X[j:k,mbrg]*betaX*logcases[j:k,]))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+X[j:k,mbrg]*betaX*logcases[j:k,])))
}
RLikelihoodRUX <- function(j,k,proposal) {
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX*logcases[j:k,]))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*betaX*logcases[j:k,])))
}
RLikelihoodRUX2 <- function(j,k,proposal) {
  prod(dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep(proposal,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1)*logcases[j:k,]))/dpois(cases[j:k,],rep(n,each=k-j+1)*exp(fe+rep( R[j:k] ,mbs)+rep(U,each=k-j+1)+X[j:k,mbrg]*rep(betaX[mbrg],each=k-j+1)*logcases[j:k,])))
}
ULikelihoodRUX <- function(j,proposal) {prod(dpois(cases[,j],n[j]*exp(fe+R+proposal+betaX*X[,mbrg[j]]*logcases[,j]))/dpois(cases[,j],n[j]*exp(fe+R+U[j]+betaX*X[,mbrg[j]]*logcases[,j])))}
XRisk <- function() {X[,mbrg]*betaX*logcases}
XRisk2 <- function() {X[,mbrg]*rep(betaX[mbrg],each=tps)*logcases}
XLikelihoodRUX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+U+x*betaX*logcases[j,])))}
XLikelihoodRUX2 <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+U+x*betaX[mbrg]*logcases[j,])))}
XLikelihoodRX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+R[j]+x*betaX*logcases[j,])))}
XLikelihoodUX <- function(j,x) {squashProd(dpois(cases[j,],n*exp(fe+U+x*betaX*logcases[j,])))}
betaXLikelihoodRUX <- function(proposal) {prod(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*proposal*logcases))/dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+rep(U,each=tps)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*betaX*logcases)))}
betaXLikelihoodRUX2 <- function(proposal,j) {prod(dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*proposal*logcases))/dpois(cases[,wch[[j]]],rep(n[wch[[j]]],each=tps)*exp(fe+rep(R,lwch[j])+rep(U[wch[[j]]],each=tps)+X[rep(j-1,each=tps)*tps+rep(1:tps,lwch[j])]*betaX[j]*logcases)))}
betaXLikelihoodRX <- function(proposal) {prod(dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*proposal*logcases))/dpois(cases,rep(n,each=tps)*exp(fe+rep(R,mbs)+X[rep(mbrg-1,each=tps)*tps+rep(1:tps,mbs)]*betaX*logcases)))}
#
# Defaults
#
RLikelihood <- RLikelihoodRUX
ULikelihood <- ULikelihoodRUX
XLikelihood <- XLikelihoodRUX
betaXLikelihood <- betaXLikelihoodRUX