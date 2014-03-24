sigmaS<-1
kS<-1
S<-matrix(rnorm(tps,0,1),1,tps)
Smu<-matrix(0,tps,tps)
SKjj<-matrix(0,1,tps)
for (j in 1:tps) {
  for (k in 0:51) {
    if (j-k>0 & j-k+51<=tps) {
      SKjj[j]<-SKjj[j]+1
      Smu[j,(j-k):(j-k+51)]<-Smu[j,(j-k):(j-k+51)]+1
    }
  }
  Smu[j,j]<-0
  Smu[j,]<--Smu[j,]/SKjj[j]
}
acceptS<-rep(0,2)
rejectS<-rep(0,2)
if (tidyup) {
  file.remove("S.txt")
  file.remove("kS.txt")
  file.remove("acceptanceS.txt")
  file.remove("sumS.txt")
  file.remove("smoothedCases.txt")
}
SSetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aS<<-1
    bS<<-10^-4
  } else if (setpriors==2) {
    aS<<-1
    bS<<-0.25
  } else if (setpriors==3) {
    aS<<-10^-4
    bS<<-10^-4
  }
}
SSumFunction <- function() {
  output<-0
  for (j in 1:(tps-51)) {
    output<-output+sum(S[j:(j+51)])^2
  }
}
SUpdate <- function(i=0) {
  # Gibb's step to update kS
  kS<<-rgamma(1,aS+(tps-51)/2,rate=(bS+SSumFunction()/2))
  for (j in 1:tps) {
    if (i%%2==0) {
      proposal<-rnorm(1,Smu[j,]*S,(1/kS/SKjj[j])^(0.5))
      ap<-SLikelihood(j,proposal)
    } else {
      proposal<-rnorm(1,S[j],sd=sigmaS)
      temp<-0
      for (k in 0:51) {
        if (j-k>0 & j-k+51<=tps) {temp<-temp+(sum(S[(j-k):(j-k+51)])-S[j]+proposal)^2-sum(S[(j-k):(j-k+51)])^2}
      }
      ap<-exp(-kS*temp/2)*SLikelihood(j,proposal)
    }
    un<-runif(1)
    if (un<=ap) {
      S[j]<<-proposal
      acceptS[1+i%%2]<<-acceptS[1+i%%2]+1
    } else {
      rejectS[1+i%%2]<<-rejectS[1+i%%2]+1
    }
  }
}
SSample <- function() {
  #fe<<-fe+mean(S)
  #S<<-S-mean(S)
  cat(c(t(S),"\n"),file="S.txt",append=TRUE,sep=" ")
  cat(kS,file="kS.txt",append=TRUE,sep="\n")
  cat(c(t(acceptS[]/(acceptS[]+rejectS[])),"\n"),file="acceptanceS.txt",append=TRUE,sep=" ")
  acceptS<<-rep(0,2)
  rejectS<<-rep(0,2)
  cat(SSumFunction(),"\n",file="sumS.txt",append=TRUE)
}
SLikelihoodRS <- function(j,proposal) {prod(dpois(cases[j,],n*exp(fe+rep(R[j]+proposal,mbs)))/dpois(cases[j,],n*exp(fe+rep(R[j]+S[j],mbs))))}
SRisk <- function() {rep(S,mbs)}
SConvergence <- function() {
  plotPairs("kS")
  plotPairs("sumS",p=F)
  plotPairs("acceptanceS",length(acceptS),F,ha=T)
}
STraces <- function() {
  S<<-plotPairs("S",tps,z=T)
}
SAnalysis <- function() {
  m<-max(-min(S),max(S))
  cols<-rainbow(6,start=2/6,end=5/6)
  ax<-1:52
  plot(ax,t="n",xlab="week",ylab="",main="Periodic Analysis",xlim=c(1,52),ylim=c(-m,m))
  for (i in 0:5) {
    lines(ax,S[ax+52*i],col=cols[1+i])
  }
  legend(0,m,2001:2006,lwd=1,col=cols)
}