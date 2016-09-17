# Defaults
sigmaW<-0.2
kW<-1
if (tidyup) {file.remove("W.txt")}
if (tidyup) {file.remove("kW.txt")}
if (tidyup) {file.remove("acceptanceW.txt")}
if (tidyup) {file.remove("sumW.txt")}
acceptW<-matrix(0,2,1)
rejectW<-matrix(0,2,1)
WInitialise <- function() {
  ws<<-max(wthr)
  cs<<-3*ws
  W<<-0.1*matrix(c(-1,0,1),1,cs)
  rl<<-list()
  mbl<<-list()
  tpl<<-list()
  cl<<-list()
  t0<-1:((2+tps)*mbs) # All days/mbs
  t1<-which(1+(t0-1)%%(tps+2)>2) # days 3 and after on the wthr scale
  t2<-which(1+(t0-1)%%(tps+2)>1)
  t2<-t2[which(1+(t2-1)%%(tps+2)<tps+2)] # days 2 to tps+1 on the wthr scale
  t3<-which(1+(t0-1)%%(tps+2)<tps+1) # days 1 to tps on the wthr scale
  for (i in 1:ws) {
    rl[[i]]<<-t1[which(wthr[t1]==i)]
    mbl[[i]]<<-1+(rl[[i]]-(1+(rl[[i]]-1)%%(tps+2)))/(tps+2) # mb list i
    tpl[[i]]<<-1+(rl[[i]]-1)%%(tps+2)-2 # time of case for list i converted to case scale
    cl[[i]]<<-tpl[[i]]+tps*(mbl[[i]]-1) # cases for list i converted to case scale
    rl[[ws+i]]<<-t2[which(wthr[t2]==i)]
    mbl[[ws+i]]<<-1+(rl[[ws+i]]-(1+(rl[[ws+i]]-1)%%(tps+2)))/(tps+2)
    tpl[[ws+i]]<<-1+(rl[[ws+i]]-1)%%(tps+2)-2+1
    cl[[ws+i]]<<-tpl[[ws+i]]+tps*(mbl[[ws+i]]-1)
    rl[[2*ws+i]]<<-t3[which(wthr[t3]==i)]
    mbl[[2*ws+i]]<<-1+(rl[[2*ws+i]]-(1+(rl[[2*ws+i]]-1)%%(tps+2)))/(tps+2)
    tpl[[2*ws+i]]<<-1+(rl[[2*ws+i]]-1)%%(tps+2)-2+2
    cl[[2*ws+i]]<<-tpl[[2*ws+i]]+tps*(mbl[[2*ws+i]]-1)
  }
}
WSetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aW<<-1
    bW<<-10^-2
  } else if (setpriors==2) {
    aW<<-1
    bW<<-0.25
  } else if (setpriors==3) {
    aW<<-10^-4
    bW<<-10^-4
  }    
}
WSumFunction<- function(W) {sum((W[1:(ws-1)]-W[2:ws])^2+(W[(ws+1):(2*ws-1)]-W[(ws+2):(2*ws)])^2+(W[(2*ws+1):(cs-1)]-W[(2*ws+2):cs])^2)}
WUpdate <- function(i=0) {
  # Gibb's step to update kW
  kW<<-rgamma(1,aW+3*(ws-1)/2,rate=(bW+WSumFunction(W)/2))
  for (j in 1:cs) {
    if (i%%2==0) {# Conditional prior proposal step to update W.
      if (j%%ws==1) {
        proposal<-rnorm(1,W[j+1],(1/kW)^(0.5))
      } else if (j%%ws==0) {
        proposal<-rnorm(1,W[j-1],(1/kW)^(0.5))
      } else {
        proposal<-rnorm(1,(W[j-1]+W[j+1])/2,(0.5/kW)^(0.5))
      }
      ap<-WLikelihood(j,proposal)
    } else {# Metropolis Hastings proposal step to update W.
      proposal<-rnorm(1,W[j],sd=sigmaW)
      if (j%%ws==1) {
        ap<-exp(-kW*((proposal-W[j+1])^2-(W[j]-W[j+1])^2)/2)*WLikelihood(j,proposal)
      } else if (j%%ws==0) {
        ap<-exp(-kW*((proposal-W[j-1])^2-(W[j]-W[j-1])^2)/2)*WLikelihood(j,proposal)
      } else {
        ap<-exp(-kW*((proposal-W[j-1])^2-(W[j]-W[j-1])^2+(proposal-W[j+1])^2-(W[j]-W[j+1])^2)/2)*WLikelihood(j,proposal)
      }
    }
    un<-runif(1)
    if (un<=ap) {
      W[j]<<-proposal
      acceptW[1+i%%2]<<-acceptW[1+i%%2]+1
    } else {
      rejectW[1+i%%2]<<-rejectW[1+i%%2]+1
    }
  }
}
WRisk <- function() {W[wthr[3:(tps+2),]]+W[ws+wthr[2:(tps+1),]]+W[2*ws+wthr[1:tps,]]}
WLikelihoodW <- function(j,proposal) {
  # calculate the likelihood ratio
  if (j<=ws) {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+proposal+W[ws+wthr[rl[[j]]-1]]+W[2*ws+wthr[rl[[j]]-2]]))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+W[j]+W[ws+wthr[rl[[j]]-1]]+W[2*ws+wthr[rl[[j]]-2]]))))
  } else if (j<=2*ws) {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+W[wthr[rl[[j]]+1]]+proposal+W[2*ws+wthr[rl[[j]]-1]]))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+W[wthr[rl[[j]]+1]]+W[j]+W[2*ws+wthr[rl[[j]]-1]]))))
  } else {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+W[wthr[rl[[j]]+2]]+W[ws+wthr[rl[[j]]+1]]+proposal))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+W[wthr[rl[[j]]+2]]+W[ws+wthr[rl[[j]]+1]]+W[j]))))
  }
}
WLikelihoodRUW <- function(j,proposal) {
  # calculate the likelihood ratio
  if (j<=ws) {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+proposal+W[ws+wthr[rl[[j]]-1]]+W[2*ws+wthr[rl[[j]]-2]]))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+W[j]+W[ws+wthr[rl[[j]]-1]]+W[2*ws+wthr[rl[[j]]-2]]))))
  } else if (j<=2*ws) {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+W[wthr[rl[[j]]+1]]+proposal+W[2*ws+wthr[rl[[j]]-1]]))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+W[wthr[rl[[j]]+1]]+W[j]+W[2*ws+wthr[rl[[j]]-1]]))))
  } else {
    return(prod(dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+W[wthr[rl[[j]]+2]]+W[ws+wthr[rl[[j]]+1]]+proposal))/dpois(cases[cl[[j]]],n[mbl[[j]]]*exp(fe+R[tpl[[j]]]+U[mbl[[j]]]+W[wthr[rl[[j]]+2]]+W[ws+wthr[rl[[j]]+1]]+W[j]))))
  }
}
WSample <- function() {
  fe<<-fe+mean(W[1:ws])
  W[1:ws]<<-W[1:ws]-mean(W[1:ws])
  fe<<-fe+mean(W[(ws+1):(2*ws)])
  W[(ws+1):(2*ws)]<<-W[(ws+1):(2*ws)]-mean(W[(ws+1):(2*ws)])
  fe<<-fe+mean(W[(2*ws+1):cs])
  W[(2*ws+1):cs]<<-W[(2*ws+1):cs]-mean(W[(2*ws+1):cs])
  cat(c(t(W),"\n"),file="W.txt",append=TRUE,sep=" ")
  cat(kW,file="kW.txt",append=TRUE,sep="\n")
  cat(c(t(acceptW[]/(acceptW[]+rejectW[])),"\n"),file="acceptanceW.txt",append=TRUE,sep=" ")
  acceptW[]<<-rep(0,2)
  rejectW[]<<-rep(0,2)
  cat(WSumFunction(W),"\n",file="sumW.txt",append=TRUE)
}
WConvergence <- function() {
  plotPairs("kW")
  plotPairs("sumW",p=F)
  plotPairs("acceptanceW",length(acceptW),F,ha=T)
}
WTraces <- function() {
  W<<-plotPairs("W",cs,z=T)
}
WAnalysis <- function() {
  m<-max(max(W),-min(W))
  plot(1:cs,t="n",xlab="category",ylab="W",main="Weather Analysis",xlim=c(1,cs),ylim=c(-m,m))
  lines(1:ws,W[1:ws])
  lines((ws+1):(2*ws),W[(ws+1):(2*ws)])
  lines((2*ws+1):cs,W[(2*ws+1):cs])
  rr<-scan('relativeriskW.txt')
  plot(1:cs,t="n",xlab="category",ylab="Relative Risk",main="Weather Analysis",xlim=c(1,cs),ylim=c(0,max(rr)))
  lines(1:ws,rr[1:ws])
  lines((ws+1):(2*ws),rr[(ws+1):(2*ws)])
  lines((2*ws+1):cs,rr[(2*ws+1):cs])
}