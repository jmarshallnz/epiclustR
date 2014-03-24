# Defaults
urs<-7
sds<-10
aB<-0.001
bB<-0.001
sigmaB<-0.5
kB<-1
B<-matrix(rnorm(urs*sds,0,1),urs,sds)
Bmax<-urs*sds
Bklink<-rep(1,Bmax)
Bks<-1
Blink<-ur+(sdi-1)*urs
Bn<-list()
acceptB<-matrix(0,2,1)
rejectB<-matrix(0,2,1)
if (tidyup) {
  file.remove("B.txt")
  file.remove("kB.txt")
  file.remove("acceptanceB.txt")
  file.remove("sumB.txt")
}
BInitialise <- function() {
  if (Bmode==0) {# Rectangular array ur x sdi
    for (j in 1:Bmax) {
      Bn[[j]]<<-c(j-urs,j+urs)
      if ((j-1)%%urs>0) {Bn[[j]]<<-c(Bn[[j]],j-1)}
      if ((j-1)%%urs<urs-1) {Bn[[j]]<<-c(Bn[[j]],j+1)}
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=Bmax]
    }
    BSumFunction <<- function() {sum((B[2:urs,]-B[1:(urs-1),])^2)+sum((B[,2:sds]-B[,1:(sds-1)])^2)}
    BRank<<-urs*sds-1
  } else if (Bmode==1) {#7 lines of 10 ur x sdi
    for (j in 1:Bmax) {
      Bn[[j]]<<-c(j-urs,j+urs)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=Bmax]
    }
    BSumFunction <<- function() {sum((B[,2:sds]-B[,1:(sds-1)])^2)}
    BRank<<-urs*(sds-1)
  } else if (Bmode==2) {#1 line for chicken farm proximity
    sdi<<-scan(paste(datapath,region,"\\Farms\\ChickenFarmDistance.txt",sep=""))
    sdi<<-matrix(1-min(sdi)+sdi,1,mbs)
    sds<<-max(sdi)
    Bmax<<-sds
    Bklink<<-rep(1,Bmax)
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:Bmax) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=Bmax]
    }
    BSumFunction <<- function() {sum((B[2:sds]-B[1:(sds-1)])^2)}
    BRank<<-sds-1
    Blink<<-sdi
  } else if (Bmode>2 & Bmode<6) {# 1 line for animal density
    if (Bmode==3) {input<-scan(paste(datapath,region,"\\Farms\\ChickenDensity.txt",sep=""))}
    if (Bmode==4) {input<-scan(paste(datapath,region,"\\Farms\\SheepDensity.txt",sep=""))}
    if (Bmode==5) {input<-scan(paste(datapath,region,"\\Farms\\CowDensity.txt",sep=""))}
    if (Bmode==5.1) {input<-scan(paste(datapath,region,"\\Farms\\DairyDensity.txt",sep=""))}
    if (Bmode==5.2) {input<-scan(paste(datapath,region,"\\Farms\\BeefDensity.txt",sep=""))}
    posput<-input[input>0]
    thresholds<-c(0,quantile(posput,c(0.25,0.5,0.75)))
    cat(thresholds,"\n")
    sdi<<-matrix(0,1,mbs)
    for (i in 1:mbs) {
      sdi[i]<<-1+length(which(thresholds<input[i]))
    }
    sdi<<-1-min(sdi)+sdi
    sds<<-max(sdi)
    Bmax<<-sds
    Bklink<<-rep(1,Bmax)
    for (i in 1:sds) {
      cat("Group",i,"has",length(which(sdi==i)),"members.\n")
    }
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:Bmax) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=Bmax]
    }
    BSumFunction <<- function() {sum((B[2:sds]-B[1:(sds-1)])^2)}
    BRank<<-sds-1
    Blink<<-sdi    
  } else if (Bmode==6) {#1 line for chicken farm proximity (rural) and 1 line for sdi (urban) 
    input<-scan(paste(datapath,region,"\\Farms\\ChickenFarmDistance.txt",sep=""))
    input<-matrix(1-min(input)+input,1,mbs)
    sdi[which(ur<7)]<<-10+input[which(ur<7)]
    sds<<-max(sdi)
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:10) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=10]
    }
    for (j in 11:sds) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=11 & Bn[[j]]<=sds]
    }
    BSumFunction <<- function() {sum((B[2:10]-B[1:9])^2)+sum((B[12:sds]-B[11:(sds-1)])^2)}
    BRank<<-sds-2
    Bmax<<-sds
    Bklink<<-rep(1,Bmax)
    Blink<<-sdi
  } else if (Bmode>6 & Bmode<12) {#1 line for animal density (rural), 1 for sdi (urban)
    if (Bmode==7) {input<-scan(paste(datapath,region,"\\Farms\\ChickenDensity.txt",sep=""))}
    if (Bmode==8) {input<-scan(paste(datapath,region,"\\Farms\\SheepDensity.txt",sep=""))}
    if (Bmode==9) {input<-scan(paste(datapath,region,"\\Farms\\CowDensity.txt",sep=""))}
    if (Bmode==10) {input<-scan(paste(datapath,region,"\\Farms\\BeefDensity.txt",sep=""))}
    if (Bmode==11) {input<-scan(paste(datapath,region,"\\Farms\\DairyDensity.txt",sep=""))}
    for (i in 1:mbs) {
      input[i]<-1+length(which(thresholds<input[i]))
    }
    sdi[which(ur<7)]<<-10+input[which(ur<7)]
    sds<<-max(sdi)
    Bmax<<-sds
    Bklink<<-rep(1,Bmax)
    for (i in 1:sds) {
      cat("Group",i,"has",length(which(sdi==i)),"members.\n")
    }
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:10) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=10]
    }
    for (j in 11:sds) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=11 & Bn[[j]]<=sds]
    }
    BSumFunction <<- function() {sum((B[2:10]-B[1:9])^2)+sum((B[12:sds]-B[11:(sds-1)])^2)}
    BRank<<-sds-2
    Blink<<-sdi
  } else if (Bmode==12) {#1 line for chicken farm proximity (rural) and 1 line for sdi (urban) 2 hyperparameters
    input<-scan(paste(datapath,region,"\\Farms\\ChickenFarmDistance.txt",sep=""))
    input<-matrix(1-min(input)+input,1,mbs)
    sdi[which(ur<7)]<<-10+input[which(ur<7)]
    sds<<-max(sdi)
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:10) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=10]
    }
    for (j in 11:sds) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=11 & Bn[[j]]<=sds]
    }
    BSumFunction <<- function() {c(sum((B[2:10]-B[1:9])^2),sum((B[12:sds]-B[11:(sds-1)])^2))}
    BRank<<-c(9,sds-11)
    Bmax<<-sds
    Bklink<<-c(rep(1,10),rep(2,sds-10))
    Bks<<-2
    Blink<<-sdi
  } else if (Bmode>12 & Bmode<18) {#1 line for animal density (rural), 1 for sdi (urban), 2 hyperparameters
    if (Bmode==13) {input<-scan(paste(datapath,region,"\\Farms\\ChickenDensity.txt",sep=""))}
    if (Bmode==14) {input<-scan(paste(datapath,region,"\\Farms\\SheepDensity.txt",sep=""))}
    if (Bmode==15) {input<-scan(paste(datapath,region,"\\Farms\\CowDensity.txt",sep=""))}
    if (Bmode==16) {input<-scan(paste(datapath,region,"\\Farms\\BeefDensity.txt",sep=""))}
    if (Bmode==17) {input<-scan(paste(datapath,region,"\\Farms\\DairyDensity.txt",sep=""))}
    for (i in 1:mbs) {
      input[i]<-1+length(which(thresholds<input[i]))
    }
    sdi[which(ur<7)]<<-10+input[which(ur<7)]
    sds<<-max(sdi)
    Bmax<<-sds
    Bklink<<-c(rep(1,10),rep(2,sds-10))
    Bks<<-2
    for (i in 1:sds) {
      cat("Group",i,"has",length(which(sdi==i)),"members.\n")
    }
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:10) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=10]
    }
    for (j in 11:sds) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=11 & Bn[[j]]<=sds]
    }
    BSumFunction <<- function() {c(sum((B[2:10]-B[1:9])^2),sum((B[12:sds]-B[11:(sds-1)])^2))}
    BRank<<-c(9,sds-11)
    Blink<<-sdi
  } else if (Bmode==18) {# single line for SDI
    Bmax<<-sds
    Bklink<<-rep(1,Bmax)
    for (i in 1:sds) {
      cat("Group",i,"has",length(which(sdi==i)),"members.\n")
    }
    B<<-matrix(rnorm(sds,0,1),1,sds)
    for (j in 1:10) {
      Bn[[j]]<<-c(j-1,j+1)
      Bn[[j]]<<-Bn[[j]][Bn[[j]]>=1 & Bn[[j]]<=10]
    }
    BSumFunction <<- function() {sum((B[2:sds]-B[1:(sds-1)])^2)}
    BRank<<-sds-1
    Blink<<-sdi            
  }
}
BSetPriors <- function(setpriors) {
  if (setpriors==0 | setpriors==1) {
    aB<<-1
    bB<<-10^-4
  } else if (setpriors==2) {
    aB<<-5
    bB<<-5/500
  } else if (setpriors==3) {
    aB<<-10^-2
    bB<<-10^-2
  }     
}
BUpdate <- function(i=0) {
  # Gibb's step to update kB
  kB<<-rgamma(Bks,aB+BRank/2,rate=(bB+BSumFunction()/2))
  # Update B
  for (j in 1:Bmax) {
    wch<-Bn[[j]]
    if (i%%2==0) {
      proposal<-rnorm(1,mean(B[wch]),(1/kB[Bklink[j]]/length(wch))^(0.5))
      ap<-BLikelihood(j,proposal,which(Blink==j))
    } else {
      proposal<-rnorm(1,B[j],sd=sigmaB)
      ap<-exp(-kB[Bklink[j]]*sum((B[wch]-proposal)^2-(B[wch]-B[j])^2)/2)*BLikelihood(j,proposal,which(Blink==j))
    }
    un<-runif(1)
    if (un<=ap) {
      B[j]<<-proposal
      acceptB[1+i%%2]<<-acceptB[1+i%%2]+1
    } else {
      rejectB[1+i%%2]<<-rejectB[1+i%%2]+1
    }
  }
}
BUpdateFlat <- function(i=0) {
  # Update B
  for (j in 1:Bmax) {
    proposal<-rnorm(1,B[j],sigmaB)
    ap<-BLikelihood(j,proposal,which(Blink==j))
    un<-runif(1)
    if (un<=ap) {
      B[j]<<-proposal
      acceptB[1+i%%2]<<-acceptB[1+i%%2]+1
    } else {
      rejectB[1+i%%2]<<-rejectB[1+i%%2]+1
    }
  }
}
BSample <- function() {
  fe<<-fe+mean(B)
  B<<-B-mean(B) 
  cat(B,"\n",file="B.txt",append=TRUE,sep=" ")
  cat(kB,"\n",file="kB.txt",append=TRUE,sep=" ")
  cat(c(t(acceptB[]/(acceptB[]+rejectB[])),"\n"),file="acceptanceB.txt",append=TRUE,sep=" ")
  acceptB<<-rep(0,2)
  rejectB<<-rep(0,2)   
  cat(BSumFunction(),"\n",file="sumB.txt",append=TRUE)
}
BRisk <- function() {rep(B[Blink],each=tps)}
BLikelihoodB <- function(j,proposal,Bmb) {prod(dpois(cases[,Bmb],rep(n[Bmb]*exp(fe+proposal),each=tps))/dpois(cases[,Bmb],rep(n[Bmb]*exp(fe+B[j]),each=tps)))}
BLikelihoodBR <- function(j,proposal,Bmb) {prod(dpois(cases[,Bmb],rep(n[Bmb],each=tps)*exp(fe+rep(R,length(Bmb))+proposal))/dpois(cases[,Bmb],rep(n[Bmb],each=tps)*exp(fe+rep(R,length(Bmb))+B[j])))}
BLikelihoodBRU <- function(j,proposal,Bmb) {prod(dpois(cases[,Bmb],rep(n[Bmb],each=tps)*exp(fe+rep(R,length(Bmb))+rep(U[Bmb],each=tps)+proposal))/dpois(cases[,Bmb],rep(n[Bmb],each=tps)*exp(fe+rep(R,length(Bmb))+rep(U[Bmb],each=tps)+B[j])))}
BAnalysis <- function() {
  m<-max(max(B),-min(B))
  if (Bmode==0 || Bmode==1) {
    plot(1:sds,t="n",xlab="SDI",ylab="B",main="UR x SDI Analysis",xlim=c(1,sds),ylim=c(-m,m))
    for (i in 1:urs) {
      lines(1:sds,B[i,],col=i)
    }
    legend(sds,0,paste("Urban/Rural",1:urs),col=c(1:urs),xjust=1,yjust=0,lty=1)
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,urs,sds)
    plot(1:sds,t="n",xlab="SDI",ylab="Relative Risk",main="UR x SDI Analysis",xlim=c(1,sds),ylim=c(0,max(rr)))
    for (i in 1:urs) {
      lines(1:sds,rr[i,],col=i)
    }
    legend(sds,0,paste("Urban/Rural",1:urs),col=c(1:urs),xjust=1,yjust=0,lty=1)
  } else if (Bmode==2) {
    ax<-0:(sds-1)
    plot(ax,t="n",xlab="Distance",ylab="B",main="Distance to chicken farm",xlim=c(0,sds-1),ylim=c(-m,m))
    lines(ax,B)
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,1,sds)
    plot(ax,t="n",xlab="Distance",ylab="Relative Risk",main="Distance to chicken farm",xlim=c(0,sds-1),ylim=c(0,max(rr)))
    lines(ax,rr)
  } else if (Bmode>=3 & Bmode<=5) {
    if (Bmode==3) {titl<-"Chicken Density"}
    if (Bmode==4) {titl<-"Sheep Density"}
    if (Bmode==5) {titl<-"Cow Density"}
    ax<-1:sds
    plot(ax,t="n",xlab="Category",ylab="B",main=titl,xlim=c(1,sds),ylim=c(-m,m))
    lines(ax,B)
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,1,sds)
    plot(ax,t="n",xlab="Category",ylab="Relative Risk",main=titl,xlim=c(1,sds),ylim=c(0,max(rr)))
    lines(ax,rr)
  } else if (Bmode==6 | Bmode==12) {
    plot(1:10,t="n",xlab="SDI",ylab="B",main="SDI Analysis",xlim=c(1,10),ylim=c(-m,m))
    lines(1:10,B[1:10])
    ax<-0:(sds-11)
    plot(ax,t="n",xlab="Distance",ylab="B",main="Distance to chicken farm",xlim=c(0,sds-11),ylim=c(-m,m))
    lines(ax,B[11:sds])
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,1,sds)
    plot(1:10,t="n",xlab="SDI",ylab="Relative Risk",main="SDI Analysis",xlim=c(1,10),ylim=c(0,max(rr[1:10])))
    lines(1:10,rr[1:10])
    plot(ax,t="n",xlab="Distance",ylab="Relative Risk",main="Distance to chicken farm",xlim=c(0,sds-11),ylim=c(0,max(rr[11:sds])))
    lines(ax,rr[11:sds])    
  } else if ((Bmode>6 & Bmode<12) | (Bmode>12 & Bmode<18)) {
    plot(1:10,t="n",xlab="SDI",ylab="B",main="SDI Analysis",xlim=c(1,10),ylim=c(-m,m))
    lines(1:10,B[1:10])
    if (Bmode==7 | Bmode==13) {titl<-"Chicken Density"}
    if (Bmode==8 | Bmode==14) {titl<-"Sheep Density"}
    if (Bmode==9 | Bmode==15) {titl<-"Cow Density"}
    if (Bmode==10 | Bmode==16) {titl<-"Beef Density"}
    if (Bmode==11 | Bmode==17) {titl<-"Dairy Density"}
    ax<-1:(sds-10)
    plot(ax,t="n",xlab="Category",ylab="B",main=titl,xlim=c(1,sds-10),ylim=c(-m,m))
    lines(ax,B[11:sds])
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,1,sds)
    plot(1:10,t="n",xlab="SDI",ylab="Relative Risk",main="SDI Analysis",xlim=c(1,10),ylim=c(0,max(rr[1:10])))
    lines(1:10,rr[1:10])
    plot(ax,t="n",xlab="Category",ylab="Relative Risk",main=titl,xlim=c(1,sds-10),ylim=c(0,max(rr[11:sds])))
    lines(ax,rr[11:sds])
  } else if (Bmode==18) {
    ax<-1:sds
    plot(ax,t="n",xlab="SDI",ylab="B",main="SDI Analysis",xlim=c(1,sds),ylim=c(-m,m))
    lines(ax,B)
    rr<-scan('relativeriskB.txt')
    rr<-matrix(rr,1,sds)
    plot(ax,t="n",xlab="SDI",ylab="Relative Risk",main="SDI Analysis",xlim=c(1,sds),ylim=c(0,max(rr)))
    lines(ax,rr)
  }  
}
BConvergence <- function() {
  plotPairs("kB",Bks)
  plotPairs("sumB",p=F)
  plotPairs("acceptanceB",length(acceptB),F,ha=T)
}
BTraces <- function() {
  if (Bmode==0 || Bmode==1) {
    B<<-matrix(plotPairs("B",Bmax,z=T),urs,sds)
  } else {
    B<<-matrix(plotPairs("B",Bmax,z=T),1,sds)
  } 
}
   