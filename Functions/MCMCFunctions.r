# Default values
incR<-FALSE
incS<-FALSE
incU<-FALSE
incV<-FALSE
incW<-FALSE
incX<-FALSE
incB<-FALSE
Bmode<-0
Xmode<-0
regionchoice<-""
tidyup<-TRUE
iters<-12000
burnin<-2000
samplefreq<-10
tps<-10
fe<--10
baseDeviance<-0
mcmcpath<-zhome
datapath<-zhome
versioncode<-""
dataset<-""
Initialise <- function(region,weather="",setpriors=0,mbdataset=dataset) {
  region<<-region
  if (region=="MidCentral") {
    mbs<<-1834
    maxne<<-17
    baseDeviance<<-21000
  } else if (region=="MidCentral07") {
    mbs<<-1834
    maxne<<-17
    baseDeviance<<-21000
    tps<<-364
  } else if (region=="MidCentral08") {
    mbs<<-1743
    maxne<<-17
    tps<<-398
    baseDeviance<<-21000    
  } else if (region=="Auckland") {
    mbs<<-3156
    maxne<<-19
  } else if (region=="NZ") {
    mbs<<-1769
    maxne<<-17
    baseDeviance<<-20000
  } else if (region=="HRC") {
    mbs<<-5014
    maxne<<-29
    if (incR && incU) {
      baseDeviance<<-103000
    } else {
      baseDeviance<<-107000
    }
  } else if (region=="Christchurch") {
    mbs<<-4313
    maxne<<-26
    if (incR && incU) {
      baseDeviance<<-103000
    } else {
      baseDeviance<<-107000
    }
  } else if (region=="AAuckland") {
    baseDeviance<<-270000  
    mbs<<-9709
    maxne<<-19
  } else if (region=="CC") {
    mbs<<-1834
    maxne<<-17
    tps<<-98
  } else if (region=="CC\\CCInc2007") {
    mbs<<-1
    tps<<-137
    maxne<<-1
  } else if (region=="ST") {
    mbs<<-1834
    maxne<<-17
    tps<<-98
  } else if (region=="STPaper") {
    mbs<<-1743
    maxne<<-17
    if (dataset=="Rum") {
      tps<<-149
    } else if (dataset=="474") {
      tps<<-147
    } else {
      tps<<-154
    }
  } else if (region=="STSpace") {
    mbs<<-1743
    maxne<<-17
    tps<<-154
  } else if (region=="ST08") {
    mbs<<-1743
    maxne<<-17
    tps<<-364-214
  } else if (region=="ST\\Inc2007") {
    mbs<<-1
    tps<<-137
    maxne<<-1
  } else if (region=="ST\\STsOnly") {
    mbs<<-1834
    maxne<<-17
    tps<<-98
  } else if (region=="ST\\OnlyInc2007") {
    mbs<<-1
    tps<<-137
    maxne<<-1
  } else if (region=="MidCentralUrban") {
    mbs<<-1292
    tps<<-364
  } else if (region=="MidCentralRural") {
    mbs<<-542
    tps<<-364
  } else if (region=="ChristchurchUrban") {
    mbs<<-3490
    tps<<-312
  } else if (region=="ChristchurchRural") {
    mbs<<-823
    tps<<-312
  } else if (region=="AAucklandUrban") {
    mbs<<-9126
    tps<<-312
  } else if (region=="AAucklandRural") {
    mbs<<-574
    tps<<-312
  } else {
    stop("Region not recognised.\n")
  }  
  # load spatial data
  if (incU | (incB & Bmode<2)) {
    input<-scan(paste(datapath,region,"\\Weights.GAL",sep=""))
    weight<<-matrix(0,mbs,maxne+1)
    i<-3
    for (j in 1:mbs) {
      k<-input[i]
      weight[k,1]<<-input[i+1]
      for (l in 1:input[i+1]) {
        weight[k,1+l]<<-input[i+1+l]
      }
      i<-i+2+input[i+1]
    }
  }
  # load case data
  input<-scan(paste(datapath,region,"\\Data",dataset,".txt",sep=""))
  cases<<-matrix(input,tps,mbs)
  # load meshblock data
  input<-scan(paste(datapath,region,"\\Meshblocks",mbdataset,".txt",sep=""))
  input<-matrix(input,ncol=mbs)
  n<<-matrix(input[2,],1,mbs)
  if (incB) {
    sdi<<-matrix(input[3,],1,mbs)
    ur<<-matrix(input[4,],1,mbs)
  }
  for (i in 1:mbs) {
    if (n[i]==0 && sum(cases[,i])>0) {
      cat("Meshblock",i,"has population zero and some cases - population assumed to be one.\n")
      n[i]<<-1
    }
    if (incB && sdi[i]==0) {
      sdi[i]<<-interpolate(sdi[weight[i,2:weight[i,1]]])
      cat("Meshblock ",i," has SDI zero - interpolated to be ",sdi[i],".\n",sep="")
    }
    if (incB && ur[i]==0) {
      ur[i]<<-interpolate(ur[weight[i,2:weight[i,1]]])
      cat("Meshblock ",i," has Urban/Rural status zero - interpolated to be ",ur[i],".\n",sep="")
    }
  }
  # load weather data
  if (weather!="") {
    incW<<-TRUE
    source(paste(mcmcpath,"Functions\\WFunctions.r",sep=""))
    wthr<<-scan(paste(datapath,region,"\\Weather\\",weather,".txt",sep=""))
    wthr<<-matrix(1-min(wthr)+wthr,tps+2,mbs)
    WInitialise()
  }
  # clean up
  if (tidyup) {
    file.remove("fixedEffects.txt")
    file.remove("deviance.txt")
    file.remove("expectedCases.txt")
  }
  # initialise
  if (incR) {source(paste(mcmcpath,"Functions\\RFunctions.r",sep=""))}
  if (incS) {source(paste(mcmcpath,"Functions\\SFunctions.r",sep=""))}
  if (incU) {source(paste(mcmcpath,"Functions\\UFunctions.r",sep=""))}
  if (incV) {source(paste(mcmcpath,"Functions\\VFunctions.r",sep=""))}
  if (incX) {
    source(paste(mcmcpath,"Functions\\XFunctions",versioncode,".r",sep=""))
    XInitialise()
  }
  if (incB) {
    source(paste(mcmcpath,"Functions\\BFunctions.r",sep=""))
    BInitialise()
  }
  if (incR) {RSetPriors(setpriors)}
  if (incS) {SSetPriors(setpriors)}
  if (incU) {USetPriors(setpriors)}
  if (incV) {VSetPriors(setpriors)}
  if (incW) {WSetPriors(setpriors)}
  if (incX) {XSetPriors(setpriors)}
  if (incB) {BSetPriors(setpriors)}
}
interpolate <- function(v) {
  if (length(which(v>0))==0) {
    cat("forced to 1.\n")
    return(1)
  } else {
    return(round(mean(v[v>0])))
  }
}
Sample <- function(i) {
  if (incR) {RSample()}
  if (incS) {SSample()}
  if (incU) {USample()}
  if (incV) {VSample()}
  if (incW) {WSample()}
  if (incX) {XSample(i)}
  if (incB) {BSample()}
  cat(fe,"\n",file="fixedEffects.txt",append=TRUE)
  cat(Deviance(),"\n",file="deviance.txt",append=TRUE)
  cat(ExpectedCases(),"\n",file="expectedCases.txt",append=TRUE)
  if (incX | incS) {cat(ExpectedCases(smoothed=TRUE),"\n",file="smoothedCases.txt",append=TRUE)}  
}
Deviance <- function() {sum(-2*log(dpois(cases[,],ECases())))-baseDeviance}
ECases <- function(smoothed=FALSE) {
  output<-rep(fe,tps*mbs)
  if (incR) {output<-output+RRisk()}
  if (incS && !smoothed) {output<-output+SRisk()}
  if (incU) {output<-output+URisk()}
  if (incV) {output<-output+VRisk()}
  if (incW) {output<-output+WRisk()}
  if (incX && !smoothed) {output<-output+XRisk()}
  if (incB) {output<-output+BRisk()}
  return(rep(n,each=tps)*exp(output))
}
ExpectedCases <- function(smoothed=FALSE) {apply(matrix(ECases(smoothed=smoothed),tps,mbs),1,sum)}
#  
plotPairs<-function(variable,components=1,posteriors=TRUE,zeroCentre=FALSE,halfCentre=FALSE,useDataCentre=FALSE,matched=TRUE) {
  input<--1
  try(input<-scan(paste(variable,".txt",sep="")),T)
  if (length(input)>1) {
    st<-1+burnin/samplefreq
    en<-length(input)/components
    input<-matrix(input,components,en)
    pmean<-matrix(0,1,components)
    rrisk<-matrix(0,1,components)
    file.remove(paste("posterior",variable,".txt",sep=""))
    if (components>1) {file.remove(paste("relativerisk",variable,".txt",sep=""))}
    if (matched) {
      m<-c(min(0,input[,st:en]),max(input[,st:en]))
      if (useDataCentre) {m<-c(min(input[,st:en]),max(input[,st:en]))}
      if (halfCentre) {m<-c(0,1)}
      if (zeroCentre) {m<-c(-max(-m[1],m[2]),max(-m[1],m[2]))}
    }
    for (i in 1:components) {
      pmean[i]<-mean(input[i,st:en])
      rrisk[i]<-mean(exp(input[i,st:en]))
      if (!matched) {
        m<-c(min(0,input[i,st:en]),max(input[i,st:en]))
        if (useDataCentre) {m<-c(min(input[i,st:en]),max(input[i,st:en]))}
        if (halfCentre) {m<-c(0,1)}
        if (zeroCentre) {m<-c(-max(-m[1],m[2]),max(-m[1],m[2]))}
      }
      hist(input[i,st:en],nclass=30,probability=TRUE,xlim=m,main=variable,xlab="")
      plot(c(st:en),input[i,st:en],t="l",main=i,ylim=m,xlab="Iteration",ylab="")
    }
    if (posteriors) {
      write.table(t(pmean),paste("posterior",variable,".txt",sep=""),row.names=F,col.names=F)
      if (components>1) {write.table(t(rrisk),paste("relativerisk",variable,".txt",sep=""),row.names=F,col.names=F)}      
    }
    return(pmean)
  }
}
Convergence <- function()
{
  try(dev.off(),T)
  pdf(paper="a4r",width=11,height=7,file="Convergence.pdf")
  par(mfrow=c(1,2))
  fe<-plotPairs("fixedEffects",z=F)
  pMeanDeviance<-plotPairs("deviance",useDataCentre=T)
  if (incR) {RConvergence()}
  if (incS) {SConvergence()}
  if (incU) {UConvergence()}
  if (incV) {VConvergence()}
  if (incW) {WConvergence()}
  if (incX) {XConvergence()}
  if (incB) {BConvergence()}
  dev.off()
  #
  pdf(paper="a4r",width=11,height=7,file="Traces.pdf")
  par(mfrow=c(3,4))
  if (incW) {WTraces()}
  if (incV) {VTraces()}
  if (incR) {RTraces()}
  if (incS) {STraces()}
  if (incU) {UTraces()}
  if (incB) {BTraces()}
  dev.off()
  #
  file.remove("DIC.txt")
  cat("Posterior Mean Deviance: ",baseDeviance+pMeanDeviance,"\n",file="DIC.txt",sep="",append=TRUE)
  cat("Effective Parameters:    ",pMeanDeviance-Deviance(),"\n",file="DIC.txt",sep="",append=TRUE)
  cat("DIC:                     ",baseDeviance+2*pMeanDeviance-Deviance(),"\n",file="DIC.txt",sep="",append=TRUE)
}

Analysis <- function()
{
  pdf(paper="a4r",width=11,height=7,file="Analysis.pdf")
  par(mfrow=c(1,1))
  if (incX) {XAnalysis()}
  if (incW) {WAnalysis()}
  if (incV) {VAnalysis()}
  if (incB) {BAnalysis()}
  if (incR) {RAnalysis()}
  if (incS) {SAnalysis()}
  if (incX) {XAnalysis2()}
  dev.off()
}
Retrospective <- function(samples) {
  file.remove("DIC.txt")
  file.remove("deviance.txt")
  file.remove("expectedCases.txt")
  iters<-samples*samplefreq
  fes<-scan("fixedEffects.txt")
  if (incR) {Rs<-matrix(scan("R.txt"),tps,samples)}
  if (incU) {Us<-matrix(scan("U.txt"),mbs,samples)}
  if (incW) {Ws<-matrix(scan("W.txt"),cs,samples)}
  for (i in 1:samples) {
    fe<<-fes[i]
    if (incR) {R<<-Rs[,i]}
    if (incU) {U<<-Us[,i]}
    if (incW) {W<<-Ws[,i]}
    cat(Deviance(),"\n",file="deviance.txt",append=TRUE)
    cat(ExpectedCases(),"\n",file="expectedCases.txt",append=TRUE)
  }   
}
Continue <- function() {
  st<<--1
  if (incR) {R<<-Upload("R",tps);kR<<-Upload("kR",1)}
  if (incU) {U<<-Upload("U",mbs);kU<<-Upload("kU",1)}
  if (incW) {W<<-Upload("W",cs);kW<<-Upload("kW",1)}
  if (incB) {
    B<<-Upload("B",Bmax)
    kB<<-Upload("kB",Bks)
    if (Bmode==0 | Bmode==1) {B<<-matrix(B,urs,sds)}
  }
  if (incX) {XContinue()}
  if (st==-2) {
    return(FALSE)
  } else {
    st<<-1+st*samplefreq
    return(TRUE)
  }    
}
Upload <- function(nam,num) {
  input<-scan(paste(nam,".txt",sep=""))
  if (st==-1) {st<<-length(input)/num}
  if (st!=-2 & st==length(input)/num) {
    input<-matrix(input,num,st)
    return(input[,st])
  } else if (st!=-2) {
    cat("File length mismatch in ",nam,".txt\n",sep="")
    st<<--2
    return(-1)
  } else {
    return(-2)
  }
}
CumulativeInitialise <- function(wk,setpriors=0) {
  tps<<-wk
  cases<<-allcases[1:wk,]#Change this to account for the 'Friday factor'.
  # initialise
  if (incR) {source(paste(mcmcpath,"Functions\\RFunctions.r",sep=""))}
  if (incS) {source(paste(mcmcpath,"Functions\\SFunctions.r",sep=""))}
  if (incU) {source(paste(mcmcpath,"Functions\\UFunctions.r",sep=""))}
  if (incV) {source(paste(mcmcpath,"Functions\\VFunctions.r",sep=""))}
  if (incX) {
    source(paste(mcmcpath,"Functions\\XFunctions",versioncode,".r",sep=""))
    XFullOutput<<-TRUE
    XInitialise()
  }
  if (incB) {
    source(paste(mcmcpath,"Functions\\BFunctions.r",sep=""))
    BInitialise()
  }
  if (incR) {RSetPriors(setpriors)}
  if (incS) {SSetPriors(setpriors)}
  if (incU) {USetPriors(setpriors)}
  if (incV) {VSetPriors(setpriors)}
  if (incW) {WSetPriors(setpriors)}
  if (incX) {XSetPriors(setpriors)}
  if (incB) {BSetPriors(setpriors)}
}
Warmup <- function(region) {
  region<<-region
  if (region=="MidCentral") {
    mbs<<-1834
    maxne<<-17
    baseDeviance<<-21000
  } else if (region=="MidCentral07") {
    mbs<<-1834
    maxne<<-17
    baseDeviance<<-21000
    tps<<-364    
  } else if (region=="Auckland") {
    mbs<<-3156
    maxne<<-19
  } else if (region=="Christchurch") {
    mbs<<-4313
    maxne<<-26
    if (incR && incU) {
      baseDeviance<<-103000
    } else {
      baseDeviance<<-107000
    }
  } else if (region=="AAuckland") {
    baseDeviance<<-270000  
    mbs<<-9709
    maxne<<-19
  } else {
    stop("Region not recognised.\n")
  }  
  # load spatial data
  input<-scan(paste(datapath,region,"\\Weights.GAL",sep=""))
  weight<<-matrix(0,mbs,maxne+1)
  i<-3
  for (j in 1:mbs) {
    k<-input[i]
    weight[k,1]<<-input[i+1]
    for (l in 1:input[i+1]) {
      weight[k,1+l]<<-input[i+1+l]
    }
    i<-i+2+input[i+1]
  }
  # load case data
  input<-scan(paste(datapath,region,"\\Data",dataset,".txt",sep=""))
  allcases<<-matrix(input,tps,mbs)
  # load meshblock data
  input<-scan(paste(datapath,region,"\\Meshblocks",dataset,".txt",sep=""))
  input<-matrix(input,5,mbs)
  n<<-matrix(input[2,],1,mbs)
  if (incB) {
    sdi<<-matrix(input[3,],1,mbs)
    ur<<-matrix(input[4,],1,mbs)
  }
  for (i in 1:mbs) {
    if (n[i]==0 && sum(allcases[,i])>0) {
      cat("Meshblock",i,"has population zero and some cases - population assumed to be one.\n")
      n[i]<<-1
    }
    if (incB && sdi[i]==0) {
      sdi[i]<<-interpolate(sdi[weight[i,2:weight[i,1]]])
      cat("Meshblock ",i," has SDI zero - interpolated to be ",sdi[i],".\n",sep="")
    }
    if (incB && ur[i]==0) {
      ur[i]<<-interpolate(ur[weight[i,2:weight[i,1]]])
      cat("Meshblock ",i," has Urban/Rural status zero - interpolated to be ",ur[i],".\n",sep="")
    }
  }
  # clean up
  if (tidyup) {
    file.remove("fixedEffects.txt")
    file.remove("deviance.txt")
    file.remove("expectedCases.txt")
  }
  # cumulative specific defaults
  span<<-13
  depth<<-13
  stwk<<-105
  enwk<<-tps
  Initialise<<-CumulativeInitialise  
}
getmcmcweeks <- function() {
  totweeks<-ceiling((enwk-stwk+1)/span)
  return(stwk-1+1:totweeks)
}
Infer <- function(wk) {
  allU<-scan("U.txt")
  samples<-allU/mbs   
  allU<-matrix(allU,mbs,samples)
  allR<-scan("R.txt")
  allR<-matrix(allR,wk,samples)
  allX<-scan("fullX.txt")
  allX<-array(allX,c(wk,rgs,samples))
  for (d in 1:depth) {
    getInitials()
    for (i in (burnin/samplefreq):samples) {
      for (j in 1:samplefreq) {
        RCumulativeUpdate(i)
        XCumulativeUpdate(i)
      }
    }   
  }     
}
    
   