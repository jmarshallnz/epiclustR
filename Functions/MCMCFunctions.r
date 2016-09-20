# Default values for model
defaults <- function() {
  zhome <- ""
  list(incR = FALSE,
       incS = FALSE,
       incU = FALSE,
       incV = FALSE,
       incW = FALSE,
       incX = FALSE,
       incB = FALSE,
       Bmode = 0,
       Xmode = 0,
       region = "",
       regionchoice = "",
       tidyup = TRUE,
       iters = 12000,
       burnin = 2000,
       samplefreq = 10,
       mbs = 10,
       maxne = 1,
       tps = 10,
       baseDeviance = 0,
       datapath = zhome,
       versioncode = "",
       dataset = ""
    )
}

Initialise <- function(region,weather="",setpriors=0) {
  params$region <<- region
  if (region == "MidCentral") {
    params$mbs <<- 1834
    params$maxne <<- 17
    params$baseDeviance <<- 21000
  } else if (region == "MidCentral07") {
    params$mbs <<- 1834
    params$maxne <<- 17
    params$baseDeviance <<- 21000
    params$tps <<- 364
  } else if (region == "MidCentral08") {
    params$mbs <<- 1743
    params$maxne <<- 17
    params$tps <<- 398
    params$baseDeviance <<- 21000    
  } else if (region == "Auckland") {
    params$mbs <<- 3156
    params$maxne <<- 19
  } else if (region == "NZ") {
    params$mbs <<- 1769
    params$maxne <<- 17
    params$baseDeviance <<- 20000
  } else if (region == "HRC") {
    params$mbs <<- 5014
    params$maxne <<- 29
    if (params$incR && params$incU) {
      params$baseDeviance <<- 103000
    } else {
      params$baseDeviance <<- 107000
    }
  } else if (region == "Christchurch") {
    params$mbs <<- 4313
    params$maxne <<- 26
    if (params$incR && params$incU) {
      params$baseDeviance <<- 103000
    } else {
      params$baseDeviance <<- 107000
    }
  } else if (region == "AAuckland") {
    params$baseDeviance <<- 270000  
    params$mbs <<- 9709
    params$maxne <<- 19
  } else if (region == "CC") {
    params$mbs <<- 1834
    params$maxne <<- 17
    params$tps <<- 98
  } else if (region == "CC\\CCInc2007") {
    params$mbs <<- 1
    params$tps <<- 137
    params$maxne <<- 1
  } else if (region == "ST") {
    params$mbs <<- 1834
    params$maxne <<- 17
    params$tps <<- 98
  } else if (region == "STPaper") {
    params$mbs <<- 1743
    params$maxne <<- 17
    if (dataset == "Rum") {
      params$tps <<- 149
    } else if (dataset == "474") {
      params$tps <<- 147
    } else {
      params$tps <<- 154
    }
  } else if (region == "STSpace") {
    params$mbs <<- 1743
    params$maxne <<- 17
    params$tps <<- 154
  } else if (region == "ST08") {
    params$mbs <<- 1743
    params$maxne <<- 17
    params$tps <<- 364-214
  } else if (region == "ST\\Inc2007") {
    params$mbs <<- 1
    params$tps <<- 137
    params$maxne <<- 1
  } else if (region == "ST\\STsOnly") {
    params$mbs <<- 1834
    params$maxne <<- 17
    params$tps <<- 98
  } else if (region == "ST\\OnlyInc2007") {
    params$mbs <<- 1
    params$tps <<- 137
    params$maxne <<- 1
  } else if (region == "MidCentralUrban") {
    params$mbs <<- 1292
    params$tps <<- 364
  } else if (region == "MidCentralRural") {
    params$mbs <<- 542
    params$tps <<- 364
  } else if (region == "ChristchurchUrban") {
    params$mbs <<- 3490
    params$tps <<- 312
  } else if (region == "ChristchurchRural") {
    params$mbs <<- 823
    params$tps <<- 312
  } else if (region == "AAucklandUrban") {
    params$mbs <<- 9126
    params$tps <<- 312
  } else if (region == "AAucklandRural") {
    params$mbs <<- 574
    params$tps <<- 312
  } else if (region == "Default") {
    cat("Called with a default region - assuming that mbs/tps/maxne etc. is set\n")
  } else {
    stop("Region not recognised.\n")
  }

  # load spatial data
  if (params$incU | (params$incB & params$Bmode<2)) {
    input<-scan(file.path(params$datapath,params$region,"Weights.GAL"))
    weight<<-matrix(0,params$mbs,params$maxne+1)
    i<-3
    for (j in 1:params$mbs) {
      k<-input[i]
      weight[k,1]<<-input[i+1]
      for (l in 1:input[i+1]) {
        weight[k,1+l]<<-input[i+1+l]
      }
      i<-i+2+input[i+1]
    }
  }

  # load case data
  load_cases <- function(params) {
    input<-scan(file.path(params$datapath,params$region,paste0("Data",params$dataset,".txt")))
    matrix(input,params$tps,params$mbs)
  }
  cases <<- load_cases(params)

  # load meshblock data
  load_spatial <- function(params, cases, neighbours) {
    input<-scan(file.path(params$datapath,params$region,"Meshblocks.txt"))
    input<-matrix(input,ncol=params$mbs)
    n <- matrix(input[2,],1,params$mbs)
    cat("n created, size", dim(n), "\n")
    sdi <- NULL; ur <- NULL;
    if (params$incB) {
      sdi <- matrix(input[3,],1,params$mbs)
      ur <- matrix(input[4,],1,params$mbs)
    }
    for (i in 1:ncol(n)) {
      if (n[i]==0 && sum(cases[,i])>0) {
        cat("Meshblock",i,"has population zero and some cases - population assumed to be one.\n")
        n[i]<-1
      }
      if (params$incB) {
        if (sdi[i]==0) {
          sdi[i] <- interpolate(sdi[neighbours[i,2:neighbours[i,1]]])
          cat("Meshblock ",i," has SDI zero - interpolated to be ",sdi[i],".\n",sep="")
        }
        if (ur[i]==0) {
          ur[i] <- interpolate(ur[neighbours[i,2:neighbours[i,1]]])
          cat("Meshblock ",i," has Urban/Rural status zero - interpolated to be ",ur[i],".\n",sep="")
        }
      }
    }
    return(list(n=n, sdi=sdi, ur=ur))
  }
  o <- load_spatial(params, cases, weights)
  n <<- o$n; sdi <<- o$sdi; ur <<- o$ur

  # load weather data
  load_weather <- function(weather, params) {
    source(file.path("Functions", "WFunctions.r"))
    input <- scan(paste(params$datapath,params$region,"\\Weather\\",weather,".txt",sep=""))
    matrix(1-min(wthr)+wthr,params$tps+2,params$mbs)
  }

  if (weather!="") {
    params$incW <<- TRUE
    wthr <<- load_weather(weather, params)
    WInitialise()
  }

  # clean up
  if (params$tidyup) {
    file.remove(file.path(params$outpath, "fixedEffects.txt"))
    file.remove(file.path(params$outpath, "deviance.txt"))
    file.remove(file.path(params$outpath, "expectedCases.txt"))
  }
  # initialise
  state <- NULL
  if (params$incR) {
    source(file.path("Functions", "RFunctions.r"))
    state <- c(state, RInitialize())
  }
  if (params$incS) { source(file.path("Functions", "SFunctions.r")) }
  if (params$incU) {
    source(file.path("Functions", "UFunctions.r"))
    state <- c(state, UInitialize())
  }
  if (params$incV) { source(file.path("Functions", "VFunctions.r")) }
  if (params$incX) {
    source(file.path("Functions", paste0("XFunctions", params$versioncode, ".r")))
    state <- c(state, XInitialise())
  }
  if (params$incB) {
    source(file.path("Functions", "BFunctions.r"))
    BInitialise()
  }
  if (params$incR) { RSetPriors(setpriors) }
  if (params$incS) { SSetPriors(setpriors) }
  if (params$incU) { USetPriors(setpriors) }
  if (params$incV) { VSetPriors(setpriors) }
  if (params$incW) { WSetPriors(setpriors) }
  if (params$incX) { XSetPriors(setpriors) }
  if (params$incB) { BSetPriors(setpriors) }
  return(state)
}

interpolate <- function(v) {
  if (length(which(v>0))==0) {
    cat("forced to 1.\n")
    return(1)
  } else {
    return(round(mean(v[v>0])))
  }
}

Sample <- function(state) {
  if (params$incR) { RSample(state) }
  if (params$incS) { SSample() }
  if (params$incU) { USample(state) }
  if (params$incV) { VSample() }
  if (params$incW) { WSample() }
  if (params$incX) { XSample(state) }
  if (params$incB) { BSample() }
  cat(state$fe,"\n",file=file.path(params$outpath, "fixedEffects.txt"),append=TRUE)
  cat(Deviance(state),"\n",file=file.path(params$outpath, "deviance.txt"),append=TRUE)
  cat(ExpectedCases(state),"\n",file=file.path(params$outpath, "expectedCases.txt"),append=TRUE)
  if (params$incX | params$incS) {cat(ExpectedCases(state,smoothed=TRUE),"\n",file=file.path(params$outpath, "smoothedCases.txt"),append=TRUE)}
}

InitAcceptance <- function(state) {
  if (params$incR) { state <- RInitAcceptance(state) }
#  if (params$incS) { SSample() }
  if (params$incU) { state <- UInitAcceptance(state) }
#  if (params$incV) { state <- VInitAcceptance(state) }
#  if (params$incW) { state <- WInitAcceptance(state) }
  if (params$incX) { state <- XInitAcceptance(state) }
#  if (params$incB) { state <- BInitAcceptance(state) }
  return(state)
}

Deviance <- function(state) {
  sum(-2*log(dpois(cases[,],ECases(state))))-params$baseDeviance
}

# expected cases per time and space
ECases <- function(state, smoothed = FALSE) {
  data <- list(cases=cases, popn=n, mbrg=mbrg, nb=weight, rgmb=wch)
  return(expected_cases(data, state, smoothed))
}

# expected number of cases per time
ExpectedCases <- function(state, smoothed = FALSE) {
  apply(matrix(ECases(state,smoothed=smoothed),params$tps,params$mbs),1,sum)
}

# This *SEEMS* to be independent of model state etc.
plotPairs<-function(variable,components=1,posteriors=TRUE,zeroCentre=FALSE,halfCentre=FALSE,useDataCentre=FALSE,matched=TRUE) {
  input <- -1
  try(input<-scan(file.path(params$outpath, paste0(variable,".txt"))),TRUE)
  if (length(input)>1) {
    st<-1+params$burnin/params$samplefreq
    en<-length(input)/components
    input<-matrix(input,components,en)
    pmean<-matrix(0,1,components)
    rrisk<-matrix(0,1,components)
    file.remove(file.path(params$outpath, paste0("posterior",variable,".txt")))
    if (components>1) { file.remove(file.path(params$outpath, paste0("relativerisk",variable,".txt"))) }
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
      write.table(t(pmean),file.path(params$outpath, paste0("posterior",variable,".txt")),row.names=F,col.names=F)
      if (components>1) {write.table(t(rrisk),file.path(params$outpath, paste0("relativerisk",variable,".txt")),row.names=F,col.names=F)}      
    }
    return(pmean)
  }
}

Convergence <- function(state)
{
  try(dev.off(),TRUE)
  pdf(paper="a4r",width=11,height=7,file=file.path(params$outpath, "Convergence.pdf"))
  par(mfrow=c(1,2))
  fe<-plotPairs("fixedEffects",z=F)
  pMeanDeviance<-plotPairs("deviance",useDataCentre=T)
  if (params$incR) { RConvergence() }
  if (params$incS) { SConvergence() }
  if (params$incU) { UConvergence(state) }
  if (params$incV) { VConvergence() }
  if (params$incW) { WConvergence() }
  if (params$incX) { XConvergence(state) }
  if (params$incB) { BConvergence() }
  dev.off()
  #
  pdf(paper="a4r",width=11,height=7,file=file.path(params$outpath, "Traces.pdf"))
  par(mfrow=c(3,4))
  if (params$incW) { WTraces() }
  if (params$incV) { VTraces() }
  if (params$incR) { RTraces() }
  if (params$incS) { STraces() }
  if (params$incU) { UTraces() }
  if (params$incB) { BTraces() }
  dev.off()
  #
  dic_file <- file.path(params$outpath, "DIC.txt")
  file.remove(dic_file)
  cat("Posterior Mean Deviance: ",params$baseDeviance+pMeanDeviance,"\n",file=dic_file,sep="",append=TRUE)
  cat("Effective Parameters:    ",pMeanDeviance-Deviance(state),"\n",file=dic_file,sep="",append=TRUE)
  cat("DIC:                     ",params$baseDeviance+2*pMeanDeviance-Deviance(state),"\n",file=dic_file,sep="",append=TRUE)
}

Analysis <- function()
{
  pdf(paper="a4r",width=11,height=7,file=file.path(params$outpath, "Analysis.pdf"))
  par(mfrow=c(1,1))
  if (params$incX) { XAnalysis() }
  if (params$incW) { WAnalysis() }
  if (params$incV) { VAnalysis() }
  if (params$incB) { BAnalysis() }
  if (params$incR) { RAnalysis() }
  if (params$incS) { SAnalysis() }
  if (params$incX) { XAnalysis2() }
  dev.off()
}

Continue <- function() {
  st <<- -1
  if (params$incR) { R <<- Upload("R",tps); kR <<- Upload("kR",1) }
  if (params$incU) { U <<- Upload("U",mbs); kU <<- Upload("kU",1) }
  if (params$incW) { W <<- Upload("W",cs); kW <<- Upload("kW",1) }
  if (params$incB) {
    B <<- Upload("B", Bmax)
    kB <<- Upload("kB", Bks)
    if (params$Bmode==0 | params$Bmode==1) { B <<- matrix(B,urs,sds) }
  }
  if (params$incX) { XContinue() }
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
