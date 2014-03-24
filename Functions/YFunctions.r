yrs<-6
tps<-52
allCases<-cases
getYear <- function(year) {
  cases<<-allCases[((year-1)*tps+1):(year*tps),]
}
openAnalysis <- function() {
  try(dev.off(),T)
  pdf(paper="a4r",width=11,height=7,file="Convergence.pdf")
  par(mfrow=c(1,2))
  file.remove("DIC.txt")
}
YAnalysis <- function(year) {
  fe<-plotPairs("fixedEffects",z=F)
  pMeanDeviance<-plotPairs("deviance",useDataCentre=T)
  if (incR) {RConvergence()}
  if (incU) {UConvergence()}
  if (incV) {VConvergence()}
  if (incW) {WConvergence()}
  if (incX) {XConvergence()}
  if (incB) {BConvergence()}
  cat("Year ",year,":\n\n",file="DIC.txt",sep="",append=TRUE)
  cat("Posterior Mean Deviance: ",baseDeviance+pMeanDeviance,"\n",file="DIC.txt",sep="",append=TRUE)
  cat("Effective Parameters:    ",pMeanDeviance-Deviance(),"\n",file="DIC.txt",sep="",append=TRUE)
  cat("DIC:                     ",baseDeviance+2*pMeanDeviance-Deviance(),"\n",file="DIC.txt",sep="",append=TRUE)
  input<-scan("U.txt")
  st<-1+burnin/samplefreq
  en<-length(input)/mbs
  input<-matrix(input,mbs,en)
  pmean<-matrix(0,1,mbs)
  rrisk<-matrix(0,1,mbs)
  file.remove(paste("posteriorU",year,".txt",sep=""))
  file.remove(paste("relativeriskU",year,".txt",sep=""))
  for (i in 1:mbs) {
    pmean[i]<-mean(input[i,st:en])
    rrisk[i]<-mean(exp(input[i,st:en]))
  }
  write.table(t(pmean),paste("posteriorU",year,".txt",sep=""),row.names=F,col.names=F)
  write.table(t(rrisk),paste("relativeriskU",year,".txt",sep=""),row.names=F,col.names=F)
  file.remove("deviance.txt")
  file.remove("expectedCases.txt")
  file.remove("fixedEffects.txt")
  file.remove("posteriordeviance.txt")
  file.remove("posteriorkU.txt")
  file.remove("posteriorfixedEffects.txt")
}

