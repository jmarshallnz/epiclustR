# this code produces a Google Earth KML file of outbreaks and their probabilities

# set the following as needed:

workingDir         <- paste(zhome,"$REGION$/RUX2_region$SUBREGION$",sep="") # directory containing this code
casesFile          <- paste(zhome,"$REGION$/Data.txt",sep="")
regionFile         <- paste(zhome,"$REGION$/Regions$SUBREGION$.txt",sep="")      # Region file (single row for each meshblock with region number)
regionNamesFile    <- paste(zhome,"$REGION$/RegionNames.txt",sep="")  # Region number:name:bigname mapping
temporalFile       <- "smoothedCases.txt"
spatioTemporalFile <- "expectedCases.txt"
outputFile         <- "outbreaks_temporal.pdf"                              # output file
startDate          <- "$BEGINTIME$"                               # date of first timepoint

setwd(workingDir)
file.remove(outputFile)

# ignore the burnin samples
burn_in_samples <- $BURNIN$/$SAMPLEFREQ$

expectedCasesRU <- read.table(temporalFile)
expectedCasesRU <- apply(expectedCasesRU[-(1:burn_in_samples),],2,mean)

num_times <- length(expectedCasesRU)

cases <- matrix(scan(casesFile), nrow=num_times)
observedCases <- apply(cases,1,sum)

expectedCasesRUX <- read.table(spatioTemporalFile)
expectedCasesRUX <- apply(expectedCasesRUX[-(1:burn_in_samples),],2,mean)

X <- read.table("posteriorX.txt")
regions <- scan(regionFile)
if (file.exists(regionNamesFile))
{
  region_names <- read.table(regionNamesFile, header=T)
} else
{
  region_names <- data.frame(number=1:num_regions, region=1:num_regions, bigregion=1:num_regions);
}
num_regions <- max(regions)

title <- "Temporal trend for $REGION$"

# find where to place our year labels
year_labels <- NULL
for (t in 1:num_times)
{
  d <- as.POSIXlt(as.Date(startDate) + 7*(t-1))
  week <- floor(d$yday/7)
  if (week == 25)
    year_labels <- rbind(year_labels,c(t,d$year+1900))
}

# find the first year
first_week <- floor(as.POSIXlt(as.Date(startDate))$yday/7)
num_years <- (num_times+52) %/% 52

pdf(outputFile, paper="a4r", width=11, height=8)
plot(1:num_times, t="n", xlab="", ylab="Cases", main=title, xlim=c(1,num_times), xaxt="n", ylim=c(0,max(observedCases, expectedCasesRU, expectedCasesRUX)))
lines(1:num_times, observedCases,col=3)
lines(1:num_times, expectedCasesRU,col=2)
lines(1:num_times, expectedCasesRUX,col=4)

axis(1, at=seq(-first_week,num_years*52-first_week,by=52), labels = FALSE)
for (i in 1:nrow(year_labels))
  mtext(year_labels[i,2], side=1, line=1, at=year_labels[i,1])

legend(x = "topright", fill = c(3,2,4), legend = c("Cases", "Trend", "Trend with outbreaks"), bg = "white" , bty = "n")

# now output each region

cases_per_region <- matrix(0, num_times, num_regions)
for (j in 1:num_regions)
{
  wch <- which(regions==j)
  if (length(wch) > 1)
  {
    cases_per_region[,j] <- apply(cases[,which(regions==j)],1,sum)
  } else
  {
    cases_per_region[,j] <- cases[,which(regions==j)]
  }
}

par(mfrow=c(2, 2))
reg <- 1:num_regions

max_cases_per_region <- round(quantile(cases_per_region[,reg], 0.99) / 3 + 0.5)*3
for (r in reg)
{
  m <- max(max_cases_per_region, cases_per_region[,r])
  title <- paste(region_names$region[r],", ",region_names$bigregion[r], sep="")
  plot(1:num_times, t="n", xlab="", xaxt="n", yaxt="n", ylab="", main=title, xlim=c(1,num_times), ylim=c(-1,1))
  axis(2, labels=c(0,0.5,1), at=c(0,0.5,1))
  segments(1:num_times, 0, 1:num_times, as.numeric(X[r,]), col="black")
  case_axis <- seq(0,max_cases_per_region,length.out=4)
  axis(4, labels=case_axis, at=-1*case_axis/m)
  wh <- which(cases_per_region[,r] > 0)
  segments(wh, 0, wh, -cases_per_region[wh,r]/m, col="green3")

  axis(1, at=seq(-first_week,num_years*52-first_week,by=52), labels = FALSE)
  axis(1, at=seq(-first_week,num_years*52-first_week,by=52/12), labels = FALSE, tck=-.02)
  for (i in 1:nrow(year_labels))
    mtext(year_labels[i,2], side=1, line=1, at=year_labels[i,1], cex=0.9)
  legend(x = "topleft", fill= c("black", "green3"), legend=c("Outbreak probability", "Observed cases"), bg = "white", bty="n")
}

dev.off()
