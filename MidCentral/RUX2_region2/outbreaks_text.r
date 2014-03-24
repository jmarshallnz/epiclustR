# outputs an outbreaks.txt file with a summary of the most likely outbreaks
#
# set the following as needed:

tolerance          <- 0.1                                             # minimum probability to pick up
workingDir         <- paste(zhome,"$REGION$/RUX2_region$SUBREGION$",sep="")      # directory containing this code
posteriorXFile     <- "posteriorX.txt"                                # output from RUX or KHR model
meshblocksFile     <- paste(zhome,"$REGION$/Meshblocks.txt",sep="")   # Meshblocks file
regionFile         <- paste(zhome,"$REGION$/Regions.txt",sep="")      # Region file (single row for each meshblock with region number)
casesFile          <- paste(zhome,"$REGION$/Cases.txt",sep="")        # cases file
idsFile            <- paste(zhome,"$REGION$/IDs.txt",sep="")          # ID file (TODO: Merge with cases file?)
startDate          <- $BEGINTIME$                                     # date of first timepoint
outputFile         <- "outbreaks.txt"                                 # output file

setwd(workingDir)
file.remove(outputFile)

# load in our case data

cases <- scan(casesFile)
cases <- matrix(cases,2,length(cases)/2) # MB, time pairs
ids <- as.matrix(read.table(idsFile))
cases <- cbind(ids, t(cases)) # id, MB, time pairs

# load in our meshblock and region data

mb <- read.table(meshblocksFile, header=F)
regions <- scan(regionFile)
num_regions <- max(regions)

# load in our posterior outbreak probabilities

X <- scan(posteriorXFile)
X <- matrix(X, ncol=num_regions)
num_times <- nrow(X)

obrks<-matrix(0,4,num_regions) # store the average, max and maxweek outbreak info per region
for (i in 1:num_regions)
{
  values <- X[(num_times-5):num_times,i]
  obrks[1,i] = mean(values)
  obrks[2,i] = max(values)
  start = which.max(values)
  cat(start,"\n")
  while (1)
  {
    if (start == 1)
      break
    if (values[start-1] < tolerance)
      break
    start <- start - 1
  }
  end = which.max(values)
  cat(end,"\n")
  while (1)
  {
    if (end > 5)
      break
    if (values[end+1] < tolerance)
      break
    end <- end + 1
  }
  obrks[3,i] = start + num_times-5-1
  obrks[4,i] = end + num_times-5-1
}

# now grab the top 10 of these over tolerance
ord <- order(obrks[2,], decreasing=TRUE)

for (i in ord[1:10])
{
  if (obrks[2,i] < tolerance)
    break
  # find the cases corresponding to this week and region
  wks<-obrks[3,i]:obrks[4,i]
  meshblocks <- mb[which(regions==i),1]
  matching_cases <- character(0)
  for (j in meshblocks)
  {
    for (k in wks)
    {
      case <- which(cases[,2]==j & cases[,3]==k)
      if (length(case) > 0)
      {
        cat("found case",case,cases[case,1],"\n")
        matching_cases<-c(matching_cases,cases[case,1])
      }
    }
  }
  if (length(matching_cases) > 0)
  {
    str<-sprintf("Potential outbreak in region %i", i);
    cat(file=outputFile,str,"\n",rep("=",times=nchar(str)),"\n\n",sep="",append=TRUE)
    if (obrks[3,i] == obrks[4,i])
      cat(file=outputFile,"Time period:       Week",obrks[3,i],"\n",append=TRUE)
    else
      cat(file=outputFile,"Time period:       Weeks",obrks[3,i],"-",obrks[4,i],"\n",append=TRUE)
    cat(file=outputFile,"Probability:      ",obrks[2,i],"\n",append=TRUE)
    cat(file=outputFile,"Cases in outbreak:", matching_cases,"\n\n\n",append=TRUE)
  }
}


