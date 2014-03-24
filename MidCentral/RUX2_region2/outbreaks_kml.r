# this code produces a Google Earth KML file of outbreaks and their probabilities

# set the following as needed:

tolerance          <- 0.1                                                   # minimum probability to pick up
workingDir         <- paste(zhome,"$REGION$/RUX2_region$SUBREGION$",sep="") # directory containing this code
posteriorXFile     <- "posteriorX.txt"                                      # output from RUX or KHR model
meshblocksFile     <- paste(zhome,"$REGION$/Meshblocks.txt",sep="")         # Meshblocks file
regionFile         <- paste(zhome,"$REGION$/Regions$SUBREGION$.txt",sep="")            # Region file (single row for each meshblock with region number)
regionNamesFile    <- paste(zhome,"$REGION$/RegionNames.txt",sep="")  # Region number:name:bigname mapping
casesFile          <- paste(zhome,"$REGION$/Cases.txt",sep="")        # cases file
idsFile            <- paste(zhome,"$REGION$/IDs.txt",sep="")          # ID file (TODO: Merge with cases file?)
startDate          <- "$BEGINTIME$"                                   # date of first timepoint
outputFile         <- "outbreaks.kml"                                 # output file

setwd(workingDir)
file.remove(outputFile)

# load in libraries our KML generation code...

library(RColorBrewer)
source(paste(zhome,"Functions/output_kml.r",sep=""))


# load in our case data

cases <- scan(casesFile)
cases <- matrix(cases,2,length(cases)/2) # MB, time pairs
ids <- as.matrix(read.table(idsFile))
cases <- cbind(ids, t(cases)) # id, MB, time pairs

# load in our meshblock and region data

mb <- read.table(meshblocksFile, header=F)
regions <- scan(regionFile)
num_regions <- max(regions)
if (file.exists(regionNamesFile))
{
	region_names <- read.table(regionNamesFile, header=T)
} else
{
	region_names <- data.frame(number=1:num_regions, region=1:num_regions, bigregion=1:num_regions);
}

# load in our posterior outbreak probabilities

X <- scan(posteriorXFile)
X <- matrix(X, ncol=num_regions)
num_times <- nrow(X)

tolerance <- quantile(X,c(0.75,0.95,0.98,0.99))

# generate some colours suitable for KML
colour <- sprintf("FF%s", gsub("#(..)(..)(..)","\\3\\2\\1",rev(brewer.pal(n = 5, name = "RdYlGn"))[2:5]))

#
# find_cases
#
# given a time and region numbers returns a list of cases found in that region at that time point
# by looking up meshblocks in the region and cases that match that meshblock and time point.
#
# returns a list of triplets (case_id, mb_number, time_number)

find_cases <- function(t, r, cases, regions, mbs_id)
{
  # find corresponding meshblocks for this region
  mb_num <- which(regions == r);
  matching_cases <- data.frame(ids = "", mbnum = 0, numcases = 0, stringsAsFactors=FALSE)[FALSE,]
  for (j in mb_num)
  {
    # find cases matching this time and this mb
    case <- which(cases[,2]==mbs_id[j] & cases[,3]==t)
    if (length(case) > 0)
    {
      ids <- paste(cases[case,1],collapse=", ")
      matching_cases[nrow(matching_cases)+1,1] <- ids;
      matching_cases[nrow(matching_cases),2:3] <- c(j,length(case))
    }
  }
  return(matching_cases)
}

write_kml_header(outputFile, "Epidemics")
for (t in 1:num_times)
{
  cat("Week", t, "\n")
  for (r in 1:num_regions)
  {
    if (X[t,r] > tolerance[1])
    {
      # find which cases this outbreak corresponds to (i.e. find cases at time t in region r)
      outbreak_cases <- find_cases(t, r, cases, regions, mb[,1])

      if (sum(outbreak_cases$numcases) > 1)
      {
        cat(sum(outbreak_cases$numcases), "outbreak cases for time", t, "region", r, "\n")

        # add a row to our output list for each of these cases
        regName <- sprintf("%s, %s",region_names$region[r], region_names$bigregion[r]);
        write_epidemic_header(outputFile, as.Date(startDate) + 7*(t-1), regName, sum(outbreak_cases$numcases), X[t,r]);

        # TODO: This will cause (spatial) dupes - is this an issue?

        for (i in 1:nrow(outbreak_cases))
        {
          case <- outbreak_cases[i,]
          cat("case", i, "is", as.character(case), "\n")
          write_epidemic_case(outputFile, as.Date(startDate) + 7*(t-1), as.Date(startDate) + 7*t-1, mb[case$mbnum,3], mb[case$mbnum,4], case$ids, mb[case$mbnum,1], X[t,r], colour[max(which(X[t,r] > tolerance))], case$numcases)
        }
        write_epidemic_footer(outputFile)
      }
    }
  }
}
write_kml_footer(outputFile)


