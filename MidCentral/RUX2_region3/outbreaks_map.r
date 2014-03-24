# this code produces a Google Earth KML file of outbreaks and their probabilities

# set the following as needed:

workingDir         <- paste(zhome,"$REGION$/RUX2_region$SUBREGION$",sep="") # directory containing this code
relRiskFile        <- "relativeriskU.txt"                                   # output from RUX or KHR model
outputFile         <- "relativerisk.pdf"                                    # output file
meshblockShapeFile <- paste(zhome,"$REGION$/SpatialUnits.shp",sep="")       # shapefile for spatial units
meshblocksFile     <- paste(zhome,"$REGION$/Meshblocks.txt",sep="")         # spatial units file

#library(splancs); library(maptools); ; library(sp); library(spatstat); 
library(RColorBrewer);
library(classInt);
library(maptools);

setwd(workingDir)
file.remove(outputFile)

# load in our relative risk data and meshblocks
rr <- scan(relRiskFile);
mbs <- read.table(meshblocksFile, header=F)

# setup our colours
pal <- brewer.pal(n = 9, name = "RdYlGn")
breaks <- c(-Inf, 0.2, 0.5, 0.8, 0.95, 1.05, 1.3, 2.0, 5.0, Inf)
textBreaks <- c("<0.2","0.2-0.5","0.5-0.8","0.8-0.95", "0.95-1.05","1.05-1.3","1.3-2.0", "2.0-5.0",">5.0")

# load in our shape file
map.shp <- NULL
try(map.shp <- readShapeSpatial(meshblockShapeFile))

if (!is.null(map.shp))
{
  # grab our meshblock ids (not all meshblocks are in our data file...)
  ids <- slot(map.shp, "data")[[1]]

  risk <- rep(1,length(ids));
  # for each mb in ids, check to see that there's a matching mb in our Meshblocks.txt file
  for (i in 1:length(ids))
  {
    w <- which(mbs[,1] == ids[i])
    if (length(w) == 1)
    {
      risk[i] <- rr[w];
    }
  }

  q <- classIntervals(risk, n = 7, style = "fixed", fixedBreaks = breaks)
  qColours <- findColours(q, rev(pal))

  xylims <- attr(map.shp, "bbox")

  ratio <- (xylims[2,2] - xylims[2,1]) / (xylims[1,2] - xylims[1,1])
  pdf(outputFile, width=8/sqrt(ratio), height=8*sqrt(ratio))

  par(pin = c(7.5/sqrt(ratio), sqrt(ratio) * 7.5), omi = c(0,0,0,0))
  plot(x = xylims[1,], y = xylims[2,], type = "n", xaxt = "n", yaxt = "n", xlab="", ylab="", xlim = xylims[1,], ylim = xylims[2,], cex.lab = 1.0)
#  plot(nzreg.shp, lwd=0.05, col="grey80", add=TRUE)
  plot(map.shp, lty=0, col = qColours, border="grey60", lwd=0.2, add=TRUE)
  legend(x = "topleft", fill = attr(qColours, "palette"), legend = textBreaks, bg = "white" , cex = 0.80, bty = "n")
  dev.off()
}
