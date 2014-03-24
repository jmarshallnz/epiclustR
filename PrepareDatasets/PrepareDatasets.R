# Code to generate datasets needed for epiclustR from shapefiles etc.

# Files output:

# 1. Meshblocks.txt   primary spatial regions.
# 2. Weights.GAL      neighbourhood file for the primary spatial regions.
# 3. Regions.txt      epidemic regions (aggregations of primary regions).

# basic idea:
#
# epiclustR requires the primary spatial units to each have at least one neighbour (no single-unit islands)
# it then requires to know which spatial unit is adjacent to which, which epidemic region (aggregation of primary regions)
# each spatial unit belongs to and the population at risk in each primary spatial region.
#
# everything for the model is in terms of numerical identifiers 1..max(field) so Meshblocks.txt contains the mapping
# from row number to its ID (eg CAU Identifier) in the first column, population at risk in second column.
#
# Weights.GAL is a standard spatial neighbourhood file used in SpatStat.

library(spatstat)
library(maptools)
library(spdep)

working_directory         <- "~/Massey/Projects/EpiclustR/PrepareDatasets/NZ"
primary_spatial_popn_file <- "Data/SpatialData/eCAUPopn.txt"
primary_spatial_shapefile <- "Data/NZ_shapefiles/AU2006.shp"
epidemic_concordance_file <- "Data/SpatialData/eConcordance.txt"
output_directory          <- "Output2"

setwd(working_directory)

source("geo_transforms.r")

# create output directory

out_dir <- file.path(working_directory, output_directory)
if (!file.exists(out_dir))
{
  dir.create(out_dir);
}

# read in the primary spatial population file
psu <- read.csv(primary_spatial_popn_file, header=T)
head(psu)

# read in the primary spatial shapefile
mbs <- readShapeSpatial(primary_spatial_shapefile)
ids <- slot(mbs, "data")[[1]]
mbs.polys <- slot(as(mbs, "SpatialPolygons"),"polygons")
mb_coords <- matrix(0, length(ids), 2)
for (i in 1:length(ids))
{
  mb_coords[i,] <- slot(mbs.polys[[i]],"labpt")
}
mb_coords_wgs84 <- nztm2000_to_wgs84(mb_coords[,2], mb_coords[,1])

# grab the neighbours of these meshblocks
nbs <- poly2nb(mbs, queen=F)

# see which of these ids match our population file
mb_list <- data.frame(CAUID=c(1), Pop=c(1), Long=c(1), Lat=c(1));
nn_list <- data.frame(id=c(1), CAUID=c(1), Pop=c(1));
nb_to_mb <- rep(0, length(nbs))
nn <- 1;
nm <- 1;
for (i in 1:length(nbs))
{
  if (length(nbs[[i]]) > 0 && nbs[[i]][1] > 0)
  {
    row <- which(psu$CAUID06 == ids[i])
    if (length(row) > 0)
    {
      mb_list[nm,] <- c(ids[i], sum(psu$Pop2010[row]), mb_coords_wgs84[i,]);
      nb_to_mb[i] <- nm;
      nm <- nm + 1;
    }
    else
    {
      cat("Region with ID", ids[i], "not found in population file\n");
    }
  }
  else
  {
    row <- which(psu$CAUID06 == ids[i])
    if (length(row) > 0)
    {
      pop <- sum(psu$Pop2010[row])
      cat("Region with ID", ids[i], "has no neighbours, but is in population file (popn", pop,")\n");
      if (pop > 0)
      {
        nn_list[nn,] <- c(i, ids[i], pop);
        nn <- nn + 1;
      }
    }
    else
    {
      cat("Region with ID", ids[i], "has no neighbours (and not found in population file)\n");
    }
  }
}
# now check which of our primary spatial units in the population file has no matching shape data
for (i in 1:nrow(psu))
{
  row <- which(mb_list$CAUID == psu$CAUID06[i]);
  if (length(row) == 0)
    cat("Region with ID", psu$CAUID06[i], "(popn", psu$Pop2010[i], ")has no shape data\n")
}
if (nn > 1)
{
  cat("Islands with population found (model currently doesn't support these).  islands.pdf produced to illustrate\n");
  col <- rep("white", length(ids));
  col[nn_list$id] <- "red";
  pdf(file=file.path(out_dir, "islands.pdf"), width=7, height=7)
  plot(mbs, col=col, lwd=0.05);
  dev.off();
}

# read in our concordance file
concord <- read.csv(epidemic_concordance_file, header=T)

u <- unique(concord$AU2006)
if (length(u) != nrow(psu))
{
  cat("WARNING: Number of primary spatial units in our concordance file differs from our primary spatial population file")
}

# generate our region list
region_map <- rep(0, nrow(mb_list));
regions <- unique(concord$TA06);
region_names <- concord[match(regions, concord$TA06),c("TA06Name","DHBName")]

for (i in 1:nrow(mb_list))
{
  row <- which(concord$AU2006 == mb_list$CAUID[i])
  if (length(row) == 0)
  {
    cat("ERROR: primary spatial unit with id", mb_list$CAUID[i], "not found in concordance file\n");
  }
  else
  {
    region_map[i] <- which(regions == concord$TA06[row[1]])
  }
}

# write our Meshblocks.txt file
write.table(mb_list, file=file.path(out_dir, "Meshblocks.txt"), row.names=F, col.names=F)

# write our Weights.GAL file
wf <- file.path(out_dir, "Weights.GAL");
cat("0", " ", nrow(mb_list),"\n", sep="", file=wf);
nm <- 1;
for (i in 1:length(nbs))
{
  if (length(nbs[[i]]) > 0 && nbs[[i]][1] > 0)
  {
    cat(nm, " ", length(nbs[[i]]),"\n", sep="", file=wf, append=T);
    cat(nb_to_mb[nbs[[i]]], sep=" ", file=wf, append=T);
    cat("\n", file=wf, append=T)
    nm <- nm + 1;
  }
}

# write our Regions.txt file
write.table(region_map, file=file.path(out_dir, "Regions.txt"), row.names=F, col.names=F)
write.table(data.frame(number = 1:length(region_names), region = region_names$TA06Name, bigregion = region_names$DHBName), file=file.path(out_dir, "RegionNames.txt"), row.names=F)
