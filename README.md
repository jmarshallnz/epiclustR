[![Build Status](https://travis-ci.org/jmarshallnz/epiclustR.svg?branch=master)](https://travis-ci.org/jmarshallnz/epiclustR)

epiclustR
=========

epiclustR is a spatio-temporal modelling tool for estimating the spatial and temporal risk of disease, and estimating the probability of an anomolous event (a potential 'outbreak') where risk is higher than would be expected given the spatial and temporal trends.

The model of disease is given by

y<sub>it</sub> = Poisson(&lambda;<sub>it</sub>)

where

&lambda;<sub>it</sub> = R<sub>t</sub> + U<sub>i</sub> + W<sub>it</sub>

where R<sub>t</sub> is a purely temporal term, U<sub>i</sub> is a purely spatial term, and W<sub>it</sub> is a spatio-temporal term.

We place Gaussian structural priors on R_t and U_i such that the risk in week t+2 is a linear extrapolation of the risk in weeks t and t+1, and that the risk in region i is the average risk in neighbouring regions.

The model for detecting outbreaks sets

W<sub>it</sub> = &beta; X<sub>it</sub>

where X<sub>it</sub> is an indicator variable, set to 1 if there's an outbreak and 0 if there isn't, that is one with probability p<sub>it</sub> and &beta; is the (average) size of each outbreak.  The prior on p<sub>it</sub> allows 1 outbreak per year per region on average.  Priors on p<sub>it</sub> then allow these outbreaks to be correlated through time (e.g. for a disease where person to person transmission occurs) or uncorrelated through time.  The default implementation is uncorrelated in time.

Data format
-----------

The code in the PrepareDatasets folder can be used to prepare the first three datasets. The data required are as follows:

Meshblocks.txt

A tab-separated file describing the spatial units.  Each row corresponds to a different spatial unit.  There are four columns: The spatial unit ID, the population at risk, and the latitude and longitude coordinates (for plotting maps).

Weights.GAL

The file describing the spatial correlation.  Each pair of rows contains the number of neighbours and then the list of neighbours (indices into the Meshblocks.txt file) for a particular spatial unit.

Regions.txt

The regions file.  Each row represents a spatial unit, with numbers according to which region that meshblock is in.  TODO: This should be merged with Meshblocks.txt.

Data.txt

The counts file.  Each row is a spatial location/time point of the number of counts in that region at that time point.  TODO: This should be generated from a CSV file of cases.

Running
-------

Run the model by placing said files into the folder for the model run and then run epiclustR.R from R.

TODO: Currently there's a bunch of manual configuration required in epiclustR.R (and other model R functions) that needs removing.

(We have a front-end for doing this, but it'd be nicer if this was done in pure R).

Data output are found in the model folder.  You can monitor progress by checking in R or checking the contents of the 'output.txt' file which will list the current iteration number.
