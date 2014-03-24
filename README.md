epiclustR
=========

epiclustR is a spatio-temporal modelling tool for estimating the spatial and temporal risk of disease, and estimating the probability of an anomolous event (a potential 'outbreak') where risk is higher than would be expected given the spatial and temporal trends.

The model of disease is given by

y_{it} = Poisson(\lambda_{it})

where

\lambda_{it} = R_t + U_i + W_{it}

where R_i is a purely temporal term, U_i is a purely spatial term, and W_{it} is a spatio-temporal term.

We place structural priors on R_t and U_i such that R_{t+2} = Normal(2R_{t+1} - R_t, \tau_R) and U_i = Normal(\sum_{j \in N(i)} U_i/|N(i)|, \tau / |N(i)|) where N(i) is the set of regions in teh neighbourhood of region i.

Data format
-----------

The data required are as follows:

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
