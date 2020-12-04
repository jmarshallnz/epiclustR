<!-- badges: start -->
[![R build status](https://github.com/jmarshallnz/epiclustR/workflows/R-CMD-check/badge.svg)](https://github.com/jmarshallnz/epiclustR/actions)
<!-- badges: end -->

epiclustR
=========

epiclustR is a spatio-temporal modelling tool for estimating the spatial and temporal risk of disease, and estimating the probability of an anomolous event (a potential 'outbreak') where risk is higher than would be expected given the spatial and temporal trends.

<a href="http://jmarshallnz.github.io/talks/video/spatial_fit.mp4">
   <img src='spatial_fit.png' width='960' height='480' />
</a>

The model of disease is given by

y<sub>it</sub> = Poisson(&lambda;<sub>it</sub>)

where

&lambda;<sub>it</sub> = R<sub>t</sub> + U<sub>i</sub> + W<sub>it</sub>

where R<sub>t</sub> is a purely temporal term, U<sub>i</sub> is a purely spatial term, and W<sub>it</sub> is a spatio-temporal term.

We place Gaussian structural priors on R<sub>t</sub> and U<sub>i</sub> such that the risk in week t+2 is a linear extrapolation of the risk in weeks t and t+1, and that the risk in region i is the average risk in neighbouring regions.

The model for detecting outbreaks sets

W<sub>it</sub> = &beta; X<sub>it</sub>

where X<sub>it</sub> is an indicator variable, set to 1 if there's an outbreak and 0 if there isn't, that is one with probability p<sub>it</sub> and &beta; is the (average) size of each outbreak.  The prior on p<sub>it</sub> allows 1 outbreak per year per region on average.  Priors on p<sub>it</sub> then allow these outbreaks to be correlated through time (e.g. for a disease where person to person transmission occurs) or uncorrelated through time.  The default implementation is uncorrelated in time.

<a href="http://jmarshallnz.github.io/talks/video/temporal_fit.mp4">
   <img src='temporal_fit.png' width='960' height='480' />
</a>

Data format
-----------

The model requires three data files:
 - A file describing spatial units and regions.
 - A spatial neighbourhood mapping.
 - A file containing cases with week and spatial unit.

Running
-------

Use `load_spatial_neighbours` to read in the spatial neighbours file. All other files can
be read in using `read.csv` or similar.

Run `check_data` to generate the model-specific data structures.

Configure the model with `init_priors` and `init_control`.

Fit the model using `fit_model`.

You can then visualise the fitted model with `plot_model`, `plot_temporal`, `plot_spatial`, and `plot_outbreaks`.

Outbreaks may also be tabulated with `table_outbreaks`.
