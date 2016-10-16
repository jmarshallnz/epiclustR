NULL

ssapply <- function(x, fun, ...) {
  simplify2array(lapply(x, function(xx, fun, ...) { simplify2array(lapply(xx, fun, ...)) }, fun = fun, ...))
}

extract_variable <- function(x, variable) {
  x[[variable]]
}

#' Plot the temporal fit of the epiclustR model
#'
#' @param mod the model to plot
#' @param data the data for the model (for plotting cases)
#' @param spatial the spatial unit(s) to plot the trend for (defaults to NULL, which is all units)
#' @param region the region(s) to plot the trend for (defaults to NULL, which is all regions)
#' @param col_cases,col_trend,col_outbreaks The colour to use for each component
#' @param lwd_cases,lwd_trend,lwd_outbreaks The line width to use for each component
#' @export
plot_temporal <- function(mod, data, spatial = NULL, region = NULL,
                          col_cases="grey80", col_trend = "black", col_outbreaks = "red",
                          lwd_cases = 1, lwd_trend = 2, lwd_outbreaks = 1) {
  u = seq_along(data$mbrg)
  regen = FALSE
  if (!is.null(spatial)) {
    # TODO: Lookup which region number(s) is given
    u = spatial
    regen = TRUE
  }
  if (!is.null(region)) {
    u = unlist(data$rgmb[region])
    regen = TRUE
  }
  if (regen) {
    ecases = ssapply(mod, cases_per_time, data=data, smoothed=FALSE, spatial=u)
    scases = ssapply(mod, cases_per_time, data=data, smoothed=TRUE, spatial=u)
  } else {
    ecases = ssapply(mod, extract_variable, 'ecases')
    scases = ssapply(mod, extract_variable, 'scases')
  }
  weeks = as.Date(rownames(data$cases))

  # grab out just the cases for the particular region
  ucases = data$cases[,u,drop=FALSE]

  plot(weeks, apply(ucases, 1, sum), type='l', col=col_cases, xaxs='i', xlab='', lwd=lwd_cases,
       ylab='Cases per week')
  lines(weeks, apply(ecases, 1, mean), col=col_outbreaks, lwd=lwd_outbreaks)
  lines(weeks, apply(scases, 1, mean), col=col_trend, lwd=lwd_trend)
}

#' Plot traces from MCMC chains to assess convergence
#' 
#' @param mod the model to plot
#' @export
plot_traces <- function(mod) {
  kR = ssapply(mod, extract_variable, 'kR')
  kU = ssapply(mod, extract_variable, 'kU')
  pX = ssapply(mod, extract_variable, 'pX')
  plot_trace <- function(x, ylab) {
    plot(NULL, xlim = c(1, nrow(x)), ylim = range(x), xlab="", ylab=ylab, type="n", xaxs="i")
    invisible(lapply(1:ncol(x), function(i) { lines(x[,i], col=i) }))
    text(nrow(x), max(x), paste("ESS =", round(ess(x),1)), adj=c(1.1, 0.5))
  }
  plot_trace(kR, expression(tau[R]))
  plot_trace(kU, expression(tau[U]))
  plot_trace(pX, expression(p))
}

#' Plot outbreaks from an epiclustR model
#' 
#' @param mod the model to plot
#' @param data the data for the model (for plotting cases)
#' @param region which region to plot. Defaults to NULL (all regions).
#' @param ... other parameters passed to the plotting functions.
#' @export
plot_region <- function(mod, data, region = NULL, ...) {
  # pull the outbreak binary variables out of our model posterior
  X = ssapply(mod, extract_variable, 'X')
  mX = apply(X, 1:2, mean)
  if (!is.null(region))
    mX = mX[,region,drop=FALSE]
  mX = apply(mX, 1, max)

  # grab the weeks from the case data
  weeks = as.Date(rownames(data$cases))

  # filter our case data down to just the regions specified
  u = seq_along(data$mbrg)
  if (!is.null(region)) {
    u = unlist(data$rgmb[region])
  }
  cases = apply(data$cases[,u,drop=FALSE], 1, sum)

  # and do the plot
  plot(weeks, mX, type="h", ylim=c(-1,1), xlim=range(weeks), xaxs="i", xlab="", ylab="", yaxt="n", ...)
  axis(2, labels=c(0,0.5,1), at=c(0,0.5,1))
  case_axis <- pretty(cases)
  axis(4, labels=case_axis, at=-1*case_axis/max(case_axis))
  segments(weeks, 0, weeks, -cases/max(case_axis), col="grey70")
  legend(x = "topleft", fill= c("black", "grey70"), legend=c("Outbreak probability", "Observed cases"), bg = "white", bty="n")
}
