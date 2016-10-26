#' @importFrom grDevices col2rgb rgb
#' @importFrom stats quantile
#' @importFrom dplyr %>%
#' @importFrom sp plot
NULL

#' Extract spatial trend(s) from an epiclustR model state
#'
#' @param state the model state
#' @param data the data used to fit the model
#' @return The spatial risk adjusted as needed for the differing means through time per period.
extract_spatial <- function(state, data) {
  #' pull out U and R from the posterior
  risk <- state$U
  for (i in seq_along(data$p2t)) {
    risk[,i] = risk[,i] + mean(state$R[data$p2t[[i]]])
  }
  risk
}

#' Plot spatial trends from an epiclustR model
#' 
#' @param mod the model to plot
#' @param data the data used to fit the model
#' @param shapefile the shape file to use for the map
#' @param period which time periods(s) to plot. Default to NULL (all times)
#' @param threshold threshold at which to show outbreaks. Any outbreak with probability over this threshold will
#' be highlighted. Defaults to 0.2
#' @param bbox The bounding box to use for the map, useful for zooming in on portions. Should be a matrix. Default to NULL
#' which uses the bounding box of the shapefile.
#' @export
plot_spatial <- function(mod, data, shapefile, period=NULL, threshold=0.2, bbox=NULL) {
  map <- maptools::readShapeSpatial(shapefile)

  spat_risk <- apply(ssapply(mod, extract_spatial, data=data), 1:2, mean)

  # average out the risk over all time periods
  spat_period <- 1:ncol(spat_risk)
  if (!is.null(period)) {
    spat_period = spat_period[period]
    spat_risk   = spat_risk[,period,drop=FALSE]
  }
  # average across the risk (weighted sum by period length)
  weights = lengths(spat_period)
  spat_risk = spat_risk %*% weights / sum(weights)

  # grab out the outbreak data
  X <- ssapply(mod, extract_variable, 'X')
  mX = apply(X, 1:2, mean)

  # filter down the period
  if (!is.null(period)) {
    # TODO: make this work with dates
    mX = mX[unlist(data$p2t[period]), ,drop=FALSE]
  }

  # filter down the threshold
  outbreaks <- data.frame(Region = 1:ncol(mX), P = apply(mX, 2, max))
  outbreaks <- outbreaks[outbreaks$P >= threshold,,drop=FALSE]

  # form the data.frame
  spat_risk <- cbind(data$spat_list, Risk=spat_risk)

  map_dat <- map@data %>%
    dplyr::left_join(spat_risk, by=c('MB06' = 'Spatial')) %>%
    dplyr::left_join(outbreaks, by='Region')

  # figure out some break-points for colours
  up_bks = quantile(map_dat$Risk[map_dat$Risk > 0], seq(1,9,by=2)/9)
  lo_bks = quantile(map_dat$Risk[map_dat$Risk < 0], seq(8,0,by=-2)/9)
  bks = round(pmax(up_bks, abs(lo_bks)), 2)

  # define the colours
  alpha <- function(col, x = 0.5) {
    rgb(t(col2rgb(col)), alpha = x * 255, maxColorValue = 255)
  }
  brewerBrBG9 <- c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A")
  cols <- alpha(brewerBrBG9, 0.7)
  vals <- cut(map_dat$Risk, breaks = c(-bks, bks))
  map_col <- cols[vals]

  if (is.null(bbox) || !is.matrix(bbox))
    bbox = sp::bbox(map)

  # plot the map
  sp::plot(map, col=map_col, lwd=0.02, border='grey80', xlim=bbox[1,], ylim=bbox[2,])

  # overlay the outbreaks
  red <- function(x) {
    rgb(1, 0, 0, ifelse(is.na(x), 0, x))
  }
  plot(map, add=TRUE, col=red(map_dat$P), border=NA)
  #  xy <- t(simplify2array(lapply(map@polygons, function(x) { x@labpt })))
  #  text(xy[!is.na(map_dat$P),1], xy[!is.na(map_dat$P),2], map_dat$Region[!is.na(map_dat$P)])
  legend('topright', legend=c(levels(vals), 'Outbreak'), fill = c(cols, red(1)), cex=0.6)
}