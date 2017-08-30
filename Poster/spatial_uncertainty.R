# spatial stuff
library(RColorBrewer)
library(classInt)
library(maptools)
library(dplyr)

brewerBrBG9 <- rev(c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A"))

breaks <- c(-Inf, 0.2, 0.5, 0.8, 0.95, 1.05, 1.3, 2.0, 5.0, Inf)

get_col <- function(c, pal) {
  func <- colorRamp(pal, space="Lab")
  rgb(func(c), maxColorValue = 255)
}

shapeFile <- "Poster/maps/AU2006.shp"

ax_col <- "grey20"

load("Poster/data/TA2006.Rdata"); TA <- new

scases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases')
dim(scases) <- c(dim(scases)[1], prod(dim(scases)[2:3]))
weeks = as.Date(rownames(TA$data$cases))

# spatial uncertainty of posterior for pre/post 2008
U <- epiclustR:::ssapply(TA$mod, epiclustR:::extract_spatial, data=TA$data)
dim(U) <- c(dim(U)[1:2], prod(dim(U)[3:4]))

library(maptools)
library(sf)
mc <- sf::st_read('Poster/maps/midcentral_phu.shp')
hb <- sf::st_read('Poster/maps/hawkes_bay.shp')
nz <- sf::st_read('Poster/maps/NZ_region-NZTM2000.shp', layer='REGION')

# convert phu to the correct projection
nzgd1949.proj = '+proj=nzmg +lat_0=-41 +lon_0=173 +x_0=2510000 +y_0=6023150 +ellps=intl +datum=nzgd49 +units=m +no_defs'
nzgd2000.proj = '+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996 +x_0=1600000 +y_0=10000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

st_crs(mc) <- nzgd1949.proj
mc <- st_transform(mc, nzgd2000.proj)

# union up the major boundaries
mc_u <- st_union(mc)
hb_u  <- st_union(hb)

# right, now compute the risk surfaces to plot
load('Poster/data/mc/TA2006.Rdata'); dat_mc <- new
load('Poster/data/hb/HNArea.Rdata'); dat_hb <- new

# combine chains into a single list of iterations
U_mc <- epiclustR:::ssapply(dat_mc$mod, epiclustR:::extract_spatial, data=dat_mc$data)
dim(U_mc) <- c(dim(U_mc)[1:2], prod(dim(U_mc)[3:4]))

U_hb <- epiclustR:::ssapply(dat_hb$mod, epiclustR:::extract_spatial, data=dat_hb$data)
dim(U_hb) <- c(dim(U_hb)[1:2], prod(dim(U_hb)[3:4]))

# break into appropriate levels
U <- c(U_mc, U_hb)
num_cols <- 101
up_bks = quantile(U[U > 0], seq(1,num_cols,by=2)/num_cols)
lo_bks = quantile(U[U < 0], seq(num_cols-1,0,by=-2)/num_cols)
bks = round(pmax(up_bks, abs(lo_bks)), 2)
breaks = c(-rev(bks), bks)

alpha <- function(col, x = 0.5) {
  rgb(t(col2rgb(col)), alpha = x * 255, maxColorValue = 255)
}
brewerBrBG9 <- c("#01665E","#35978F","#80CDC1","#C7EAE5","#F5F5F5","#F6E8C3","#DFC27D","#BF812D","#8C510A")
#cols <- alpha(brewerBrBG9, 0.7)
cols <- colorRampPalette(brewerBrBG9, space='Lab')(num_cols)

num_alpha <- 10
breaks_alpha <- seq(0,1,length.out=num_alpha+1)

uncertain_func <- function(x) {
  y <- 1-(max(sum(x > 0), sum(x < 0)) / length(x)*2 - 1)
  # those more uncertain, make more so (y == 1)
  sqrt(y)
}

unc_mc <- apply(U_mc, 1:2, uncertain_func)
unc_mc <- apply(unc_mc, 2, function(x) { as.numeric(cut(x, breaks=breaks_alpha, include.lowest = TRUE)) })

unc_hb <- apply(U_hb, 1:2, uncertain_func)
unc_hb <- apply(unc_hb, 2, function(x) { as.numeric(cut(x, breaks=breaks_alpha, include.lowest = TRUE)) })

cols_unc <- lapply(cols, function(x) { colorRampPalette(c(x, 'grey50'), space='Lab')(num_alpha+1) })
alpha_cols <- simplify2array(cols_unc)

# mean risk
risk_mc <- apply(U_mc, 1:2, mean)
risk_hb <- apply(U_hb, 1:2, mean)

spat_mc <- cbind(dat_mc$data$spat_list, Risk=risk_mc, Uncertainty=unc_mc)
spat_mc$Spatial = as.numeric(as.character(spat_mc$Spatial))

mc_risk <- mc %>%
  dplyr::left_join(spat_mc, by=c('MB06' = 'Spatial'))

spat_hb <- cbind(dat_hb$data$spat_list, Risk=risk_hb, Uncertainty=unc_hb)
spat_hb$Spatial = as.numeric(as.character(spat_hb$Spatial))

hb_risk <- hb %>%
  dplyr::left_join(spat_hb, by=c('MB2013' = 'Spatial'))

plot_map <- function(map, risk_var, uncertainty_var, ...) {
  col_val <- cut(map[[risk_var]], breaks = breaks)
  map_col <- alpha(alpha_cols[cbind(map[[uncertainty_var]], col_val)], 1)
  
  plot(map[1], col=map_col, ...)
}

havelock <- st_union(hb_risk %>% filter(Region == "Havelock"))

# Doing the plot, ya'll
xlim <- 1695000 + c(0,330000)
ylim <- 5400000 + c(0,330000 * 46.8/33.1)
pdf("Poster/figures/test.pdf")
par(mar=c(0,0,0,0))
plot(nz["REGION"], col='grey70', border=NA, main='')
rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border='black', lwd=2)
dev.off()

xlim <- 1700000 + c(0,330000)
ylim <- 5405000 + c(0,330000 * 46.8/33.1)
pdf("Poster/figures/main.pdf", width=33.1, height=46.8)
par(mar=c(0,0,0,0))
plot(nz["REGION"], col='grey70', border=NA, xlim=xlim, ylim=ylim, main='')
plot_map(mc_risk, 'Risk.2', 'Uncertainty.2', lwd=0.02, border=NA, add=TRUE)
plot_map(hb_risk, 'Risk', 'Uncertainty', lwd=0.02, border=NA, add=TRUE)

plot(mc_u, col=NA, border='black', lwd=2, add=TRUE)
plot(hb_u, col=NA, border='black', lwd=2, add=TRUE)
plot(havelock, col=NA, border='red3', lwd=3, add=TRUE)
xlim <- 1920000 + c(0, 20000)
ylim <- 5599000 + c(0, 25000)
rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border='black', lwd=2)

xlim <- 1817500 + c(0,9000)
ylim <- 5525000 + c(0,9000)

rect(xlim[1], ylim[1], xlim[2], ylim[2], col=NA, border='black', lwd=2)

dev.off()

# now do some zoomed in ones for the rest of the poster

# First MidCentral
xlim <- 1817500 + c(0,9000)
ylim <- 5525000 + c(0,9000)
pdf("Poster/figures/mc2006.pdf", width=7.5, height=7.5)
par(mar=c(0,0,0,0))
plot_map(mc_risk, 'Risk.1', 'Uncertainty.1', lwd=0.02, border='grey80', xlim=xlim, ylim=ylim, main='')
dev.off()
pdf("Poster/figures/mc2016.pdf", width=7.5, height=7.5)
par(mar=c(0,0,0,0))
plot_map(mc_risk, 'Risk.2', 'Uncertainty.2', lwd=0.02, border='grey80', xlim=xlim, ylim=ylim, main='')
dev.off()

# Now Hawkes Bay
xlim <- c(1920791, 1942290)
ylim <- c(5597710, 5625823)
plot_map(hb_risk, 'Risk', 'Uncertainty', lwd=0.02, border='grey80', xlim=xlim, ylim=ylim, main='')

pdf("Poster/figures/hn_hastings.pdf", width=7.5, height=10)
par(mar=c(0,0,0,0))
xlim <- 1920000 + c(0, 18750)
ylim <- 5599000 + c(0, 25000)
plot_map(hb_risk, 'Risk', 'Uncertainty', lwd=0.02, border='grey80', xlim=xlim, ylim=ylim, main='')
plot(havelock, col=NA, border='red3', lwd=3, add=TRUE)
dev.off()

# Legend, you bloody legend!
super_sat <- alpha_cols[1,seq(1,101,by=10)]
pdf("Poster/figures/legend.pdf", width=10, height=2)
par(mar=c(3,1.1,0.2,7), mgp=c(2,0.5,0.2), tcl=-0.3, cex=1.3)
plot(c(0,101),c(0,11),type='n', axes=FALSE, xlab="", ylab="", xaxs='i', yaxs='i')

lapply(1:11, function(x) { rect(0:100,x-1,1:101,x,col=alpha_cols[x,],border=alpha_cols[x,]) })
sw1 <- seq(1,41,by=10)
sw2 <- seq(62,102,by=10)
br <- c(breaks[sw1],0,breaks[sw2])
myround <- function(x) {
  nc <- ifelse(x < 0.05, 3, ifelse(x < 0.6, 2, ifelse(x > 2, 0, 1)))
  unlist(lapply(seq_along(x), function(y) { as.character(round(x[y], nc[y])) }))
}
axis(1, labels=myround(exp(br)), at=c(sw1-1,50.5,sw2-1))
title(xlab="Relative Risk")
axis(4, labels=c("Low uncertainty", "High uncertainty"), at=c(0.5,10.5), line=-0.2, pos=NA, tick=FALSE, las=1)
#title(ylab="Uncertainty")
dev.off()

#pdf("Poster/figures/hn_napier.pdf", width=6, height=6)
#ar(mar=c(0,0,0,0))
#xlim <- c(1928962, 1938802) + c(0,-1000)
#ylim <- c(5616509, )
#plot_map(hb_risk, 'Risk', 'Uncertainty', lwd=0.02, border='grey80', xlim=xlim, ylim=ylim, main='')
#dev.off()

# Now the temporal outbreak graphs
plot_temporal <- function(weeks, scases, ecases, week, ax_col='grey30') {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, ecases, ylim=c(0,12), type='l', col="red", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, scases, col="black", lwd=2)
  rect(weeks[week]-14, 0, weeks[week]+14, 12, col=alpha("steelblue", 0.7), border=NA)
}

scases_mc = apply(epiclustR:::ssapply(dat_mc$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecases_mc = apply(epiclustR:::ssapply(dat_mc$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
weeks_mc = as.Date(rownames(dat_mc$data$cases))

# Find the Raw Milk outbreak and Pahiatua outbreaks
wch_rawmilk <- which(ecases_mc - scases_mc > 0.3)[5]
wch_pahiatua <- which(ecases_mc - scases_mc > 0.4)[1]

ecases_mc[wch_rawmilk] <- scases_mc[wch_rawmilk]+1 # Dodgy extra highlight???

# do a plot with these
pdf("Poster/figures/mc_temporal.pdf", width=22, height=5)
plot_temporal(weeks_mc, scases_mc, ecases_mc, week=NA)
axis(2)
dev.off()

# Now Hawke's Bay
scases_hb = apply(epiclustR:::ssapply(dat_hb$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecases_hb = apply(epiclustR:::ssapply(dat_hb$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
weeks_hb = as.Date(rownames(dat_hb$data$cases))

extract_temporal <- function(mod, var) {
  v = epiclustR:::ssapply(mod, epiclustR:::extract_variable, var)
  dim(v) <- c(dim(v)[1], prod(dim(v)[2:3]))
  v
}
scases_hb = extract_temporal(dat_hb$mod, 'scases')
ecases_hb = apply(extract_temporal(dat_hb$mod, 'ecases'),1,median)
# extract out the median trend
df <- data.frame(Week = weeks_hb, Cases = apply(dat_hb$data$cases, 1, sum), t(apply(scases_hb, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))), Outbreaks=ecases_hb)
names(df)[1:5+2] <- c("Min", "LQ", "Median", "UQ", "Max")

library(ggplot2)

pdf("Poster/figures/hb_temporal.pdf", width=17, height=5)
cols = c(C= '#0000001f', B= '#0000002f', D='red')
ggplot(df) + 
  geom_line(aes(x=Week,y=Cases, col='A')) +
  geom_ribbon(aes(x=Week,ymin=Min,ymax=Max, fill='C')) +
  geom_ribbon(aes(x=Week,ymin=LQ,ymax=UQ, fill='B')) +
#  geom_line(aes(x=Week,y=Outbreaks), col='red') +
  geom_ribbon(aes(x=Week,ymin=Median,ymax=Outbreaks, fill='D'), size=0.5, col='red') +
  geom_line(aes(x=Week,y=Median, col='B'), size=1) +
  #  geom_ribbon(aes(x=Week,ymin=Median,ymax=Median,fill='A', col='A'), size=1) +
  theme_bw(base_size=20) + scale_x_date(name="", expand=c(0,0), breaks=as.Date(paste0(2010:2016,'-01-01')), labels=2010:2016) +
  scale_y_continuous(name="Cases per week", limits=c(0,15)) +
  scale_colour_manual(name="", labels=c("Cases   ", "Expected cases"), values=c(A='#0000002f', B='black'),
                      guide=guide_legend(order=1)) +
  scale_fill_manual(name="",labels=c("50% CI   ", "95% CI   ", "Outbreaks"), values=cols,
                    guide=guide_legend(override.aes=list(col=NA), order=2)) +
  #  guides(colour=guide_legend(override.aes = list(fill=cols[1:3]))) +
  theme(axis.title.x=element_blank(),
        legend.position = c(0.5,0.95), legend.direction='horizontal', legend.box='horizontal')
#  guides(fill=guide_legend(keywidth=0.5, keyheight=0.1, default.unit="inch"))
dev.off()

scases_mc = extract_temporal(dat_mc$mod, 'scases')
ecases_mc = apply(extract_temporal(dat_mc$mod, 'ecases'),1,median)
# extract out the median trend
df <- data.frame(Week = weeks_mc, Cases = apply(dat_mc$data$cases, 1, sum), t(apply(scases_mc, 1, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))), Outbreaks=ecases_mc)
names(df)[1:5+2] <- c("Min", "LQ", "Median", "UQ", "Max")

library(ggplot2)

pdf("Poster/figures/mc_temporal.pdf", width=23.8, height=5)
cols = c(C= '#0000001f', B= '#0000002f', D='red')
ggplot(df) + 
  geom_line(aes(x=Week,y=Cases, col='A')) +
  geom_ribbon(aes(x=Week,ymin=Min,ymax=Max, fill='C')) +
  geom_ribbon(aes(x=Week,ymin=LQ,ymax=UQ, fill='B')) +
  #  geom_line(aes(x=Week,y=Outbreaks), col='red') +
  geom_ribbon(aes(x=Week,ymin=Median,ymax=Outbreaks, fill='D'), size=0.5, col='red') +
  geom_line(aes(x=Week,y=Median, col='B'), size=1) +
  #  geom_ribbon(aes(x=Week,ymin=Median,ymax=Median,fill='A', col='A'), size=1) +
  theme_bw(base_size=20) + scale_x_date(name="", expand=c(0,0), breaks=as.Date(paste0(2006:2016,'-01-01')), labels=2006:2016) +
  scale_y_continuous(name="Cases per week", limits=c(0,15)) +
  scale_colour_manual(name="", labels=c("Cases   ", "Expected cases"), values=c(A='#0000002f', B='black'),
                      guide=guide_legend(order=1)) +
  scale_fill_manual(name="",labels=c("50% CI   ", "95% CI   ", "Outbreaks"), values=cols,
                    guide=guide_legend(override.aes=list(col=NA), order=2)) +
  #  guides(colour=guide_legend(override.aes = list(fill=cols[1:3]))) +
  theme(axis.title.x=element_blank(),
        legend.position = c(0.5,0.95), legend.direction='horizontal', legend.box='horizontal')
#  guides(fill=guide_legend(keywidth=0.5, keyheight=0.1, default.unit="inch"))
dev.off()

# Outbreak attribution
library(forcats)
outbreaks <- rbind(read.csv("Poster/data/pahiatua.csv") %>% mutate(Outbreak="Pahiatua"),
                   read.csv("Poster/data/rawmilk.csv") %>% mutate(Outbreak="Raw milk")) %>%
  mutate(Source = fct_relevel(Source, "Poultry", "Cattle", "Sheep", "Other"))

library(ggplot2)
pdf("Poster/figures/outbreak_attribution_pahiatua.pdf", width=6, height=6)
ggplot(outbreaks %>% filter(Outbreak == "Pahiatua")) +
  geom_violin(aes(Source, p, fill=Source), scale='width') +
  theme_bw(base_size=15) +
  scale_fill_manual(values = c("plum4", "steelblue", "steelblue2", "brown"), guide = "none") +
  xlab("") +
  coord_flip() +
  scale_y_continuous(name = "Attribution of human cases", labels = scales::percent) +
  scale_x_discrete(limits=rev(levels(outbreaks$Source)))
dev.off()

pdf("Poster/figures/outbreak_attribution_rawmilk.pdf", width=6, height=6)
ggplot(outbreaks %>% filter(Outbreak == "Raw milk")) +
  geom_violin(aes(Source, p, fill=Source), scale='width') +
  theme_bw(base_size=15) +
  scale_fill_manual(values = c("plum4", "steelblue", "steelblue2", "brown"), guide = "none") +
  xlab("") +
  coord_flip() +
  scale_y_continuous(name = "Attribution of human cases", labels = scales::percent) +
  scale_x_discrete(limits=rev(levels(outbreaks$Source)))
dev.off()




















# This code is used to union up to the various outbreak regions
au_reg = phu@data %>% left_join(AU$data$spat_list, by=c("MB06"="Spatial"))
au_shp = unionSpatialPolygons(phu, sprintf("%02d", au_reg$Region))
ta_reg = phu@data %>% left_join(TA$data$spat_list, by=c("MB06"="Spatial"))
ta_shp = unionSpatialPolygons(phu, sprintf("%02d", ta_reg$Region))






# outbreak plot...
mean_case_rate <- function(mod, data) {
  ecases <- matrix(0, nrow(data$cases), ncol(data$cases))
  rownames(ecases) <- rownames(data$cases)
  colnames(ecases) <- colnames(data$cases)
  for (i in seq_along(mod)) {
    cat("up to chain", i, "\n")
    for (j in seq_along(mod[[i]])) {
      ecases <- ecases + epiclustR::log_case_rate(data, mod[[i]][[j]], smoothed=TRUE)
    }
  }
  ecases <- ecases / (length(mod) * length(mod[[1]]))
  ecases - mean(ecases)
}

Uta <- mean_case_rate(TA$mod, TA$data)
Uau <- mean_case_rate(AU$mod, AU$data)

# grab outreak data
roll_outbreak_probs <- function(mod, window=20) {
  X <- epiclustR:::ssapply(mod, epiclustR:::extract_variable, 'X')
  mX = apply(X, 1:2, mean)
  
  # make it a rolling window
  mR <- mX
  for (i in 1:window)
    mR[-(1:i),] <- mR[-(1:i),] + mX[1:(nrow(mX)-i),] * exp(-5/window*i)
  mR
}
Rta <- roll_outbreak_probs(TA$mod)
Rau <- roll_outbreak_probs(AU$mod)


# plot these through time
plot_ob_map <- function(sl, u, x, bbox = sp::bbox(phu)) {
  spat_risk <- cbind(sl, Risk=u)
  outbreaks <- data.frame(Region = 1:length(x), P=x)

  map_dat <- phu@data %>%
    dplyr::left_join(spat_risk, by=c('MB06' = 'Spatial')) %>%
    dplyr::left_join(outbreaks, by="Region")

  alpha2 <- function(col, x = 0.5) {
    rgb(t(col2rgb(col))*x + (1-x) * 255, maxColorValue = 255)
  }

  vals <- cut(map_dat$Risk, breaks = breaks)
  map_col <- alpha2(cols[vals], 0.7)

  red <- function(x) {
    red_func <- colorRamp(brewer.pal(9, "Reds"), space="Lab")
    rgb(red_func(pmin(x,1)), maxColorValue = 255)
  }
  sp::plot(phu, col=red(map_dat$P), lwd=0.02, border='grey80', xlim=bbox[1,], ylim=bbox[2,])
#  sp::plot(phu, add=TRUE, col=red(map_dat$P), border=NA)
}


# plot these through time
plot_ob_map2 <- function(map, x, bbox = sp::bbox(phu)) {
  red <- function(x) {
    red_func <- colorRamp(brewer.pal(9, "Reds")[-1], space="Lab")
    rgb(red_func(pmin(x,1)), maxColorValue = 255)
  }
  sp::plot(map, col=red(x), lwd=0.5, border='grey50', xlim5623977=bbox[1,], ylim=bbox[2,])
  #  sp::plot(phu, add=TRUE, col=red(map_dat$P), border=NA)
}

plot_temporal <- function(weeks, scases, ecases, week) {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, ecases, ylim=c(0,12), type='l', col="red", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, scases, col="black", lwd=2)
  rect(weeks[week]-14, 0, weeks[week]+14, 12, col=alpha("steelblue", 0.7), border=NA)
}

library(maptools)
phu <- readShapeSpatial('maps/midcentral_phu')
AU$data$spat_list$Region <- as.numeric(as.factor(AU$data$spat_list$Region))

num_cols <- 101
up_bks = quantile(U[U > 0], seq(1,num_cols,by=2)/num_cols)
lo_bks = quantile(U[U < 0], seq(num_cols-1,0,by=-2)/num_cols)
bks = round(pmax(up_bks, abs(lo_bks)), 3)
breaks = c(-rev(bks), bks)

library(dplyr)

for (i in 1:nrow(Rta)) {#seq_len(nrow(U))) {
  png(sprintf("png/outbreak_map_nt_%04d.png", i), width=960, height=480, bg="transparent")
  layout(matrix(1:4, 2, 2), heights=c(4,1), widths=c(1,1))
  par(mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", nrow(Rta), "\n")
  # tail mX off over a bunch of frames
  plot_ob_map2(ta_shp, Rta[i,])
  par(mai=c(0.5,0.416,0,0.416))
  plot_temporal(weeks, scasesTA, ecasesTA, i)
  par(mai=c(0.1,0.5,0.1,0.5))
  plot_ob_map2(au_shp, Rau[i,], bbox=matrix(c(2717946,6078900,2754889,6102322), 2))
  par(mai=c(0.5,0.416,0,0.416))
  plot_temporal(weeks, scasesAU, ecasesAU, i)
  # plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
  dev.off()
}

system("~/ffmpeg -y -r 10 -i png/outbreak_map_nt_%04d.png -b:v 1M video/outbreak_map.mp4")

system("convert -delay 10 -loop 0 -dispose background outbreak_map_nt_*.png outbreak_map2.gif")
system("convert outbreak_map2.gif -transparent white outbreak_map.gif")

animation::saveVideo(
for (i in 1:nrow(Uta)) {#seq_len(nrow(U))) {
#  png(sprintf("outbreak_map_nt_%04d.png", i), width=960, height=480, bg="transparent")
  layout(matrix(1:4, 2, 2), heights=c(4,1), widths=c(1,1))
  par(mai=c(0,0,0,0), bg="#FFFFFF")
  cat("plotting", i, "of", nrow(U), "\n")
  # tail mX off over a bunch of frames
  plot_ob_map(TA$data$spat_list, Uta[i,], Rta[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesTA, ecasesTA, i)
  par(mai=c(0,0,0,0))
  plot_ob_map(AU$data$spat_list, Uau[i,], Rau[i,])
  par(mai=c(0.5,0.5,0,0.5))
  plot_temporal(weeks, scasesAU, ecasesAU, i)
  # plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
#  dev.off()
}, video.name="outbreak_map.mp4", ffmpeg="avconv", ani.width=960, ani.height=480, interval=1/24)


system("convert -delay 10 -loop 0 -dispose background outbreak_map_nt_*.png outbreak_map2.gif")
system("convert outbreak_map2.gif -transparent white outbreak_map.gif")


# TODO: could do just one, combined plot using orange + purple instead of red
i <- 137
png(sprintf("test.png", i), width=960, height=480, bg="transparent")
layout(matrix(c(1,3,2,3), 2, 2), heights=c(4,1))
par(mai=c(0,0,0,0), bg="#FFFFFF")
cat("plotting", i, "of", nrow(U), "\n")
# tail mX off over a bunch of frames
plot_ob_map(TA$data$spat_list, Uta[i,], Rta[i,])
plot_ob_map(AU$data$spat_list, Uau[i,], Rau[i,])
par(mai=c(0.5,0.5,0,0.5))
plot_temporal2(weeks, scasesTA, ecasesTA, ecasesAU, i)
# plot_map(U[i,], bbox=matrix(c(2727946,6086900,2734889,6094322), 2))
dev.off()

plot_temporal2 <- function(weeks, scases, ecases1, ecases2, week) {
  par(mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  PuOr <- brewer.pal(11, "PuOr")[c(2,10)]
  plot(weeks, ecases1, ylim=c(0,12), type='n', col=PuOr[1], xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='', axes=FALSE, col.lab=ax_col)
  #  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  segments(weeks, y0=scases, y1=ecases1, col=PuOr[1], lwd=2)
  segments(weeks, y0=scases, y1=ecases2, col=PuOr[2], lwd=2)
#  lines(weeks, ecases2, col=puor[2], lwd=1)
  lines(weeks, scases, col="black", lwd=2)
  rect(weeks[week]-15, 0, weeks[week]+15, 12, col=alpha("steelblue", 0.7), border=NA)
}

# TODO: Add the temporal plot underneath this...
scasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesTA = apply(epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases'), 1, median)
scasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'scases'), 1, median)
ecasesAU = apply(epiclustR:::ssapply(AU$mod, epiclustR:::extract_variable, 'ecases'), 1, median)

plot_temporal(weeks, scasesTA, ecasesTA)
plot_temporal(weeks, scasesAU, ecasesAU)



# compute the expected cases for all observations
scases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'scases')
ecases = epiclustR:::ssapply(TA$mod, epiclustR:::extract_variable, 'ecases')
dim(scases) <- c(dim(scases)[1], prod(dim(scases)[2:3]))
dim(ecases) <- c(dim(ecases)[1], prod(dim(ecases)[2:3]))
weeks = as.Date(rownames(TA$data$cases))

# generate some interpolated random normals
set.seed(3)
d = sample(1:ncol(scases), 20)
d = c(d, d[1])

o_s = apply(scases[,d], 1, function(x) { spline(seq_along(x), x, method="periodic", xout=seq(1, length(x), length.out=length(x)*12+1))$y[-1]})
o_e = apply(ecases[,d], 1, function(x) { spline(seq_along(x), x, method="periodic", xout=seq(1, length(x), length.out=length(x)*12+1))$y[-1]})

med <- apply(scases, 1, median)
lci <- apply(scases, 1, quantile, 0.25)
uci <- apply(scases, 1, quantile, 0.75)

med_e <- apply(ecases, 1, median)

for (i in seq_len(nrow(o))) {
  png(sprintf("ob_fit%04d.png", i), width=960, height=480)
  par(mar=c(3,3,0.5,1), mgp=c(2,.7,0), tck=-.015, bg="#FFFFFF")
  plot(weeks, apply(TA$data$cases, 1, sum), ylim=c(0,20), type='l', col="grey70", xaxs='i', yaxs='i', xlab='', lwd=1,
       ylab='Cases', axes=FALSE, col.lab=ax_col)
  axis(2, col=ax_col, col.axis=ax_col, las=1, cex.axis=0.8)
  axis(1, col=ax_col, col.axis=ax_col, at=as.Date(paste0(2006:2017,"-01-01")), labels=rep("",12))
  mtext(2006:2016, side=1, col=ax_col, at=as.Date(paste0(2006:2016,"-07-01")), line=0.5)
  lines(weeks, o_e[i,], lwd=1, col=alpha("red",0.5))
  lines(weeks, o_s[i,], lwd=1, col=alpha("black",1))
  lines(weeks, med_e, col="red", lwd=1)
  lines(weeks, med, col="black", lwd=2)
#  polygon(c(weeks, rev(weeks)), c(lci, rev(uci)), col=alpha(col[1], 0.3), border=NA)
  dev.off()
}

system("convert -delay 4 -loop 0 -dispose background ob_fit*.png outbreak_fit2.gif")
system("convert outbreak_fit2.gif -transparent white outbreak_fit.gif")
