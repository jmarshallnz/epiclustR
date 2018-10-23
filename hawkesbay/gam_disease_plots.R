#----Packages----
library("tidyverse")
library("sf")
library("maptools")
library("mgcv")
library("visreg")
library("broom")
library("patchwork")

#----Functions----

#put into friday -> report week sum
to_week_ending_friday <- function(x) {
  x[wday(x) == 7] <- x[wday(x) == 7] + 6 # Saturday to the next week
  wday(x) <- 6 # to Friday
  x
}

#Subset by disease
subset_disease <- function(x, disease){
  x %>%
    ungroup() %>%
    filter(Disease == disease) %>%
    mutate(ReportWeek = as.Date(ReportWeek), Outbreak = ifelse(Outbreak.count.total >= 1, 1, 0))
}

#Subset by disease
subset_region <- function(x, TA.2013){
  start.date <- as.Date("10/01/1997", "%d/%m/%Y")
  x %>% filter(TA2013 == TA.2013) %>%
  mutate(Days = as.numeric(ReportWeek-start.date))}

#run gam model
enteric.gam <- function(data, knots){
  gam(Total_Cases ~ s(Days, k = knots), offset = log(Total_Population),
                family = quasipoisson, data = data)}

#produce output df for plotting
enteric.output <- function(x, y){
  #x = subsetted region df
  #y = augmented gam df
  x %>%
    bind_cols(y) %>%
    mutate(lwr = .se.fit*(-1.96), upr = .se.fit*1.96, fit.case = 100000*(exp(.fitted)), CI.up = 100000*(exp(.fitted + upr)),
           CI.lw = 100000*(exp(.fitted + lwr)))}


count_outbreaks <- function(ob) {
  # run through and accumulate up
  in_ob = FALSE
  ob_num = numeric(length(ob))
  cur_ob = 1
  for (i in seq_along(ob)) {
    if (in_ob) {
      if (ob[i] == 0) {
        # finishing an outbreak, increment the outbreak counter for the next one
        cur_ob = cur_ob + 1
        in_ob = FALSE
      }
    } else {
      if (ob[i] > 0) {
        # starting an outbreak
        in_ob = TRUE
      }
    }
    ob_num[i] <- ifelse(in_ob, cur_ob, 0)
  }
  ob_num
}

#GAM and plot -by-disease
plot_gam_disease <- function(all.data, enteric.disease, central.k, hastings.k, havelock.k, napier.k, wairoa.k,
                             y.lim, smooth=TRUE) {
  
  x <- subset_disease(all.data, enteric.disease)
  
  #subset to each region and create "days" column
  central <- subset_region(x, "Central Hawke's Bay District")
  hastings <- subset_region(x, "Hastings District")
  havelock <- subset_region(x, "Havelock North")
  napier <- subset_region(x, "Napier City")
  wairoa <- subset_region(x, "Wairoa District")
  
  #separate gam on each region
  central.gam <- enteric.gam(central, knots = central.k)
  hastings.gam <- enteric.gam(hastings, knots = hastings.k)
  havelock.gam <- enteric.gam(havelock, knots = havelock.k)
  napier.gam <- enteric.gam(napier, knots = napier.k)
  wairoa.gam <- enteric.gam(wairoa, knots = wairoa.k)
  
  #augment to extract gam info
  central.gam.aug <- augment(central.gam)
  hastings.gam.aug <- augment(hastings.gam)
  havelock.gam.aug <- augment(havelock.gam)
  napier.gam.aug <- augment(napier.gam)
  wairoa.gam.aug <- augment(wairoa.gam)
  
  central.disease.output <- enteric.output(central, central.gam.aug)
  hastings.disease.output <- enteric.output(hastings, hastings.gam.aug)
  havelock.disease.output <- enteric.output(havelock, havelock.gam.aug)
  napier.disease.output <- enteric.output(napier, napier.gam.aug)
  wairoa.disease.output <- enteric.output(wairoa, wairoa.gam.aug)
  
  output <- rbind(central.disease.output, hastings.disease.output, havelock.disease.output,
                  napier.disease.output, wairoa.disease.output)
  
  outbreaks = output %>% mutate(ReportWeek = as.Date(ReportWeek)) %>%
    group_by(TA2013) %>% arrange(ReportWeek) %>% mutate(outbreakNum = count_outbreaks(Outbreak.count.total)) %>%
    filter(outbreakNum > 0) %>%
    group_by(TA2013, outbreakNum) %>% summarise(n(), start=min(ReportWeek), end=max(ReportWeek)+7)

  g1 = ggplot() +
    scale_x_date(labels = date_format("%Y-%m-%d"), date_labels = "%Y", expand = c(0,0), date_breaks = "1 year") +
    labs(x = "Year", y = "Cases/100,000 Population") +
    geom_line(data = output, aes(ReportWeek, case_by_population), colour = "grey65", size=0.4)
  if (nrow(outbreaks) > 0) {
    g1 = g1 + geom_rect(data=outbreaks, mapping=aes(xmin=start-3, ymin=0, xmax=end+3, ymax=y.lim), col=NA, fill="red", alpha=0.4)
  }
  g1 + geom_ribbon(data=output, aes(x=ReportWeek, ymin = CI.lw, ymax = CI.up), fill = "steelblue", alpha = 0.5) +
    geom_line(data = output, aes(x=ReportWeek, fit.case), size = 0.8) +
    coord_cartesian(ylim = c(0, y.lim), expand = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~TA2013, ncol= 1) + ggtitle(paste(enteric.disease, "notifications"))
}




#GAM and plot - by-region
plot_gam_region <- function(all.data, enteric.region, campy.k, crypto.k, giardia.k, lepto.k, salmonella.k, ecoli.k,
                            y.lim) {
  
  x <- subset_region(all.data, enteric.region)
  
  #subset to each disease and create "days" column
  campy <- subset_disease(x, "Campylobacteriosis")
  crypto <- subset_disease(x, "Cryptosporidiosis")
  giardia <- subset_disease(x, "Giardiasis")
  lepto <- subset_disease(x, "Leptospirosis")
  salmonella <- subset_disease(x, "Salmonellosis")
  ecoli <- subset_disease(x, "VTEC/STEC infection")
  
  #separate gam on each region
  campy.gam <- enteric.gam(campy, campy.k)
  crypto.gam <- enteric.gam(crypto, crypto.k)
  giardia.gam <- enteric.gam(giardia, giardia.k)
  lepto.gam <- enteric.gam(lepto, lepto.k)
  salmonella.gam <- enteric.gam(salmonella, salmonella.k)
  ecoli.gam <- enteric.gam(ecoli, ecoli.k)
  
  #augment to extract gam info
  campy.gam.aug <- augment(campy.gam)
  crypto.gam.aug <- augment(crypto.gam)
  giardia.gam.aug <- augment(giardia.gam)
  lepto.gam.aug <- augment(lepto.gam)
  salmonella.gam.aug <- augment(salmonella.gam) 
  ecoli.gam.aug <- augment(ecoli.gam)
  
  campy.disease.output <- enteric.output(campy, campy.gam.aug)
  crypto.disease.output <- enteric.output(crypto, crypto.gam.aug)
  giardia.disease.output <- enteric.output(giardia, giardia.gam.aug)
  lepto.disease.output <- enteric.output(lepto, lepto.gam.aug)
  salmonella.disease.output <- enteric.output(salmonella, salmonella.gam.aug)
  ecoli.disease.output <- enteric.output(ecoli, ecoli.gam.aug)
  
  output <- rbind(campy.disease.output, crypto.disease.output, giardia.disease.output,
                  lepto.disease.output, salmonella.disease.output, ecoli.disease.output)
  
  outbreaks = output %>% mutate(ReportWeek = as.Date(ReportWeek)) %>%
    group_by(Disease) %>% arrange(ReportWeek) %>% mutate(outbreakNum = count_outbreaks(Outbreak.count.total)) %>%
    filter(outbreakNum > 0) %>%
    group_by(Disease, outbreakNum) %>% summarise(n(), start=min(ReportWeek), end=max(ReportWeek)+7)

  g1 = ggplot() +
    scale_x_date(labels = date_format("%Y-%m-%d"), date_labels = "%Y", expand = c(0,0), date_breaks = "1 year") +
    labs(x = "Year", y = "Cases/100,000 Population") +
    geom_line(data = output, aes(ReportWeek, case_by_population), colour = "grey65", size=0.4)
  if (nrow(outbreaks) > 0) {
    g1 = g1 + geom_rect(data=outbreaks, mapping=aes(xmin=start-3, ymin=0, xmax=end+3, ymax=y.lim), col=NA, fill="red", alpha=0.4)
  }
  g1 + geom_ribbon(data=output, aes(x=ReportWeek, ymin = CI.lw, ymax = CI.up), fill = "steelblue", alpha = 0.5) +
    geom_line(data = output, aes(ReportWeek, fit.case), size = 0.8) +
    coord_cartesian(ylim = c(0, y.lim), expand = FALSE) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    facet_wrap(~Disease, ncol= 1) + ggtitle(paste(enteric.region, "notifications"))
}
  
  
  
  
  
