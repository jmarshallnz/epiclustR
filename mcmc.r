library(epiclustR)
library(dplyr)
library(tidyr)

# the priors and MCMC control parameters
prior <- init_priors()

control <- init_control(thinning=50, samples = 1000, burnin=20)

spat <- read.table("Meshblocks.txt", header=TRUE)

nb   <- load_spatial_neighbours("Weights.GAL")
popn <- spat$Population
mbrg <- spat$Zone7

case <- read.csv("CaseData/cases.csv", stringsAsFactors = FALSE)

case <- case %>% filter(Meshblock06 %in% spat$Meshblock06) %>%
  mutate(ReportWeek = as.Date(ReportWeek)) %>%
  filter(Year >= 2006 & Year <= 2015)

# now compute the number of cases per meshblock per time
weeks <- seq(min(case$ReportWeek), max(case$ReportWeek), by = 7)
cases <- matrix(0, length(weeks), length(mbrg))
rownames(cases) <- as.character(weeks)
colnames(cases) <- spat$Meshblock06

tab <- case %>% group_by(ReportWeek, Meshblock06) %>% summarize(Cases = n())
week_by_mb <- tab %>% spread(Meshblock06, Cases, fill = 0)
cases[as.character(week_by_mb$ReportWeek),names(week_by_mb)[-1]] <- as.matrix(week_by_mb[,-1])

# assemble data
data <- check_data(list(cases=cases, popn=popn, mbrg=mbrg, nb=nb))
data$spat_list <- spat %>% select(ID = ID, Spatial = Meshblock06, Region = Zone7)
data$case_list <- case %>% select(CaseID = EpiSurvNumber, Spatial = Meshblock06, ReportWeek)

# fit the model
print(system.time({
mod <- fit_model(data, prior, control, seed = 1)
}))

# save model
zone7 = list(data=data, mod=mod)
save(list = c('zone7'), file='zone7.RData')
