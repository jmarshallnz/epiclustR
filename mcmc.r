library(epiclustR)
library(dplyr)
library(tidyr)

# the priors and MCMC control parameters
prior <- init_priors()

control <- init_control(thinning=50, samples = 1000, burnin=20)

spat <- read.table("Meshblocks.txt", header=TRUE)

library(meshblocknz)
spat <- spat %>% left_join(mb2006, by=c("Meshblock06" = "MB2006")) %>%
  mutate(TA2006 = as.numeric(TA2006), TA2006 = ifelse(TA2006 == 43, 42, TA2006) - 38)

nb   <- load_spatial_neighbours("Weights.GAL")
popn <- spat$Population
mbrg <- spat$TA2006

case <- read.csv("CaseData/cases.csv", stringsAsFactors = FALSE)

case <- case %>% filter(Meshblock06 %in% spat$Meshblock06) %>%
  mutate(ReportWeek = as.Date(ReportWeek)) %>%
  filter(Year >= 2006 & Year <= 2016)

# now compute the number of cases per meshblock per time
weeks <- seq(min(case$ReportWeek), max(case$ReportWeek), by = 7)
cases <- matrix(0, length(weeks), length(mbrg))
rownames(cases) <- as.character(weeks)
colnames(cases) <- spat$Meshblock06

tab <- case %>% group_by(ReportWeek, Meshblock06) %>% summarize(Cases = n())
week_by_mb <- tab %>% spread(Meshblock06, Cases, fill = 0)
cases[as.character(week_by_mb$ReportWeek),names(week_by_mb)[-1]] <- as.matrix(week_by_mb[,-1])

# assemble data
time_map <- ifelse(weeks < "2008-01-01", 1, 2)
data <- check_data(list(cases=cases, popn=popn, mbrg=mbrg, nb=nb, t2p=time_map))
data$spat_list <- spat %>% select(ID = ID, Spatial = Meshblock06, Region = TA2006)
data$case_list <- case %>% select(CaseID = EpiSurvNumber, Spatial = Meshblock06, ReportWeek)

# fit the model
print(system.time({
mod <- fit_model(data, prior, control, seed = 1)
}))

# save model
new = list(data=data, mod=mod)
save(list = c('new'), file='TA2006.Rdata')
