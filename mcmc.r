library(epiclustR)
library(dplyr)
library(tidyr)

# the priors and MCMC control parameters
prior <- init_priors()
control <- init_control(thinning=50, samples = 1000, burnin=20)

spat <- read.table("Meshblocks.txt", header=TRUE)

nb   <- load_spatial_neighbours("Weights.GAL")
popn <- spat$Population
mbrg <- spat$CAU

case <- read.csv("CaseData/cases.csv")

case <- case %>% filter(Meshblock06 %in% spat$Meshblock06)

# now compute the number of cases per meshblock per time
cases <- matrix(0, max(case$WeekFromStart), length(mbrg))
rownames(cases) <- 1:max(case$WeekFromStart)
colnames(cases) <- spat$Meshblock06

tab <- case %>% group_by(WeekFromStart, Meshblock06) %>% summarize(Cases = n())
week_by_mb <- tab %>% spread(Meshblock06, Cases, fill = 0)
cases[as.character(week_by_mb$WeekFromStart),names(week_by_mb)[-1]] <- as.matrix(week_by_mb[,-1])

# assemble data
data <- check_data(list(cases=cases, popn=popn, mbrg=mbrg, nb=nb))

# fit the model
print(system.time({
posterior <- fit_model(data, prior, control, seed = 1)
}))

# do analysis
ssapply <- function(x, fun, ...) {
  simplify2array(lapply(x, function(xx, fun, ...) { simplify2array(lapply(xx, fun, ...)) }, fun = fun, ...))
}

extract_variable <- function(x, variable) {
  x[[variable]]
}

# now dump out our posterior chains
R = ssapply(posterior, extract_variable, 'R')
U = ssapply(posterior, extract_variable, 'U')
betaX = ssapply(posterior, extract_variable, 'betaX')
fe = ssapply(posterior, extract_variable, 'fe')
kR = ssapply(posterior, extract_variable, 'kR')
kU = ssapply(posterior, extract_variable, 'kU')
X = ssapply(posterior, extract_variable, 'X')
pX = ssapply(posterior, extract_variable, 'pX')
ecases = ssapply(posterior, cases_per_time, data=data, smoothed=FALSE)
scases = ssapply(posterior, cases_per_time, data=data, smoothed=TRUE)
plot(apply(ecases, 1, median), type='l')
lines(apply(scases, 1, median), col='red')
#hist(apply(R, 1, ess))

library(ggplot2)
pX %>% data.frame %>%
  mutate(iteration=seq_len(n())) %>%
  gather('chain', 'value', -iteration) %>%
  ggplot(aes(x=iteration, y=value)) + geom_line(aes(colour=chain))

#write.table(t(R), 'MidCentral/current_version/R.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(U), 'MidCentral/current_version/U.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(betaX), 'MidCentral/current_version/betaX.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(fe), 'MidCentral/current_version/fixedEffects.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(ecases), 'MidCentral/current_version/expectedCases.txt', row.names=FALSE, col.names=FALSE)
#write.table(t(scases), 'MidCentral/current_version/smoothedCases.txt', row.names=FALSE, col.names=FALSE)

compare_var <- function(variable) {
  ref <- scan(file.path('MidCentral/RUX2_region2', paste0(variable, '.txt')))
  new <- scan(file.path('MidCentral/current_version', paste0(variable, '.txt')))
  print(all.equal(ref, new))
}

# compare versions
compare_var('fixedEffects')
compare_var('R')
compare_var('U')
compare_var('betaX')
compare_var('expectedCases')
compare_var('smoothedCases')
