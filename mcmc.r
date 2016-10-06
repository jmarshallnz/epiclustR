library(epiclustR)

# the priors and MCMC control parameters
prior <- init_priors()
control <- init_control(thinning=50, samples = 1000, burnin=20)

load_spatial <- function(path, file = "Meshblocks.txt") {
  input<-read.table(file.path(path, file), header=FALSE)
  n <- input[,2]
  return(n)
}

load_regions <- function(path, file) {
  regions <- read.table(file.path(path, file), header=FALSE)
  mbrg <- regions[,1]
  return(mbrg)
}

load_cases <- function(path, file) {
  input<-scan(file.path(path, file))
  tps <- 312
  matrix(input, tps)
}

nb   <- load_spatial_neighbours(file.path("MidCentral", "Weights.GAL"))
popn <- load_spatial("MidCentral", "Meshblocks.txt")
mbrg <- load_regions("MidCentral", "Regions2.txt")
cases <- load_cases("MidCentral", "Data.txt")

# assemble data
data <- list(cases=cases, popn=popn, mbrg=mbrg, nb=nb)

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
