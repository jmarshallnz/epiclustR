library(epiclustR)

# the priors and MCMC control parameters
prior <- init_priors()
control <- init_control(thinning=50, samples = 1000, burnin=20)

# load spatial data
load_spatial_neighbours <- function(path, file = "Weights.GAL") {
  input<-scan(file.path(path, file))
  n <- input[2]
  spatial_list <- list()
  i<-3
  for (j in 1:n) {
    k <- input[i] # name of spatial location
    nb <- input[i+1+1:input[i+1]] # it's neighbours
    i <- i + 2 + length(nb)
    spatial_list[[k]] <- nb
  }
  if (n != length(spatial_list))
    stop("Invalid spatial neighbours: Incorrect length")
  length_neighbours <- unlist(lapply(spatial_list, length))
  max_neighbours <- max(length_neighbours)
  spatial_matrix <- matrix(0, n, max_neighbours + 1)
  spatial_matrix[,1] <- length_neighbours
  for (j in 1:n)
    spatial_matrix[j,1 + seq_along(spatial_list[[j]])] <- spatial_list[[j]]
  spatial_matrix
}

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

init_region_lut <- function(mbrg) {
  lut <- list()
  rgs <- max(mbrg)
  for (j in 1:rgs) {
    lut[[j]] <- which(mbrg==j)
  }
  lut
}

nb   <- load_spatial_neighbours("MidCentral", "Weights.GAL")
popn <- load_spatial("MidCentral", "Meshblocks.txt")
mbrg <- load_regions("MidCentral", "Regions2.txt")
rgmb <- init_region_lut(mbrg)

# load the cases
load_cases <- function(path, file) {
  input<-scan(file.path(path, file))
  tps <- 312
  matrix(input, tps)
}
cases <- load_cases("MidCentral", "Data.txt")

# sanity checking...
check_popn <- function(n, cases) {
  wch <- which(n == 0 & apply(cases, 2, sum) > 0)
  if (length(wch > 0)) {
    cat("Setting population to 1 in spatial locations ", wch, " as we have cases there and popn=0\n")
  }
  n[wch] = 1
  n
}

popn <- check_popn(popn, cases)

# construct the state
set.seed(1)
pX = 0.1
state <- list(fe = -10,
              R  = rnorm(nrow(cases),0,1),
              U = rnorm(length(popn),0,1),
              X = matrix(rbinom(nrow(cases)*length(rgmb),1,pX),nrow(cases),length(rgmb)),
              betaX = rep(0.2,length(rgmb)),
              pX = 0.1,
              kU = 1,
              kR = 1,
              acceptR = rep(0, 1+4), # TODO: Move these somewhere else...
              rejectR = rep(0, 1+4),
              acceptU = rep(0,2),
              rejectU = rep(0,2),
              acceptX = 0,
              rejectX = 0)

# data
data <- list(cases=cases, popn=popn, mbrg=mbrg, nb=nb, rgmb=rgmb)

# fit the model
print(system.time({
posterior <- fit_model(data, state, prior, control)
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
