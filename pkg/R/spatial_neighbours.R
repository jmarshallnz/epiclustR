#' Load spatial neighbours data
#' 
#' Loads a spatial neighbours file, of the format
#' 
#' 0 6
#' 1 5
#' 2 3 4 5 6
#' 2 2 
#' 1 3
#' 
#' where even rows indicate the spatial unit being considered, and
#' it's number of neighbours, and odd rows indicating what the neighbours are.
#' i.e. there are 6 spatial units in the
#' file above, and unit 1 has 5 neighbours which are units 2 through 5,
#' while unit 2 has 2 neighbours (1 and 3).
#' 
#' @param file the file to load.
#' @return A matrix of spatial neighbours, the first column of which is the number of neighbours
#' @export
load_spatial_neighbours <- function(file) {
  input <- scan(file)
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

#' Initialize a reverse lookup table from a lookup vector
#' 
#' @param luv A factor mapping items to a grouping variable
#' @return a list for reverse lookups from levels to vector entries in `luv`
init_lut <- function(luv) {
  lut <- list()
  luv <- factor(luv)
  n <- levels(luv)
  for (j in n) {
    lut[[j]] <- which(luv==j)
  }
  lut
}
