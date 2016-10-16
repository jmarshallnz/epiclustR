#' Computes the scaled eigen-vectors of a matrix A
#'
#' @param A the matrix to compute scaled eigen-vectors on
#' @return the scaled eigenvectors of A
eigen_scale <- function(A) {
  eig <- eigen(A, symmetric=TRUE)
  eig$vectors %*% diag(sqrt(pmax(eig$values,0)))
}

#' block_precision_order_2
#' 
#' Computes the precision matrix for a random walk of order 2
#' of the given block size
#'
#' @param n the block size
#' @return The precision matrix K corresponding to a random walk of order 2.
block_precision_order_2 <- function(n) {
  # Generates the structure matrix for a random walk of order 2
  K<-6*diag(2 + n + 2) # need 2 on either side for ends
  for (j in 1:nrow(K)) {
    if (j>2) {K[j,j-2]<-1}
    if (j>1) {K[j,j-1]<--4}
    if (j+1<=ncol(K)) {K[j,j+1]<--4}
    if (j+2<=ncol(K)) {K[j,j+2]<-1}
  }
  K[1,1:2] <- K[n+4,n+4:3] <- c(1, -2)
  K[2,1:2] <- K[n+3,n+4:3] <- c(-2, 5)
  K
}

#' construct_block_matrices
#'
#' Construct the cholesky/eigen matrix for a block update of the given length as well
#' as pre and post matrices (left and right of main block)
#'
#' @param block_length the length of the block update
#' @return a list of constructed block matrices containing K_eigen, K_bef and K_aft
construct_block_matrices <- function(block_length) {
  K_eigen <- list()
  K_before <- list()
  K_after <- list()
  for (o in 0:4) { # offset
    co <- o+1
    K <- block_precision_order_2(block_length)
    invK <- solve(K[1:block_length+o, 1:block_length+o])
    K_eigen[[co]] <- eigen_scale(invK)
    Kbef <- K[1:block_length+o,seq_len(min(o,2))+max(o,2)-2, drop=FALSE]
    Kaft <- K[1:block_length+o,o + block_length + seq_len(min(4-o,2)), drop=FALSE]
    K_before[[co]] <- -invK %*% Kbef
    K_after[[co]] <-  -invK %*% Kaft
  }
  list(K_eigen = K_eigen,
       K_before = K_before,
       K_after = K_after)
}
