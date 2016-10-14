
five_colour_graph <- function(A) {
  if (!(is.matrix(A) && is.logical(A))) {
    stop("Not an adjacency matrix\n")
  }

  # helper function to find the vertex of degree 5 with neighbours
  # of degree no more than 7 that are not adjacent, which must
  # exist for a planar graph
  find_degree_five <- function(A) {
    deg = rowSums(A)
    wch = which(deg == 5)
    for (i in seq_along(wch)) {
      nb = which(A[wch[i],])
      # find 2 neighbours N1, N2 such that they are not neighbours of each other, and such that each has
      # degree at most 7, and such that they're not neighbours of each other
      nb = nb[deg[nb] <= 7]
      if (length(nb) >= 2) {
        # run through and find the neighbours of each pair. There is likely a more
        # efficient way to do this!
        l = lapply(seq_along(nb), function(x) { setdiff(nb[-x], which(A[nb[x],])) })
        nb_wch = which(unlist(lapply(l, function(x) { length(x) > 0 })))
        if (length(nb_wch) > 0) {
          return(c(wch[i], nb[nb_wch[1]], l[[ nb_wch[1] ]][1]))
        }
      }
    }
    stop("Can't find a 5 degree vertex with non-adjacent neighbours of degree no more than 7. Non planer graph?\n")
  }

  deg = rowSums(A)
  cols = numeric(length(deg))

  if (length(deg) == 1) {
    # empty graph - we're done, use any colour we like
    cols[1] = 1
  } else if (any(deg < 5)) {
    # remove the vertex, and colour the previous
    v = which(deg < 5)[1]
    cols[-v] = five_colour_graph(A[-v,-v,drop=FALSE])
    col = setdiff(1:5, cols[A[v,]])
    if (length(col) == 0) {
      stop("Error colouring graph\n")
    }
    cols[v] = col[1]
  } else if (any (deg == 5)) {
    v = find_degree_five(A)
    # add new vertex N1/N2 linking to their neighbours
    linkto = A[v[1],] | A[v[2],]
    A_new = rbind(cbind(A, linkto), c(linkto, FALSE))

    # remove old vertices x, n1, n2 and colour
    col_sub = five_colour_graph(A_new[-v, -v, drop=FALSE])
    cols[-v] <- col_sub[-length(col_sub)]
    cols[v[2:3]] <- col_sub[length(col_sub)]
    col = setdiff(1:5, cols[A[v[1],]])
    if (length(col) == 0) {
      stop("Error colouring graph\n")
    }
    cols[v[1]] = col[1]
  } else {
    stop("Error colouring graph - not a planer graph?\n")
  }
  return(cols)
}

build_region_adjacency <- function(mbrg, nb) {
  # build the LUT
  rgmb = init_region_lut(mbrg)
  adj = matrix(FALSE, length(rgmb), length(rgmb))
  for (i in 1:length(rgmb)) {
    spat_in  = rgmb[[i]]
    spat_all = unique(as.numeric(nb[spat_in,-1]))
    spat_nb  = setdiff(spat_all, spat_in)
    rg_nb = unique(mbrg[spat_nb])
    adj[i, rg_nb] = TRUE
  }
  # check:
  if (!isTRUE(all.equal(adj, t(adj)))) {
    stop("Adjacency matrix isn't symmetric, so something has gone wrong\n")
  }
  adj
}

#' Colour regions on a map defined by groupings of spatial units using the five colour algorithm
#' 
#' @param mbrg The mapping from fine spatial units to regions
#' @param nb   The neighbourhood mapping of spatial units to each other
#' @return A vector containing the numbers 1 through 5 as to which colour each fine spatial unit
#' should be.
#' @export
five_colour_map <- function(mbrg, nb) {
  # build the adjacency matrix
  a <- build_region_adjacency(mbrg, nb)
  # colour the graph
  col_index = five_colour_graph(a)
  # return the coloured regions
  col_index[mbrg]
}
