#' @importFrom igraph %>% neighbors degree delete_vertices add_vertices add_edges graph_from_adjacency_matrix
NULL

find_degree_five <- function(g) {
  wch = which(igraph::degree(g) == 5)
  for (i in seq_along(wch)) {
    nb = igraph::neighbors(g, wch[i])
    # find 2 neighbours N1, N2 such that they are not neighbours of each other, and such that each has
    # degree at most 7, and such that they're not neighbours of each other
    nb = nb[igraph::degree(g, nb) <= 7]
    if (length(nb) >= 2) {
      # run through and find the neighbours of each pair
      l = lapply(seq_along(nb), function(x) { setdiff(nb[-x], igraph::neighbors(g, nb[x])) })
      nb_wch = which(unlist(lapply(l, function(x) { length(x) > 0 })))
      if (length(nb_wch) > 0) {
        return(list(x = wch[i], n1 = nb[nb_wch[1]], n2 = V(g)[l[[nb_wch[1]]][1]]))
      }
    }
  }
  return(NULL)
}

five_colour_graph <- function(g) {
  d = igraph::degree(g)
  cols = numeric(length(d))

  if (length(d) == 1) {
    # empty graph - we're done, use any colour we like
    cols[1] = 1
  } else if (any(d < 5)) {
    # remove the vertex, and colour the previous
    v = which(d < 5)[1]
    cols[-v] = five_colour_graph(igraph::delete_vertices(g, v))
    nb_cols = cols[igraph::neighbors(g, v)]
    col = setdiff(1:5, nb_cols)
    if (length(col) == 0) {
      stop("Error colouring graph\n")
    }
    cols[v] = col[1]
  } else if (any (d == 5)) {
    xn = find_degree_five(g)
    if (is.null(xn)) {
      stop("Can't find a 5 degree vertex with non-adjacent neighbours of degree no more than 7. Non planer graph?\n")
    }
    # now we remove node x, n1, n2 from the graph, and add a new node N1/N2 that
    # links previous N1,N2 to the others
    linkto = union(igraph::neighbors(g, xn$n1), igraph::neighbors(g, xn$n2))
    el = rep(length(V(g))+1, 2*length(linkto))
    el[seq(1,length(el), by=2)] = linkto
    g_del = igraph::add_vertices(g, 1) %>%
      igraph::add_edges(el) %>%
      igraph::delete_vertices(c(xn$x, xn$n1, xn$n2))
    col_sub = five_colour_graph(g_del)
    cols[-c(xn$x, xn$n1, xn$n2)] <- col_sub[-length(col_sub)]
    cols[c(xn$n1, xn$n2)] <- col_sub[length(col_sub)]
    nb_cols = cols[igraph::neighbors(g, xn$x)]
    col = setdiff(1:5, nb_cols)
    if (length(col) == 0) {
      stop("Error colouring graph\n")
    }
    cols[xn$x] = col[1]
  } else {
    stop("Error colouring graph - not a planer graph?\n")
  }
  return(cols)
}

build_region_adjacency <- function(mbrg, nb) {
  # build the LUT
  rgmb = init_region_lut(mbrg)
  adj = matrix(0, length(rgmb), length(rgmb))
  for (i in 1:length(rgmb)) {
    spat_in  = rgmb[[i]]
    spat_all = unique(as.numeric(nb[spat_in,-1]))
    spat_nb  = setdiff(spat_all, spat_in)
    rg_nb = unique(mbrg[spat_nb])
    adj[i, rg_nb] = 1
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
  # convert to a graph
  g <- igraph::graph_from_adjacency_matrix(a, mode='undirected')
  # colour the graph
  col_index = five_colour_graph(g)
  # return the coloured regions
  col_index[mbrg]
}
