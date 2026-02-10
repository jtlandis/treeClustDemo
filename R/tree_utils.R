as_tree_vctr <- function(hc) {
  new_tree_vctr(
    node = order(hc$order),
    which_tree = 1L,
    tree = list(hc)
  )
}

with_descendants <- function(tree) {
  if (is.null(tree$descendants)) {
    m <- tree$merge
    n <- nrow(m) + 1L
    is_merge <- m < 0
    # descendants are indexed by observation
    m[is_merge] <- -m[is_merge]
    m[!is_merge] <- m[!is_merge] + n
    node_end <- n + nrow(m)
    children <- vector("list", length = node_end)
    children[1:n] <- as.list(seq_len(n))
    for (i in seq_len(n - 1)) {
      ind <- m[i, ]
      children[[i + n]] <- unlist(children[ind])
    }
    tree$descendants <- children
  }
  tree
}


inner_nodes <- function(tree) {
  n <- length(tree$order)
  (n + 1L):(n + nrow(tree$merge))
}

#' @export
is_inner <- function(x) {
  with_tree_vctr(x, function(node, tree) {
    node > length(tree$order)
  })
}

#' @export
is_tree_vctr <- function(x) {
  inherits(x, "tree_vctr")
}


phylo_to_hclust <- function(x) {
  n <- length(x$tip.label)
  if (n == 1) {
    stop("needs n >= 2 observations for a classification")
  }
  is_tip <- x$edge[, 2] <= n
  order <- x$edge[is_tip, 2]
  if (n == 2) {
    m <- matrix(c(-1L, -2L), 1, 2)
    bt <- x$edge.length[1]
  } else {
    x$node.label <- NULL
    bt <- ape::branching.times(x)
    N <- n - 1L
    x <- reorder(x, "postorder")
    m <- matrix(x$edge[, 2], N, 2, byrow = TRUE)
    anc <- x$edge[c(TRUE, FALSE), 1]
    bt <- bt[as.character(anc)]
    bt <- c(sort(bt[-N]), bt[N])
    o <- match(names(bt), anc)
    m <- m[o, ]
    TIPS <- m <= n
    m[TIPS] <- -m[TIPS]
    oldnodes <- as.numeric(names(bt))[-N]
    m[match(oldnodes, m)] <- 1:(N - 1)
    names(bt) <- NULL
  }
  obj <- list(
    merge = m, height = 2 * bt, order = order, labels = x$tip.label,
    call = match.call(), method = "unknown"
  )
  class(obj) <- "hclust"
  obj
}
