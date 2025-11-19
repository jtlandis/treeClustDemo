#' @export
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
