wt <- function(x) {
  vctrs::field(x, "which_tree")
}

#' @export
tree_id <- function(x) {
  as.integer(wt(x))
}

#' @export
`tree_id<-` <- function(x, value) {
  fct <- wt(x)
  value <- vctrs::new_factor(value, levels = levels(fct))
  vctrs::field(x, "which_tree") <- value
  x
}

#' @export
node_id <- function(x) {
  vctrs::field(x, "node")
}

#' @export
`node_id<-` <- function(x, value) {
  vctrs::field(x, "node") <- value
  x
}

#' @export
trees <- function(x) {
  attr(x, "tree", exact = TRUE) %||% stop("malformed `tree_vctr` object.")
}

#' @export
tree <- function(x, which = NULL) {
  trees <- trees(x)
  if (is.null(which)) {
    if (length(trees) != 1L) {
      stop("`tree` expects there to be only one tree object on `tree_vctr`")
    }
    return(trees[[1L]])
  }
  which <- as.integer(which)
  trees[[which]]
}


node_descendants <- function(tree) {
  with_descendants(tree)$descendants
}

#' reports the descendants of each node according
#' to the current encoding of the tree_vctr. For
#' nodes that are leaves, this should be identical
#' to the node_id() funciton.
#' @export
descendants <- function(x) {
  with_tree_vctr(
    x,
    \(node, tree) node_descendants(tree)[node]
  )
}
