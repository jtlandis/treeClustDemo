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
tree <- function(x) {
  trees <- trees(x)
  if (length(trees) != 1L) {
    stop("`tree` expects there to be only one tree object on `tree_vctr`")
  }
  trees[[1L]]
}


#' @export
node_encoding <- function(x) {
  out <- attr(x, "node_encoding", exact = TRUE) %||%
    stop("malformed `tree_vctr` object.")
  if (!out %in% c("order", "observation")) {
    stop("malformed `tree_vctr` object.")
  }
  out
}




swap_node_encoding <- function(x) {
  node_encoding <- node_encoding(x)
  fn <- switch(node_encoding,
    observation = function(node, tree) {
      n <- length(tree$order)
      is_leaf <- node <= n
      swap <- match(node, tree$order)
      node[is_leaf] <- swap[is_leaf]
      node
    },
    order = function(node, tree) {
      n <- length(tree$order)
      is_leaf <- node <= n
      swap <- match(
        node,
        match(seq_len(length(tree$order)), tree$order)
      )
      node[is_leaf] <- swap[is_leaf]
      node
    }
  )
  new_nodes <- with_tree_vctr(x, fn)
  new_tree_vctr(
    node = new_nodes,
    which_tree = tree_id(x),
    tree = trees(x),
    node_encoding = switch(node_encoding,
      observation = "order",
      order = "observation"
    )
  )
}

#' @export
encode_order <- function(x) {
  switch(node_encoding(x),
    observation = swap_node_encoding(x),
    order = x
  )
}

#' @export
encode_obsv <- function(x) {
  switch(node_encoding(x),
    observation = x,
    order = swap_node_encoding(x)
  )
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
  switch(node_encoding(x),
    observation = with_tree_vctr(
      x,
      function(node, tree) {
        node_descendants(tree)[node]
      }
    ),
    order = with_tree_vctr(
      x,
      function(node, tree) {
        desc <- node_descendants(tree)
        inner_seq <- inner_nodes(tree)
        desc[inner_seq] <- lapply(
          desc[inner_seq], match,
          table = tree$order
        )
        desc[node]
      }
    )
  )
}
