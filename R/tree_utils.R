as_tree_vctr <- function(hc) {
  new_tree_vctr(
    node = order(hc$order),
    which_tree = 1L,
    tree = list(hc)
  )
}

#' utility function
#' @param tree a hclust/phylo obj
#' @noRd
with_descendants <- function(tree) {
  if (is.null(tree$descendants)) {
    tree <- UseMethod("with_descendants")
  }
  tree
}

#' @noRd
#' @exportS3Method
with_descendants.hclust <- function(tree) {
  m <- tree$merge
  n <- nrow(m) + 1L
  is_merge <- m < 0
  # descendants are indexed by node...
  m[is_merge] <- match(-m[is_merge], tree$order)
  m[!is_merge] <- m[!is_merge] + n
  node_end <- n + nrow(m)
  children <- vector("list", length = node_end)
  children[1:n] <- as.list(seq_len(n))
  for (i in seq_len(n - 1)) {
    ind <- m[i, ]
    children[[i + n]] <- unlist(children[ind])
  }
  tree$descendants <- children
  tree
}

#' @noRd
#' @exportS3Method
with_descendants.phylo <- function(tree) {
  e <- tree$edge
  ansc <- e[, 1]
  desc <- e[, 2]
  n <- min(ansc) - 1L
  node_end <- n + tree$Nnode
  children <- vector("list", length = node_end)
  seq_n <- seq_len(n)
  known <- logical(node_end)
  children[seq_n] <- as.list(seq_n)
  known[seq_n] <- TRUE
  ans_ind <- vctrs::vec_group_loc(ansc)
  ans_i <- nrow(ans_ind)
  while (ans_i > 0) {
    flag_rm <- logical(ans_i)
    for (i in seq_len(ans_i)) {
      target <- ans_ind$key[i]
      these_desc <- desc[ans_ind$loc[[i]]]
      if (any(!known[these_desc])) {
        next
      }
      flag_rm[i] <- TRUE
      known[target] <- TRUE
      children[[target]] <- unlist(children[these_desc])
    }
    ans_ind <- ans_ind[!flag_rm, ]
    ans_i <- nrow(ans_ind)
  }
  tree$descendants <- children
  tree
}

node_inner <- function(tree) {
  UseMethod("node_inner", tree)
}

#' @noRd
#' @exportS3Method
node_inner.hclust <- function(tree) {
  n <- length(tree$order)
  (n + 1L):(n + nrow(tree$merge))
}

#' @noRd
#' @exportS3Method
node_inner.phylo <- function(tree) {
  ansc <- tree$edge[, 1]
  min(ansc):max(ansc)
}

tree_top_node <- function(tree) {
  which.max(lengths(node_descendants(tree)))
}

tree_leaf_max_node <- function(tree) {
  UseMethod("tree_leaf_max_node")
}

#' @noRd
#' @exportS3Method
tree_leaf_max_node.hclust <- function(tree) {
  length(tree$order)
}

#' @noRd
#' @exportS3Method
tree_leaf_max_node.phylo <- function(tree) {
  min(tree$edge[, 1]) - 1L
}

tree_n_inner_node <- function(tree) {
  UseMethod("tree_n_inner_node")
}

#' @noRd
#' @exportS3Method
tree_n_inner_node.hclust <- function(tree) {
  nrow(tree$merge)
}

#' @noRd
#' @exportS3Method
tree_n_inner_node.phylo <- function(tree) {
  tree$Nnode
}


#' @export
is_inner <- function(x) {
  with_tree_vctr(x, function(node, tree) {
    node > tree_leaf_max_node(tree)
  })
}

#' @export
is_leaf <- function(x) {
  !is_inner(x)
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
    #
    labels <- x$node.label %||% sprintf("node-%i", n + seq_len(x$Nnode))
    x$node.label <- NULL
    bt <- ape::branching.times(x)
    .nms <- names(bt)
    bt <- data.frame(bt = bt, labs = labels)
    rownames(bt) <- .nms
    N <- n - 1L
    x <- reorder(x, "postorder")
    m <- matrix(x$edge[, 2], N, 2, byrow = TRUE)
    anc <- x$edge[c(TRUE, FALSE), 1]
    bt <- bt[as.character(anc), ]
    bt_last <- bt[N, ]
    bt <- bt[-N, ]
    bt <- rbind(bt[order(bt$bt, decreasing = FALSE), ], bt_last)
    o <- match(row.names(bt), anc)
    m <- m[o, ]
    if (!all(is.na(bt$labs))) {
      row.names(m) <- bt$labs
    }
    TIPS <- m <= n
    m[TIPS] <- -m[TIPS]
    .nms <- row.names(bt)
    bt <- bt$bt
    names(bt) <- .nms
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
