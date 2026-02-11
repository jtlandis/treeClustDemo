#' @export
root_nodes <- function(tree_vec) {
  fct <- wt(tree_vec)
  hashes <- levels(fct)
  all_trees <- trees(tree_vec)
  trees <- Map(
    function(tree, which_tree) {
      new_pure_tree_vctr(
        node = vctrs::vec_cast(tree_top_node(tree), integer()),
        which_tree = vctrs::new_factor(which_tree, levels = hashes),
        tree = all_trees
      )
    },
    tree = all_trees,
    which_tree = seq_along(all_trees)
  )

  vctrs::vec_c(!!!trees, .ptype = tree_vec)
}

#' @export
inner_nodes <- function(tree_vec) {
  fct <- wt(tree_vec)
  hashes <- levels(fct)
  all_trees <- trees(tree_vec)
  trees <- Map(
    function(tree, which_tree) {
      node <- vctrs::vec_cast(node_inner(tree), integer())
      new_pure_tree_vctr(
        node = node,
        which_tree = vctrs::new_factor(
          vctrs::vec_rep(which_tree, length(node)),
          levels = hashes
        ),
        tree = all_trees
      )
    },
    tree = all_trees,
    which_tree = seq_along(all_trees)
  )
  vctrs::vec_c(!!!trees, .ptype = tree_vec)
}

#' @export
leaf_nodes <- function(tree_vec) {
  fct <- wt(tree_vec)
  hashes <- levels(fct)
  all_trees <- trees(tree_vec)
  trees <- Map(
    function(tree, which_tree) {
      len <- tree_leaf_max_node(tree)
      node <- seq_len(len)
      new_pure_tree_vctr(
        node = node,
        which_tree = vctrs::new_factor(
          vctrs::vec_rep(which_tree, len),
          levels = hashes
        ),
        tree = all_trees
      )
    },
    tree = all_trees,
    which_tree = seq_along(all_trees)
  )
  vctrs::vec_c(!!!trees, .ptype = tree_vec)
}


#' @export
all_nodes <- function(tree_vec) {
  fct <- wt(tree_vec)
  hashes <- levels(fct)
  all_trees <- trees(tree_vec)
  trees <- Map(
    function(tree, which_tree) {
      len <- tree_leaf_max_node(tree)
      N <- tree_n_inner_node(tree)
      node <- seq_len(len + N)
      new_pure_tree_vctr(
        node = node,
        which_tree = vctrs::new_factor(
          vctrs::vec_rep(which_tree, len + N),
          levels = hashes
        ),
        tree = all_trees
      )
    },
    tree = all_trees,
    which_tree = seq_along(all_trees)
  )
  vctrs::vec_c(!!!trees, .ptype = tree_vec)
}


#' creates a data.frame that maps inner_nodes
#' to the indices of the input tree_vctr. Note
#' that only unique values of tree_vctr are
#' considered. If leaf nodes are missing in the
#' input, then they will not appear in the output.
#' @export
generate_inner_slice <- function(
  tree_vec, node_level = NULL,
  only_inner = TRUE
) {
  ## TO DO - check functionality
  ## after refactoring
  trees <- trees(tree_vec)
  nodes <- node_id(tree_vec)
  ind <- vctrs::vec_group_loc(tree_id(tree_vec))
  ind <- vctrs::vec_slice(ind, order(ind$key))
  table_fn <- function(nodes, tree) {
    match(nodes, match(seq_len(length(tree$order)), tree$order))
  }

  desc_id <- if (only_inner) {
    lapply(trees, node_inner)
  } else {
    lapply(trees, function(tree) {
      seq_along(node_descendants(tree))
    })
  }
  desc <- mapply(
    function(nodes, tree, ind, desc_id) {
      # desc are now encoded by other nodes
      desc <- node_descendants(tree)
      # convert nodes to observation encoding
      table <- data.frame(
        node = nodes,
        ind = ind
      )
      # match descendants obsv encoding to what we pulled
      lapply(desc[desc_id],
        \(desc, table) {
          vctrs::vec_slice(
            table$ind,
            which(table$node %in% desc)
          )
        },
        table = table
      )
    },
    nodes = vctrs::vec_chop(nodes, ind$loc),
    tree = trees,
    ind = ind$loc,
    desc_id = desc_id,
    SIMPLIFY = FALSE
  )
  sizes <- lengths(desc)
  tbl <- tibble::tibble(
    nodes = new_tree_vctr(
      node = vctrs::list_unchop(desc_id),
      which_tree = rep(seq_along(trees), sizes),
      tree = trees
    ),
    children = vctrs::list_unchop(desc)
  )
  if (!is.null(node_level)) {
    if (!is_tree_vctr(node_level) ||
      !any(m <- vctrs::vec_in(tbl$nodes, node_level))) {
      cli::cli_warn(
        c(
          "unable to locate inner slice node levels.",
          "`node_levels` should be a tree_vctr of inner nodes",
          "returning all possible levels"
        )
      )
      return(tbl)
    }
    tbl <- tbl[m, ]
  }
  tbl
}

node_child <- function(tree, nodes, times = 1L, remove_leaf = TRUE) {
  UseMethod("node_child", tree)
}

#' @noRd
#' @exportS3Method
node_child.hclust <- function(tree, nodes, times = 1L, remove_leaf = TRUE) {
  merge <- tree$merge
  n <- length(tree$order)
  is_leaf <- merge < 0
  node_is_leaf <- function(node) {
    node <= n | is.na(node)
  }
  # assuming order encoding
  merge[is_leaf] <- match(-merge[is_leaf], tree$order)
  merge[!is_leaf] <- merge[!is_leaf] + n
  nodes <- as.list(nodes)
  while (times > 0) {
    times <- times - 1L
    nodes <- lapply(nodes, function(node) {
      is_leaf <- node_is_leaf(node)
      nodes <- as.list(node)
      if (remove_leaf && any(is_leaf)) {
        nodes[is_leaf] <- NA_integer_
      }
      nodes[!is_leaf] <- vctrs::vec_chop(
        merge[node[!is_leaf] - n, , drop = FALSE]
      )
      vctrs::list_unchop(lapply(nodes, as.vector))
    })
  }
  nodes
}

#' @noRd
#' @exportS3Method
node_child.phylo <- function(tree, nodes, times = 1L, remove_leaf = TRUE) {
  edge <- tree$edge
  ansc <- edge[, 1]
  desc <- edge[, 2]
  data <- vctrs::vec_split(desc, ansc)
  n <- tree_leaf_max_node.phylo(tree)
  node_is_leaf <- function(node) {
    node <= n | is.na(node)
  }
  nodes <- as.list(nodes)
  while (times > 0) {
    times <- times - 1L
    nodes <- lapply(nodes, function(node) {
      is_leaf <- node_is_leaf(node)
      nodes <- as.list(node)
      if (remove_leaf && any(is_leaf)) {
        nodes[is_leaf] <- NA_integer_
      }
      nodes[!is_leaf] <- data$val[match(as.vector(nodes[!is_leaf]), data$key)]
      vctrs::vec_unchop(nodes)
    })
  }
  nodes
}

#' @export
child_nodes <- function(tree, times = 1L, remove_leaf = TRUE) {
  new_nodes <- with_tree_vctr(
    tree,
    function(node, tree, times, remove_leaf) {
      node_child(tree, nodes = node, times = times, remove_leaf = remove_leaf)
    },
    time = times,
    remove_leaf = remove_leaf
  )
  lns <- lengths(new_nodes)
  out <- new_pure_tree_vctr(
    node = vctrs::list_unchop(new_nodes),
    which_tree = rep(wt(tree), lns),
    tree = trees(tree)
  )
  unique(out[!is.na(node_id(out))])
}

node_parent <- function(tree, nodes, times = 1L, remove_top = TRUE) {
  UseMethod("node_parent")
}

#' @noRd
#' @exportS3Method
node_parent.hclust <- function(tree, nodes, times = 1L, remove_top = TRUE) {
  merge <- tree$merge
  n <- length(tree$order)
  n_inner <- n - 1L
  N <- n + n_inner
  is_leaf <- merge < 0
  node_is_gparent <- function(node) {
    node >= N | is.na(node)
  }
  merge[is_leaf] <- match(-merge[is_leaf], tree$order)
  merge[!is_leaf] <- merge[!is_leaf] + n
  while (times > 0) {
    times <- times - 1L
    is_gparent <- node_is_gparent(node)
    if (remove_top && any(is_gparent)) {
      node[is_gparent] <- NA_integer_
    }
    node[!is_gparent] <- (match(node[!is_gparent], merge) - 1L) %%
      n_inner + 1L + n
  }
  node
}

#' @noRd
#' @exportS3Method
node_parent.phylo <- function(tree, nodes, times = 1L, remove_top = TRUE) {
  edge <- tree$edge
  ansc <- edge[, 1]
  desc <- edge[, 2]
  data <- vctrs::vec_split(desc, ansc)
  N <- tree_top_node(tree)
  node_is_gparent <- function(node) {
    node == N | is.na(node)
  }
  while (times > 0) {
    times <- times - 1L
    is_gparent <- node_is_gparent(nodes)
    if (remove_top && any(is_gparent)) {
      nodes[is_gparent] <- NA_integer_
    }
    # get parent from desc....
    nodes[!is_gparent] <- ansc[match(nodes[!is_gparent], desc)]
  }
  nodes
}

#' @export
parent_nodes <- function(tree, times = 1L, remove_top = TRUE) {
  new_nodes <- with_tree_vctr(
    tree,
    function(node, tree, times, remove_top) {
      node_parent(tree, nodes = node, times = times, remove_top = remove_top)
    },
    time = times,
    remove_top = remove_top
  )
  new_pure_tree_vctr(
    node = vctrs::vec_cast(new_nodes, integer()),
    which_tree = wt(tree),
    tree = trees(tree)
  )
}
