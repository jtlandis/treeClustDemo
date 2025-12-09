#' @export
inner <- function(tree_vec) {
  trees <- lapply(
    trees(tree_vec),
    function(tree) {
      new_tree_vctr(
        node = inner_nodes(tree),
        which_tree = 1L,
        tree = list(tree),
        node_encoding = node_encoding(tree_vec)
      )
    }
  )
  vctrs::vec_c(!!!trees, .ptype = tree_vec)
}


#' creates a data.frame that maps inner_nodes
#' to the indices of the input tree_vctr. Note
#' that only unique values of tree_vctr are
#' considered. If leaf nodes are missing in the
#' input, then they will not appear in the output.
#' @export
generate_inner_slice <- function(tree_vec, node_level = NULL) {
  trees <- trees(tree_vec)
  nodes <- node_id(tree_vec)
  ind <- vctrs::vec_group_loc(tree_id(tree_vec))
  ind <- vctrs::vec_slice(ind, order(ind$key))
  inner <- lapply(trees, inner_nodes)
  table_fn <- switch(node_encoding(tree_vec),
    observation = function(nodes, tree) nodes,
    order = function(nodes, tree) {
      match(nodes, match(seq_len(length(tree$order)), tree$order))
    }
  )
  desc <- mapply(
    function(nodes, tree, ind, inner) {
      # desc are always encoded as observations
      desc <- node_descendants(tree)
      inner <- inner_nodes(tree)
      # convert nodes to observation encoding
      table <- data.frame(
        node = table_fn(nodes, tree),
        ind = ind
      )
      # match descendants obsv encoding to what we pulled
      lapply(desc[inner],
        \(desc, table) {
          vctrs::vec_slice(
            table,
            which(table$node %in% desc)
          )
        },
        table = table
      ) |>
        # slice the indices by this mapping.
        lapply(`[[`, "ind")
    },
    nodes = vctrs::vec_chop(nodes, ind$loc),
    tree = trees,
    ind = ind$loc,
    inner = inner,
    SIMPLIFY = FALSE
  )
  inner_sizes <- vapply(trees, \(tree) nrow(tree$merge), 1L)
  tbl <- tibble::tibble(
    nodes = new_tree_vctr(
      node = vctrs::list_unchop(inner),
      which_tree = rep(seq_along(trees), inner_sizes),
      tree = trees,
      node_encoding = node_encoding(tree_vec)
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


child_nodes <- function(tree, times = 1L, remove_leaf = TRUE) {
  merge_fn <- switch(node_encoding(tree),
    observation = function(merge, order) {
      merge
    },
    order = function(merge, order) {
      match(merge, order)
    }
  )
  new_nodes <- with_tree_vctr(
    tree,
    function(node, tree, ...) {
      merge <- tree$merge
      n <- length(tree$order)
      is_leaf <- merge < 0
      node_is_leaf <- function(node) {
        node <= n | is.na(node)
      }
      merge[is_leaf] <- merge_fn(-merge[is_leaf], tree$order)
      merge[!is_leaf] <- merge[!is_leaf] + n
      nodes <- as.list(node)
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
  )
  lns <- lengths(new_nodes)
  out <- new_pure_tree_vctr(
    node = vctrs::list_unchop(new_nodes),
    which_tree = rep(wt(tree), lns),
    tree = trees(tree),
    node_encoding = node_encoding(tree)
  )
  unique(out[!is.na(node_id(out))])
}


parent_nodes <- function(tree, times = 1L, remove_top = TRUE) {
  merge_fn <- switch(node_encoding(tree),
    observation = function(merge, order) {
      merge
    },
    order = function(merge, order) {
      match(merge, order)
    }
  )
  new_nodes <- with_tree_vctr(
    tree,
    function(node, tree, ...) {
      merge <- tree$merge
      n <- length(tree$order)
      n_inner <- n - 1L
      N <- n + n_inner
      is_leaf <- merge < 0
      node_is_gparent <- function(node) {
        node >= N | is.na(node)
      }
      merge[is_leaf] <- merge_fn(-merge[is_leaf], tree$order)
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
  )
  new_pure_tree_vctr(
    node = new_nodes,
    which_tree = wt(tree),
    tree = trees(tree),
    node_encoding = node_encoding(tree)
  )
}
