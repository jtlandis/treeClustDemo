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
generate_inner_slice <- function(tree_vec) {
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
      table <- table_fn(nodes, tree)
      # match descendants obsv encoding to what we pulled
      lapply(desc[inner], match, table = table, nomatch = 0L) |>
        # slice the indices by this mapping.
        lapply(vctrs::vec_slice, x = ind)
    },
    nodes = vctrs::vec_chop(nodes, ind$loc),
    tree = trees,
    ind = ind$loc,
    inner = inner,
    SIMPLIFY = FALSE
  )
  inner_sizes <- vapply(trees, \(tree) nrow(tree$merge), 1L)
  tibble::tibble(
    nodes = new_tree_vctr(
      node = vctrs::list_unchop(inner),
      which_tree = rep(seq_along(trees), inner_sizes),
      tree = trees,
      node_encoding = node_encoding(tree_vec)
    ),
    children = vctrs::list_unchop(desc)
  )
}
