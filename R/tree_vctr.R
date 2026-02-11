#' a tree_vctr object
#'
#' This object allows a way to encode many heirachical objects along with
#' the rows of a matrix-like object. Each element of a `tree_vctr` object
#' will have a node id and a tree id, which together identify a specific
#' observation within the `hclust` object stored in the `tree` attribute.
#' DEFUNCT: The node id may be encoded in two ways, which mainly provide convience
#' with sorting this vector according to either the order of observations
#' or an order that would cluster the observations per hclust object.
#'
#' @param node an integer vector represeting some encoded node id
#' @param which_tree identifies which element of `node` maps to which hclust
#' tree in the `tree` attribute
#' @param tree a list of `hclust` objects
#' @export
new_tree_vctr <- function(
  node = integer(),
  which_tree = integer(),
  tree
) {
  node <- vctrs::vec_cast(node, integer())
  which_tree <- vctrs::vec_recycle(
    vctrs::vec_cast(which_tree, integer()),
    size = vctrs::vec_size(node)
  )
  if (is_tree_obj(tree)) {
    tree <- list(tree)
  }

  stopifnot(
    "`tree` should be a list or an hclust/phylo object" = is.list(tree),
    "all elements of `tree` should be an hclust/phylo object" = all(vapply(
      tree,
      is_tree_obj,
      FUN.VALUE = FALSE
    ))
  )
  tree <- lapply(tree, with_descendants)
  hashes <- vctrs::vec_unique(vapply(tree, digest::digest, ""))
  new_pure_tree_vctr(
    node = node,
    which_tree = vctrs::new_factor(which_tree, levels = hashes),
    tree = tree
  )
}

tree_vctr_class <- c("tree_vctr", "vctrs_rcrd", "vctrs_vctr")

is_tree_obj <- function(obj) {
  inherits(obj, "hclust") || inherits(obj, "phylo")
}

new_pure_tree_vctr <- function(
  node = integer(),
  which_tree = factor(),
  tree = list()
) {
  out <- list(
    which_tree = which_tree,
    node = node
  )
  # attributes(out) <- list(
  #   names = c("which_tree", "node"),
  #   tree = tree, class = tree_vctr_class
  # )
  attr(out, "tree") <- tree
  class(out) <- tree_vctr_class
  out
}


#' @export
with_tree_vctr <- function(vec, .f, ...) {
  trees <- trees(vec)
  nodes <- node_id(vec)
  if (length(trees) == 1L) {
    .f(node = nodes, tree = trees[[1L]], ...)
  } else {
    tree_ind <- vctrs::vec_group_loc(tree_id(vec))$loc
    out <- vector("list", length(trees))
    for (i in seq_along(trees)) {
      tree <- trees[[i]]
      ind <- tree_ind[[i]]
      out[[i]] <- .f(node = nodes[ind], tree = tree, ...)
    }
    vctrs::list_unchop(out, indices = tree_ind)
  }
}


#' @export
format.tree_vctr <- function(x, ...) {
  which_tree <- vctrs::field(x, "which_tree")
  out <- sprintf("%i:%i", as.integer(which_tree), vctrs::field(x, "node"))
  out[is.na(which_tree)] <- "<NA>"
  out
}

#' @export
obj_print_footer.tree_vctr <- function(x, ...) {
  n_trees <- length(trees(x))
  cat(
    sep = "",
    "# n tree(s): ", n_trees, "\n"
  )
}


tree_vctr_drop <- function(x) {
  wt <- wt(x)
  old_lvls <- levels(wt)
  # 4, 2, 3, 4, 2
  ints <- as.integer(wt)
  # 4, 2, 3
  unique_tree_id <- vctrs::vec_unique(ints)
  # 2, 3, 1
  unique_ord <- vctrs::vec_order(unique_tree_id)
  # 2, 3, 4
  unique_tree_id <- vctrs::vec_slice(
    unique_tree_id,
    unique_ord[!is.na(unique_ord)]
  )
  to_trees <- trees(x)
  if (!all(is.na(unique_tree_id))) {
    to_trees <- vctrs::vec_slice(trees(x), unique_tree_id)
    old_lvls <- vctrs::vec_slice(old_lvls, unique_tree_id)
  }
  new_pure_tree_vctr(
    node = node_id(x),
    # 3, 1, 2, 3, 1
    which_tree = vctrs::new_factor(
      vctrs::vec_match(ints, unique_tree_id),
      old_lvls
    ),
    tree = to_trees
  )
}

#' @export
tree_vctr <- function(x) {
  UseMethod("tree_vctr", x)
}

#' @export
tree_vctr.tree_vctr <- function(x) {
  tree_vctr_drop(x)
}

#' @export
tree_vctr.hclust <- function(x) {
  as_tree_vctr(x)
}

#' @export
tree_vctr.phylo <- function(x) {
  desc <- x$edge[, 2]
  new_tree_vctr(
    node = desc[desc <= tree_leaf_max_node.phylo(x)],
    which_tree = 1L,
    tree = x
  )
}

#' @export
tree_vctr.default <- function(x) {
  tree_vctr(as.hclust(x))
}


#' @importFrom vctrs vec_ptype2
#' @export
vec_ptype2.tree_vctr.tree_vctr <- function(x, y, ...) {
  xwt <- levels(vctrs::field(x, "which_tree"))
  ywt <- levels(vctrs::field(y, "which_tree"))
  ykeep <- which(!vctrs::vec_in(ywt, xwt))
  xtree <- attr(x, "tree")
  if (length(ykeep)) {
    xwt <- vctrs::vec_c(xwt,
      vctrs::vec_slice(ywt, ykeep),
      .ptype = character()
    )
    xtree <- vctrs::vec_c(xtree,
      vctrs::vec_slice(attr(y, "tree"), ykeep),
      .ptype = list()
    )
  }
  # opting for speed here
  which_tree <- integer()
  attributes(which_tree) <- list(
    levels = xwt,
    class = "factor"
  )
  new_pure_tree_vctr(
    which_tree = which_tree,
    tree = xtree
  )
}

#' @importFrom vctrs vec_cast
#' @export
vec_cast.tree_vctr.tree_vctr <- function(x, to, ...) {
  vctrs::field(x, "which_tree") <- vctrs::allow_lossy_cast(
    vctrs::vec_cast(
      vctrs::field(x, "which_tree"),
      vctrs::field(to, "which_tree")
    )
  )
  x
}


tree_labels <- function(x) {
  with_tree_vctr(x, function(node, tree) {
    n <- length(tree$order)
    labels <- sprintf("node%i", node)
    is_leaf <- node <= n
    tree_labels <- tree$labels
    if (!is.null(tree_labels)) {
      tree_labels <- tree_labels[tree$order]
      labels[is_leaf] <- tree_labels[node][is_leaf]
    }
    labels
  })
}

leaf_labels <- function(tree, use_labels = TRUE) {
  vctrs::vec_unchop(
    mapply(
      \(tree, i) {
        n <- length(tree$order)
        labels <- tree$labels[tree$order]
        if (is.null(labels) || !use_labels) {
          labels <- sprintf("node%i", seq_len(n))
        }
        sprintf("%s-%i", labels, i)
      },
      tree = trees(tree),
      i = seq_along(trees(tree)),
      SIMPLIFY = FALSE
    )
  )
}

#' @export
as.character.tree_vctr <- function(x) {
  sprintf("%s-%i", tree_labels(x), tree_id(x))
}

#' @export
as.factor.tree_vctr <- function(x) {
  factor(as.character(x), levels = leaf_labels(x))
}

#' @export
is.infinite.tree_vctr <- function(x) {
  logical(vctrs::vec_size(x))
}

#' @export
vec_ptype2.tree_vctr.double <- function(x, y, ...) {
  double()
}

#' @export
vec_ptype2.double.tree_vctr <- function(x, y, ...) {
  double()
}

#' @export
vec_cast.double.tree_vctr <- function(x, to, ...) {
  as.double(as.factor(x))
}

#' @export
vec_cast.tree_vctr.integer <- function(x, to, ...) {
  mtch <- leaf_labels(to, use_labels = FALSE)[x]
  mtch <- sub("(node|obsv)", "", x = mtch)
  split <- strsplit(mtch, "-")
  nodes <- vapply(split, `[[`, FUN.VALUE = "", 1L) |> as.integer()
  if (any(is_na <- is.na(mtch))) {
    wt <- nodes
    wt[!is_na] <- vapply(split[!is_na], `[[`, FUN.VALUE = "", 2L) |>
      as.integer()
  } else {
    wt <- vapply(split, `[[`, FUN.VALUE = "", 2L) |>
      as.integer()
  }

  new_pure_tree_vctr(
    node = nodes,
    which_tree = vctrs::new_factor(wt, levels(wt(to))),
    tree = trees(to)
  )
}

#' @export
vec_cast.tree_vctr.double <- function(x, to, ...) {
  vec_cast(
    vec_cast(x, integer()),
    to,
    ...
  )
}

#' @export
vec_ptype2.tree_vctr.integer <- function(x, y, ...) {
  integer()
}

#' @export
vec_ptype2.integer.tree_vctr <- function(x, y, ...) {
  integer()
}

#' @export
vec_cast.integer.tree_vctr <- function(x, to, ...) {
  as.integer(as.factor(x))
}
