#' a tree_vctr object
#'
#' This object allows a way to encode many heirachical objects along with
#' the rows of a matrix-like object. Each element of a `tree_vctr` object
#' will have a node id and a tree id, which together identify a specific
#' observation within the `hclust` object stored in the `tree` attribute.
#' The node id may be encoded in two ways, which mainly provide convience
#' with sorting this vector according to either the order of observations
#' or an order that would cluster the observations per hclust object.
#'
#' @param node an integer vector represeting some encoded node id
#' @param which_tree identifies which element of `node` maps to which hclust
#' tree in the `tree` attribute
#' @param tree a list of `hclust` objects
#' @param node_encoding character scalar of how this object's node is encoded.
#' The default, "order", will encode nodes according to how observations are
#' ordered in the respective `hclust` object. The alternative, "observation",
#' will encode nodes according to the original observation indices used to
#' generate the hclust object. These options only affect leaf nodes, as
#' internal nodes will always be encoded by the height of the hclust.
#' @export
new_tree_vctr <- function(
  node = integer(),
  which_tree = integer(),
  tree,
  node_encoding = c("order", "observation")
) {
  node <- vctrs::vec_cast(node, integer())
  which_tree <- vctrs::vec_recycle(
    vctrs::vec_cast(which_tree, integer()),
    size = vctrs::vec_size(node)
  )
  if (inherits(tree, "hclust")) {
    tree <- list(tree)
  }

  stopifnot(
    "`tree` should be a list or an hclust object" = is.list(tree),
    "all elements of `tree` should be an hclust object" = all(vapply(
      tree,
      inherits,
      what = "hclust",
      FUN.VALUE = FALSE
    ))
  )
  tree <- lapply(tree, with_descendants)
  hashes <- vctrs::vec_unique(vapply(tree, digest::digest, ""))
  class <- "tree_vctr"

  vctrs::new_rcrd(
    fields = list(
      which_tree = vctrs::new_factor(which_tree, levels = hashes),
      node = node
    ),
    tree = tree,
    node_encoding = match.arg(node_encoding, c("order", "observation")),
    class = class
  )
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
  encode <- node_encoding(x)
  cat(
    sep = "",
    "# encoding by ", encode, "\n# n tree(s): ", n_trees, "\n"
  )
}

.tree_registry <- new.env(parent = emptyenv())

register_tree <- function(tree) {
  tree_hash <- digest::digest(tree)
  if (is.null(.tree_registry[[tree_hash]])) {
    .tree_registry[[tree_hash]] <- tree
  }
  tree_hash
}

unregister_hash <- function(hash) {
  val <- .tree_registry[[hash]]
  if (is.null(val)) {
    stop(
      sprintf(
        "cannot find hash %s. was it possible these were called out of order?",
        hash
      )
    )
  }
  rm(list = hash, envir = .tree_registry)
  val
}

new_tree_vctr_ptype <- function(learned = character(), encoding = "order") {
  vctrs::new_vctr(
    integer(),
    learned = learned,
    encoding = encoding,
    class = "tree_vctr_ptype"
  )
}


into_tree_vctr_ptype <- function(tree_vec) {
  new_tree_vctr_ptype(
    learned = vapply(trees(tree_vec), register_tree, ""),
    encoding = node_encoding(tree_vec)
  )
}

#' @importFrom vctrs vec_ptype
#' @export
vec_ptype.tree_vctr <- function(x, ...) {
  into_tree_vctr_ptype(x)
}

#' @importFrom vctrs vec_ptype_finalise
#' @export
vec_ptype_finalise.tree_vctr_ptype <- function(x, ...) {
  learned <- vctrs::vec_unique(attr(x, "learned"))
  tree <- lapply(learned, unregister_hash)
  new_tree_vctr(
    tree = tree,
    node_encoding = attr(x, "encoding")
  )
}

#' @importFrom vctrs vec_restore
#' @export
vec_restore.tree_vctr <- function(x, to, ...) {
  #
  if (vctrs::vec_size(x)) {
    # 4, 2, 3, 4, 2
    ints <- as.integer(x$which_tree)
    # 4, 2, 3
    unique_tree_id <- vctrs::vec_unique(ints)
    # 2, 3, 1
    unique_ord <- vctrs::vec_order(unique_tree_id)
    # 2, 3, 4
    unique_tree_id <- vctrs::vec_slice(
      unique_tree_id,
      unique_ord[!is.na(unique_ord)]
    )
    to_trees <- trees(to)
    if (!all(is.na(unique_tree_id))) {
      to_trees <- vctrs::vec_slice(trees(to), unique_tree_id)
    }
    new_tree_vctr(
      node = x$node,
      # 3, 1, 2, 3, 1
      which_tree = vctrs::vec_match(ints, unique_tree_id),
      tree = to_trees,
      node_encoding = node_encoding(to)
    )
  } else {
    new_tree_vctr(
      node = x$node,
      which_tree = as.integer(x$which_tree),
      tree = trees(to),
      node_encoding = node_encoding(to)
    )
  }
}

#' @importFrom vctrs vec_ptype2
#' @export
vec_ptype2.tree_vctr_ptype.tree_vctr_ptype <- function(x, y, ...) {
  new_tree_vctr_ptype(
    learned = vctrs::vec_c(
      attr(x, "learned"), attr(y, "learned"),
      .ptype = character()
    ),
    encoding = attr(x, "encoding")
  )
}

#' @importFrom vctrs vec_cast
#' @export
vec_cast.tree_vctr.tree_vctr <- function(x, to, ...) {
  #
  vctrs::field(x, "which_tree") <- vctrs::allow_lossy_cast(
    vctrs::vec_cast(
      vctrs::field(x, "which_tree"),
      vctrs::field(to, "which_tree")
    )
  )

  if (!identical(node_encoding(x), node_encoding(to))) {
    x <- swap_node_encoding(x)
  }
  x
}


#' @export
as.character.tree_vctr <- function(x) {
  fn <- switch(node_encoding(x),
    observation = function(node, tree) {
      n <- length(tree$order)
      labels <- sprintf("obsv%i", node)
      is_leaf <- node <= n
      tree_labels <- tree$labels
      if (!is.null(tree_labels)) {
        labels[is_leaf] <- tree_labels[node][is_leaf]
      }
      labels
    },
    order = function(node, tree) {
      n <- length(tree$order)
      labels <- sprintf("node%i", node)
      is_leaf <- node <= n
      tree_labels <- tree$labels
      if (!is.null(tree_labels)) {
        tree_labels <- tree_labels[tree$order]
        labels[is_leaf] <- tree_labels[node][is_leaf]
      }
      labels
    }
  )
  labels <- with_tree_vctr(x, fn)
  sprintf("%s-%i", labels, tree_id(x))
}












# .on_load <- function(ns) {
#   box::register_S3_method(name = "format", class = "tree_vctr", format.tree_vctr)
#   box::register_S3_method(name = "obj_print_footer", class = "tree_vctr", obj_print_footer.tree_vctr)
#   box::register_S3_method(name = "vec_restore", "tree_vctr", vec_restore.tree_vctr)
#   # box::register_S3_method(name = "vec_ptype2", "tree_vctr.tree_vctr", vec_ptype2.tree_vctr.tree_vctr)
#   box::register_S3_method(name = "vec_cast", "tree_vctr.tree_vctr", vec_cast.tree_vctr.tree_vctr)
#   # box::register_S3_method(name = "vec_cast", "tree_vctr.partial_tree_vctr", vec_cast.tree_vctr.partial_tree_vctr)
#   # box::register_S3_method(name = "vec_cast", "partial_tree_vctr.tree_vctr", vec_cast.partial_tree_vctr.tree_vctr)
#   box::register_S3_method(
#     name = "group_by_inner", "data.frame", group_by_inner.data.frame
#   )
#   box::register_S3_method(
#     name = "as.character", "tree_vctr", as.character.tree_vctr
#   )
#   box::register_S3_method(
#     name = "vec_ptype", "tree_vctr", vec_ptype.tree_vctr
#   )
#   box::register_S3_method(
#     name = "vec_ptype_finalise", "tree_vctr_ptype", vec_ptype_finalise.tree_vctr_ptype
#   )
#   box::register_S3_method(
#     name = "vec_ptype2", "tree_vctr_ptype.tree_vctr_ptype", vec_ptype2.tree_vctr_ptype.tree_vctr_ptype
#   )
#   # box::register_S3_method(
#   #   name = "vec_ptype2", "tree_vctr.partial_tree_vctr", vec_ptype2.tree_vctr.partial_tree_vctr
#   # )
#   # box::register_S3_method(
#   #   name = "vec_ptype2", "partial_tree_vctr.tree_vctr", vec_ptype2.partial_tree_vctr.tree_vctr
#   # )
#   # box::register_S3_method(
#   #   name = "vec_proxy_equal", "tree_vctr",
#   #   vec_proxy_equal.tree_vctr
#   # )
#   # box::register_S3_method(
#   #   name = "vec_proxy_compare", "tree_vctr",
#   #   vec_proxy_compare.tree_vctr
#   # )
# }

# if (is.null(box::name())) {
#   box::use(. / `__test_tree_vctr__`)
# }
