#######################
#  Old code, needs refactor
#######################

tree_mapped <- function(tree) {
  new_tree_mapped(
    vctrs::vec_cast(tree, double()),
    tree = tree
  )
}

new_tree_mapped <- function(data = double(), tree = list()) {
  vctrs::new_vctr(
    .data = data,
    tree_vctr = tree,
    class = "tree_mapped",
    inherit_base_type = TRUE
  )
}

# #' @export
# vec_restore.tree_mapped <- function(x, to, ...) {
#   # browser()
#   if (vctrs::vec_size(x)) {
#     shouldb <- tryCatch(
#       {
#         i <- eval.parent(quote((\() i)()))
#         FALSE
#       },
#       error = function(e) TRUE
#     )
#     if (shouldb) {
#       browser()
#     }
#   } else {
#     i <- 0L
#   }

#   new_tree_mapped(
#     x,
#     tree = vctrs::vec_slice(attr(to, "tree_vctr"), i)
#   )
# }

vec_ptype2.tree_mapped.tree_mapped <- function(x, y, ...) {
  new_tree_mapped(
    tree = vec_ptype2(attr(x, "tree_vctr"), attr(y, "tree_vctr"))
  )
}

vec_cast.tree_mapped.tree_mapped <- function(x, to, ...) {
  new_tree_mapped(
    vctrs::vec_data(x),
    tree = vec_cast(
      attr(x, "tree_vctr"),
      attr(to, "tree_vctr"),
      ...
    )
  )
}


#' @export
vec_ptype2.tree_mapped.double <- function(x, y, ...) {
  double()
}

#' @export
vec_ptype2.double.tree_mapped <- function(x, y, ...) {
  double()
}

#' @export
vec_cast.double.tree_mapped <- function(x, to, ...) {
  vctrs::vec_data(x)
}


# tree_mapped_vctr_class <- c(
#   "tree_mapped", "tree_vctr", "vctrs_rcrd", "vctrs_vctr"
# )

# tree_mapped <- function(tree_vec) {
#   new_tree_mapped(vec_cast(tree_vec, integer()), tree_vec)
# }

# #' @export
# vec_proxy.tree_mapped <- function(x, ...) {
#   class(x) <- tree_vctr_class
#   vec_proxy(x)
# }

# #' @export
# vec_restore.tree_mapped <- function(x, to, ...) {
#   new_tree_mapped(
#     x$position,
#     new_pure_tree_vctr(
#       node = x$node,
#       which_tree = x$which_tree,
#       tree = trees(to),
#       node_encoding = node_encoding(to)
#     )
#   )
# }

# new_tree_mapped <- function(position, tree) {
#   tree <- unclass(tree)
#   tree$position <- position
#   class(tree) <- tree_mapped_vctr_class
#   tree
# }

tree_mapped_arith <- function(op, x, y) {
  y <- vec_recycle(y, vec_size(x))
  op_fn <- getExportedValue("base", op)
  op_fn(
    vec_data(x),
    vec_data(y)
  ) |> new_tree_mapped(tree = attr(x, "tree_vctr"))
}

#' @export
#' @method vec_arith tree_mapped
vec_arith.tree_mapped <- function(op, x, y, ...) {
  UseMethod("vec_arith.tree_mapped", y)
}

#' @export
vec_arith.tree_mapped.default <- function(op, x, y, ...) {
  vctrs::stop_incompatible_op(op, x, y)
}

#' @export
#' @method vec_arith.tree_mapped tree_mapped
vec_arith.tree_mapped.tree_mapped <- function(op, x, y, ...) {
  switch(op,
    "+" = ,
    "-" = ,
    "*" = ,
    "/" = tree_mapped_arith(op, x, y),
    vctrs::stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @importFrom vctrs vec_arith.numeric
#' @method vec_arith.numeric tree_mapped
vec_arith.numeric.tree_mapped <- function(op, x, y, ...) {
  switch(op,
    "+" = ,
    "-" = ,
    "*" = ,
    "/" = tree_mapped_arith(op, x, y),
    vctrs::stop_incompatible_op(op, x, y)
  )
}

#' @export
#' @method vec_arith.tree_mapped double
vec_arith.tree_mapped.double <- function(op, x, y, ...) {
  switch(op,
    "+" = ,
    "-" = ,
    "*" = ,
    "/" = tree_mapped_arith(op, x, y),
    vctrs::stop_incompatible_op(op, x, y)
  )
}


# #' @export
# vec_ptype2.tree_mapped.double <- function(x, y, ...) {
#   double()
# }

# #' @export
# vec_ptype2.double.tree_mapped <- function(x, y, ...) {
#   double()
# }

# #' @export
# vec_cast.double.tree_mapped <- function(x, to, ...) {
#   as.double(vctrs::field(x, "position"))
# }

# #' @export
# vec_cast.tree_mapped.integer <- function(x, to, ...) {
#   mtch <- leaf_labels(to, use_labels = FALSE)[x]
#   mtch <- sub("(node|obsv)", "", x = mtch)
#   split <- strsplit(mtch, "-")
#   nodes <- vapply(split, `[[`, FUN.VALUE = "", 1L) |> as.integer()
#   if (any(is_na <- is.na(mtch))) {
#     wt <- nodes
#     wt[!is_na] <- vapply(split[!is_na], `[[`, FUN.VALUE = "", 2L) |>
#       as.integer()
#   } else {
#     wt <- vapply(split, `[[`, FUN.VALUE = "", 2L) |>
#       as.integer()
#   }

#   new_pure_tree_vctr(
#     node = nodes,
#     which_tree = vctrs::new_factor(wt, levels(wt(to))),
#     tree = trees(to),
#     node_encoding = node_encoding(to)
#   )
# }

# #' @export
# vec_cast.tree_vctr.double <- function(x, to, ...) {
#   vec_cast(
#     vec_cast(x, integer()),
#     to,
#     ...
#   )
# }

# #' @export
# vec_ptype2.tree_vctr.integer <- function(x, y, ...) {
#   integer()
# }

# #' @export
# vec_ptype2.integer.tree_vctr <- function(x, y, ...) {
#   integer()
# }

# #' @export
# vec_cast.integer.tree_vctr <- function(x, to, ...) {
#   as.integer(as.factor(x))
# }
