#' @export
climb <- function(.data, method = c("strict", "ordered_pval"), ...) {
  method <- match.arg(method, c("strict", "ordered_pval"))
  f <- switch(
    method,
    strict = climb_strict,
    ordered_pval = climb_pvalue_arrange
  )

  f(.data, ...)
}

#' @description
#' climb an expanded results data with a tree.
#' @returns logical vector of which nodes to keep
climb_strict <- function(tree_vec, pvals, which = NULL, ...) {
  .clust <- tree(tree_vec, which = which) # should error if not length 1
  stopifnot(
    "`tree_vec` and `pvals` should be same size" = length(tree_vec) ==
      length(pvals),
    "expected an hclust `tree_vec`" = inherits(.clust, "hclust")
  )
  ord <- order(node_id(tree_vec))
  m <- .clust$merge
  o <- .clust$order
  n <- nrow(m) + 1
  is_merge <- m < 0
  m[is_merge] <- match(-m[is_merge], o)
  m[!is_merge] <- m[!is_merge] + n
  # pvals <- .data$pvalue
  keep <- rep(NA, length(pvals))
  is_na <- is.na(pvals)
  pvals[is_na] <- 1
  pvals <- pvals[ord]

  for (i in seq_len(nrow(m))) {
    nodes <- m[i, ]

    parent_node <- i + n
    # Before we check p-values, check if we need to keep
    # climbing. Specifically if the parent is FALSE
    # we do not need to compare anything.
    dont_use <- vapply(keep[nodes], isFALSE, logical(1))
    if (any(dont_use)) {
      keep[parent_node] <- FALSE
      # some data may be newly merged and havent been tested
      # we don't need to compare p-values for these if we know
      # this merge index wont be used
      keep[nodes[!dont_use]] <- TRUE
      next
    }

    node_p <- pvals[nodes]
    if (all(node_p > pvals[parent_node])) {
      keep[nodes] <- FALSE
      keep[parent_node] <- TRUE
    } else {
      keep[nodes] <- TRUE
      keep[parent_node] <- FALSE
    }
  }

  # .data$pvalue[is_na] <- NA_real_
  keep[order(ord)]
}

#' @param .data a result object
#' @param ... unused
climb_pvalue_arrange <- function(
  tree_vec,
  pvals,
  which = NULL,
  ...,
  merge_alpha = 1
) {
  n_leafs <- tree_leaf_max_node(tree(tree_vec, which = which))

  # leafs <- filter(res, n_children==1) |>
  #   tidyr::unnest(descendants)
  desc <- numeric()
  o_pval <- order(pvals)
  keep <- logical(length(tree_vec))
  pvals <- pvals[o_pval]
  tree_vec <- tree_vec[o_pval]
  i <- 0L
  theDesc <- descendants(tree_vec)
  nchild <- vapply(theDesc, length, integer(1))
  while (length(desc) < n_leafs) {
    i <- i + 1L
    desc2 <- theDesc[[i]]
    simi <- any(desc2 %in% desc)
    if (!simi && (nchild[i] == 1L || pvals[i] < merge_alpha)) {
      desc <- c(desc, desc2)
      keep[i] <- TRUE
    }
  }
  keep[order(o_pval)]
}
