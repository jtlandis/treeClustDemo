# box::use(S7[new_generic, S7_dispatch, `method<-`])
# box::use(
#   . / bp[...],
#   ggplot2[...],
#   mods / util[vec_rep],
#   mods / tree[dendro_data],
#   stats[dist, hclust],
# )
# box::use(ggside[geom_ysidesegment]) |>
#   suppressMessages()
# ## classes ----


## Generics ----

#' @importFrom S7 S7_dispatch
#' @importFrom S7 new_generic
NULL

#' generates hclust from data
#' @param data an object
#' @param similarity_method one of 'sorensen', 'basepair',
#' or 'counts'.
#' @param dist_method method to provide to `stats::dist`
#' @param clust_method method to provide to `stats:hclust`
#'
#' @export
cluster_data <- new_generic(
  "cluster_data",
  "data",
  fun = function(data,
                 similarity_method = c("basepair_sorensen", "basepair_max_union", "basepair", "exon_sorensen", "exon"),
                 dist_method = "euclidean",
                 clust_method = "complete", name = NULL) {
    S7_dispatch()
  }
)


#' @export
plot_cluster <- new_generic(
  "plot_cluster",
  "data",
  fun = function(data,
                 similarity_method = c("basepair_sorensen", "basepair_max_union", "basepair", "exon_sorensen", "exon"),
                 dist_method = "euclidean",
                 clust_method = "complete", name = NULL) {
    S7_dispatch()
  }
)

Matrix_outer <- function(M, Fun = "*") {
  Fun <- match.fun(Fun)
  Fun(M, Matrix::t(M))
}

diag2 <- function(mat) {
  l <- nrow(mat)
  mat[(1:l) * l - l + 1:l]
}

sorensen_scale <- function(mat) {
  diag_vec <- diag2(mat)

  # non-zero locations
  logic_mat <- mat > 0
  n_per_col <- Matrix::colSums(logic_mat)

  mat_dims <- mat
  mat_dims[logic_mat] <- vctrs::vec_rep(diag_vec, n_per_col)

  mat[logic_mat] <- 2 * mat[logic_mat] /
    (Matrix_outer(mat_dims, `+`)[logic_mat])

  mat
}

max_scale <- function(mat) {
  diag_vec <- diag2(mat)

  logic_mat <- mat > 0
  n_per_col <- Matrix::colSums(logic_mat)

  mat_dims <- mat
  mat_dims[logic_mat] <- vctrs::vec_rep(diag_vec, n_per_col)

  mat[logic_mat] <- mat[logic_mat] /
    (Matrix_outer(mat_dims, function(x, y) {
      x_smaller_y <- x < y
      x[x_smaller_y] <- y[x_smaller_y]
      x
    })[logic_mat])

  mat
}

with_exons <- function(list_obj, with) {
  n_transcripts <- length(list_obj)
  chrom_exons <- plyxp::list_unchop(unname(list_obj))
  exon_mat <- with(chrom_exons)
  n_exons <- length(chrom_exons)
  mat <- Matrix::sparseMatrix(
    seq_len(n_exons),
    vctrs::vec_rep(seq_len(n_transcripts), lengths(list_obj)),
    x = 1L,
    dimnames = list(NULL, names(list_obj))
  )
  Matrix::t(mat) %*% exon_mat %*% mat
}

compute_simi <- function(object, with) {
  if (methods::is(object, "list")) {
    mat <- with_exons(object, with = with)
  } else {
    mat <- with(object)
  }
  dimnames(mat) <- list(names(object), names(object))
  mat
}

bp_simi_width <- function(object) {
  compute_simi(object = object, with = compute_bp_simi_width)
}

bp_sorensen_width <- function(object) {
  sorensen_scale(
    bp_simi_width(object)
  )
}

bp_max_union_width <- function(object) {
  max_scale(
    bp_simi_width(object)
  )
}

bp_simi_percent <- function(object) {
  compute_simi(object = object, with = compute_bp_simi_percent)
}

bp_sorensen_percent <- function(object) {
  sorensen_scale(
    bp_simi_percent(object)
  )
}

bp_max_union_percent <- function(object) {
  max_scale(
    bp_simi_percent(object)
  )
}


similarity_opts <- c(
  "basepair",
  "basepair_sorensen",
  "basepair_max-union",
  "percent",
  "percent_sorensen",
  "percent_max-union"
)

parse_method_opt <- function(method_name) {
  nn <- length(method_name)
  if (nn == length(similarity_opts) &&
    all.equal(method_name, similarity_opts)) {
    return(method_name[1])
  }
  if (nn != 1L) {
    stop("`similarity_method` should be length 1")
  }
  method_type <- match.arg(
    unique(sub("_.*$", "", method_name)),
    c("basepair", "percent")
  )

  method_scale <- sub("[^_]+", "", method_name)
  if (nchar(method_scale) != 0L) {
    method_scale <- match.arg(
      method_scale,
      c("_sorensen", "_max-union"),
    )
  }

  sprintf("%s%s", method_type, method_scale)
}

select_method <- function(method_name) {
  switch(nm <- parse_method_opt(method_name),
    basepair = bp_simi_width,
    basepair_sorensen = bp_sorensen_width,
    `basepair_max-union` = bp_max_union_width,
    percent = bp_simi_percent,
    percent_sorensen = bp_sorensen_percent,
    `percent_max-union` = bp_max_union_percent,
    stop(sprintf("no method for `%s`", nm))
  )
}

cluster_data_default <- function(
  data,
  similarity_method = c(
    "basepair",
    "basepair_sorensen",
    "basepair_max-union",
    "percent",
    "percent_sorensen",
    "percent_max-union"
  ),
  dist_method = "euclidean",
  clust_method = "complete", name = NULL
) {
  simi_fun <- select_method(similarity_method)
  dist_mat <- dist(simi_fun(data), method = dist_method)
  clust <- hclust(dist_mat, method = clust_method)
  clust
}

S7::method(
  cluster_data,
  S7::new_union(GR, GRL, class_simple_range, class_nested_range)
) <- cluster_data_default

S7::method(
  cluster_data,
  RSE
) <- function(
  data,
  similarity_method = c(
    "basepair",
    "basepair_sorensen",
    "basepair_max-union",
    "percent",
    "percent_sorensen",
    "percent_max-union"
  ),
  dist_method = "euclidean",
  clust_method = "complete", name = NULL
) {
  cluster_data_default(
    SummarizedExperiment::rowRanges(data),
    similarity_method = similarity_method,
    dist_method = dist_method,
    clust_method = clust_method
  )
}


# method(plot_cluster, GR) <-
#   function(data,
#            similarity_method = c("basepair_sorensen", "basepair_max_union", "basepair", "exon_sorensen", "exon"),
#            dist_method = "euclidean",
#            clust_method = "complete", name = NULL) {
#     clust <- gr_cluster_data(
#       data = data, similarity_method = similarity_method,
#       dist_method = dist_method, clust_method = clust_method
#     )
#     similarity_method <- parse_method_opt(similarity_method)

#     data_clust <- dendro_data(clust)

#     plot_data <- tibble::tibble(
#       start = IRanges::start(data),
#       end = IRanges::end(data),
#       seqnames = as.character(GenomicRanges::seqnames(data))
#     )
#     if (!is.null(name)) {
#       if (anyDuplicated(clust$labels)) {
#         clust$labels <- make.unique(clust$labels, sep = "_")
#       }
#       plot_data$cluster_order <- factor(
#         x = clust$labels, levels = clust$labels[clust$order]
#       )
#     } else {
#       plot_data <- plot_data[clust$order, ]
#       plot_data$cluster_order <- seq_len(nrow(plot_data))
#       plot_data <- plot_data[order(clust$order), ]
#     }

#     p <- ggplot(plot_data, aes(x = start, y = cluster_order)) +
#       geom_segment(aes(xend = end, yend = cluster_order)) +
#       geom_ysidesegment(
#         aes(x = y, y = x, xend = yend, yend = xend),
#         data = data_clust,
#         inherit.aes = FALSE
#       ) +
#       theme_bw() +
#       labs(
#         x = "Genome Position",
#         y = "Cluster Label",
#         subtitle = sprintf(
#           "Simularity: %s | dist: %s | clust: %s",
#           similarity_method,
#           dist_method, clust_method
#         )
#       )
#     if ("exons" %in% names(data@elementMetadata)) {
#       exon_ <- unlist(data$exons)
#       exon_data <- tibble::tibble(
#         start = exon_@ranges@start,
#         end = start + exon_@ranges@width - 1L,
#         seqnames = as.character(exon_@seqnames),
#         cluster_order = vec_rep(plot_data$cluster_order, vapply(data$exons, length, 1L))
#       )
#       p <- p +
#         geom_rect(
#           aes(
#             xmin = start, ymin = as.integer(cluster_order) - .25,
#             xmax = end, ymax = as.integer(cluster_order) + .25
#           ),
#           data = exon_data,
#           fill = "blue", inherit.aes = F
#         )
#     }
#     p
#   }

# method(plot_cluster, GRL) <-
#   function(data,
#            similarity_method = c("basepair_sorensen", "basepair_max_union", "basepair", "exon_sorensen", "exon"),
#            dist_method = "euclidean",
#            clust_method = "complete", name = NULL) {
#     clust <- gr_cluster_data(
#       data = data, similarity_method = similarity_method,
#       dist_method = dist_method, clust_method = clust_method
#     )

#     data_clust <- dendro_data(clust)
#     similarity_method <- parse_method_opt(similarity_method)
#     plot_data <- tibble::tibble(
#       start = min(IRanges::start(data)),
#       end = max(IRanges::end(data)),
#       seqnames = as.character(unique(GenomicRanges::seqnames(data)))
#     )
#     if (!is.null(name)) {
#       if (anyDuplicated(clust$labels)) {
#         clust$labels <- make.unique(clust$labels, sep = "_")
#       }
#       plot_data$cluster_order <- factor(
#         x = clust$labels, levels = clust$labels[clust$order]
#       )
#     } else {
#       plot_data <- plot_data[clust$order, ]
#       plot_data$cluster_order <- seq_len(nrow(plot_data))
#       plot_data <- plot_data[order(clust$order), ]
#     }

#     p <- ggplot(plot_data, aes(x = start, y = cluster_order)) +
#       geom_segment(aes(xend = end, yend = cluster_order)) +
#       geom_ysidesegment(
#         aes(x = y, y = x, xend = yend, yend = xend),
#         data = data_clust,
#         inherit.aes = FALSE
#       ) +
#       theme_bw() +
#       labs(
#         x = "Genome Position",
#         y = "Cluster Label",
#         subtitle = sprintf(
#           "Simularity: %s | dist: %s | clust: %s",
#           similarity_method,
#           dist_method, clust_method
#         )
#       )
#     exons <- unlist(data)
#     exon_data <- tibble::tibble(
#       start = IRanges::start(exons),
#       end = IRanges::end(exons),
#       seqnames = as.character(GenomicRanges::seqnames(exons)),
#       cluster_order = vec_rep(plot_data$cluster_order, lengths(data))
#     )
#     p <- p +
#       geom_rect(
#         aes(
#           xmin = start, ymin = as.integer(cluster_order) - .25,
#           xmax = end, ymax = as.integer(cluster_order) + .25
#         ),
#         data = exon_data,
#         fill = "blue", inherit.aes = FALSE
#       )
#     p
#   }

# method(plot_cluster, RSE) <- function(
#   data,
#   similarity_method = c("basepair_sorensen", "basepair_max_union", "basepair", "exon_sorensen", "exon"),
#   dist_method = "euclidean",
#   clust_method = "complete", name = NULL
# ) {
#   plot_cluster(
#     SummarizedExperiment::rowRanges(data),
#     similarity_method, dist_method, clust_method, name
#   )
# }
