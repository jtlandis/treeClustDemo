#' @export
as.data.frame.tree_vctr <- function(x, ...) {
  data <- with_tree_vctr(
    encode_obsv(x),
    function(nodes, tree) {
      merge <- tree$merge
      n <- length(tree$order)
      is_leaf <- merge < 0
      merge[is_leaf] <- -merge[is_leaf]
      merge[!is_leaf] <- merge[!is_leaf] + n
      N <- n + n - 1L
      x_pos <- double(N)
      x_last <- matrix(NA_real_,
        nrow = N, ncol = 2,
        dimnames = list(NULL, c("min", "max"))
      )
      y_pos <- double(N)
      y_next <- rep(NA_real_, N)
      y_pos[(n + 1):N] <- tree$height
      x_pos[seq_len(n)] <- order(tree$order)
      for (i in seq_len(nrow(merge))) {
        node <- n + i
        is <- merge[i, ]
        xx <- x_pos[is]
        x_last[node, ] <- xx
        x_pos[node] <- mean(xx)
        y_next[is] <- y_pos[node]
      }
      data <- data.frame(
        x = x_pos,
        y = y_pos,
        y_end = y_next,
        x_start = x_last[, "min"],
        x_end = x_last[, "max"]
      )
      data[nodes, ]
    }
  )
  dplyr::bind_cols(nodes = x, data, ...)
}


#' @export
plot_tree_vctr_data <- function(data, ...) {
  ggplot2::ggplot(data, ggplot2::aes(x, y, ...)) +
    ggplot2::geom_point() +
    ggplot2::geom_segment(
      ggplot2::aes(yend = y_end),
      ~ subset(.x, !is.na(y_end))
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = x_start, xend = x_end),
      ~ subset(.x, !is.na(x_end))
    )
}

#' @importFrom ggplot2 autoplot
#' @export
autoplot.tree_vctr <- function(object, ..., .other_data = list()) {
  plot_tree_vctr_data(
    as.data.frame(object, !!!.other_data), ...
  ) + ggplot2::facet_wrap(~ tree_id(nodes), scales = "free_x")
}
