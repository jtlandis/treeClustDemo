ScaleDiscreteTree <- ggplot2::ggproto(
  "ScaleDiscreteTree",
  ggplot2::ScaleDiscretePosition,
  train = function(self, x) {
    if (inherits(x, "tree_vctr")) {
      x <- as.factor(x)
    }
    parent <- ggplot2::ggproto_parent(ggplot2::ScaleDiscretePosition, self)
    parent$train(x = x)
  },
  map = function(self, x, limits = self$get_limits()) {
    original_x <- x
    if (is_tree <- inherits(x, "tree_vctr")) {
      self$tree_ptype <- vctrs::vec_ptype2(self$tree_ptype, x)
      original_x <- x <- vctrs::vec_cast(x, self$tree_ptype)
      x <- as.factor(x)
    }
    parent <- ggplot2::ggproto_parent(ggplot2::ScaleDiscretePosition, self)
    mapped <- parent$map(x = x, limits = limits)
    if (is_tree) {
      mapped <- new_tree_mapped(mapped, original_x)
    }
    mapped
  },
  tree_ptype = NULL
)

#' @export
scale_y_tree <- function(
  name = ggplot2::waiver(), ..., palette = seq_len, expand = ggplot2::waiver(),
  guide = ggplot2::waiver(), position = "left", sec.axis = ggplot2::waiver(),
  continuous.limits = NULL
) {
  sc <- ggplot2::discrete_scale(
    aesthetics = ggplot2:::ggplot_global$y_aes, name = name,
    palette = palette, ..., expand = expand, guide = guide,
    position = position, super = ScaleDiscreteTree
  )
  sc$range_c <- scales::ContinuousRange$new()
  sc$continuous_limits <- continuous.limits
  ggplot2:::set_sec_axis(sec.axis, sc)
}


#' @export
scale_x_tree <- function(
  name = ggplot2::waiver(), ..., palette = seq_len, expand = ggplot2::waiver(),
  guide = ggplot2::waiver(), position = "bottom", sec.axis = ggplot2::waiver(),
  continuous.limits = NULL
) {
  sc <- ggplot2::discrete_scale(
    aesthetics = ggplot2:::ggplot_global$x_aes, name = name,
    palette = palette, ..., expand = expand, guide = guide,
    position = position, super = ScaleDiscreteTree
  )
  sc$range_c <- scales::ContinuousRange$new()
  sc$continuous_limits <- continuous.limits
  ggplot2:::set_sec_axis(sec.axis, sc)
}

#' @importFrom ggplot2 scale_type
#' @export
scale_type.tree_vctr <- function(x) "tree"


StatHcluster <- ggplot2::ggproto(
  "StatHcluster",
  ggplot2::Stat,
  required_aes = c("xcluster", "ycluster", "cluster"),
  default_aes = ggplot2::aes(
    x = after_stat(xcluster),
    y = after_stat(ycluster)
  ),
  setup_data = function(data, params) {
    # if (!"xcluster" %in% names(data)) {
    #   if (!"x" %in% names(data)) {
    #     stop("StatHcluster requires an x aesthetic mapping.", call. = FALSE)
    #   }
    # }
    data$group <- ggplot2:::add_group(data[!names(data) %in% c("group", "xcluster", "ycluster")])$group
    data
  },
  compute_group = function(
    data, scales, distance_method = "euclidean", p = 2, linkage_method = "complete"
  ) {
    wide <- tidyr::pivot_wider(data,
      id_cols = xcluster, names_from = ycluster, values_from = cluster
    )
    mat <- as.matrix(wide[-1L])
    if (is.null(scales$x)) {
      rownames(mat) <- wide[[1L]]
    } else {
      rownames(mat) <- scales$x$range$range[wide[[1L]]]
    }
    if (!is.null(scales$y)) {
      colnames(mat) <- scales$y$range$range[as.integer(colnames(mat))]
    }
    tree_x <- hclust(dist(mat, method = distance_method, p = p),
      method = linkage_method
    ) |> as_tree_vctr()
    data$xcluster <- tree_x[vctrs::vec_match(
      data$xcluster,
      tree_labels(tree_x)
    )]
    tree_y <- hclust(dist(t(mat), method = distance_method, p = p),
      method = linkage_method
    ) |> as_tree_vctr()
    data$ycluster <- tree_y[vctrs::vec_match(
      data$ycluster,
      tree_labels(tree_y)
    )]
    data
  },
  compute_layer = function(self, data, params, layout) {
    browser()
    parent <- ggplot2::ggproto_parent(ggplot2::Stat, self)
    data <- parent$compute_layer(data = data, params = params, layout = layout)
    data[c("xcluster", "ycluster")] <- lapply(
      data[c("xcluster", "ycluster")],
      \(tree) new_tree_mapped(vec_cast(tree, double()), tree)
    )
    data
  }
)

geom_heatmap <- function(
  mapping = NULL, data = NULL, stat = "hcluster", position = "identity",
  ..., lineend = "butt", linejoin = "mitre", na.rm = FALSE,
  show.legend = NA, inherit.aes = TRUE
) {
  ggplot2::layer(
    mapping = mapping, data = data, geom = "tile", stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = rlang::list2(
      na.rm = na.rm, lineend = lineend, linejoin = linejoin,
      ...
    )
  )
}

#' mat <- matrix(0, 6, 5)
#' mat[1:3, 1:3] <- rnorm(9, 10)
#' mat[1:2, 4:5] <- rnorm(4, 4)
#' mat[4:6,1:2] <- rnorm(6, 7)
#' mat[4:6, 3:5] <- rnorm(9, 2)
#' mat[3, 4:5] <- rnorm(2, 15)
#' dimnames(mat) <- list(sprintf("rows%i", 1:6), sprintf("cols%i", 1:5))
#' df <- tibble::as_tibble(mat, rownames = "rows") |>
#'  dplyr::mutate(grow = rep(1:2, each = 3)) |>
#'  tidyr::pivot_longer(cols = -c(rows, grow), names_to = "cols") |>
#'  dplyr::mutate(gcol = rep(c(1,1,1,2,2), 6))
