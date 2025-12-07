GeomDendro <- ggplot2::ggproto("GeomDendro", ggplot2::GeomSegment,
  required_aes = c("x|y"),
  default_aes = ggplot2::aes(
    colour = from_theme(colour %||% ink),
    linewidth = from_theme(linewidth),
    linetype = from_theme(linetype),
    shape = from_theme(pointshape),
    size = from_theme(pointsize),
    stroke = from_theme(borderwidth),
    fill = from_theme(fill %||% NA),
    alpha = NA
  ),
  draw_panel = function(
    self, data, panel_params, coord, arrow = NULL,
    arrow.fill = NULL, lineend = "butt", linejoin = "round",
    na.rm = FALSE, add_points = TRUE
  ) {
    parent <- ggplot2::ggproto_parent(GeomDendro, self)
    # vertical bars
    parent$draw_panel(
      data = subset(data, !is.na(y_end)),
      panel_params = panel_params,
      coord = coord, arrow = arrow,
      arrow.fill = arrow.fill, lineend = lineend,
      linejoin = linejoin, na.rm = na.rm
    )
    parent$draw_panel(
      data = subset(data, !is.na(x_end)),
      panel_params = panel_params,
      coord = coord, arrow = arrow,
      arrow.fill = arrow.fill, lineend = lineend,
      linejoin = linejoin, na.rm = na.rm
    )
    if (add_points) {
      ggplot2::GeomPoint$draw_panel(
        data = data,
        panel_params = panel_params,
        coord = coord,
        na.rm = na.rm
      )
    }
  },
  extra_params = c("na.rm", "add_points"),
  draw_key = function(data, params, size) {
    ggplot2::draw_key_path(data, params, size)
    if (params$add_points) {
      ggplot2::draw_key_point(data, params, size)
    }
  },
  setup_data = function(data, params) {
    browser()
    flip <- "y" %in% names(data)
    vec <- data[[if (flip) "y" else "x"]]
    stopifnot(inherits(vec, "tree_mapped"))
    tree <- attr(vec, "tree_vctr")
    stopifnot(!is.null(tree), length(tree) == length(vec))
    unique_trees <- vctrs::vec_unique_loc(tree)
    data <- vctrs::vec_slice(data, unique_trees)
    expand <- vctrs::vec_unique(vctrs::vec_c(tree, inner(tree)))
    expand
  }
)

geom_dendro <- function(
  mapping = NULL, data = NULL, stat = "identity",
  position = "identity", ..., arrow = NULL, arrow.fill = NULL,
  lineend = "butt", linejoin = "round",
  na.rm = FALSE, show.legend = NA, inherit.aes = TRUE,
  add_points = TRUE
) {
  ggplot2::layer(
    mapping = mapping, data = data, geom = GeomDendro, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = rlang::list2(
      na.rm = na.rm, arrow = arrow, arrow.fill = arrow.fill,
      lineend = lineend, linejoin = linejoin, add_points = add_points, ...
    )
  )
}

# ggplot2::geom_point() +
# ggplot2::geom_segment(
#   ggplot2::aes(yend = y_end),
#   ~ subset(.x, !is.na(y_end))
# ) +
# ggplot2::geom_segment(
#   ggplot2::aes(x = x_start, xend = x_end),
#   ~ subset(.x, !is.na(x_end))
# )
