#' @export
group_by_inner <- function(x, .by, node_level = NULL) {
  UseMethod("group_by_inner")
}

#' @export
group_by_inner.data.frame <- function(x, .by, node_level = NULL) {
  by_tree <- dplyr::select(x, {{ .by }})
  stopifnot(
    "`.by` should only select one variable" = length(by_tree) == 1L,
    "`.by` should be a 'tree_vctr'" = inherits(by_tree[[1]], "tree_vctr")
  )
  groups <- generate_inner_slice(by_tree[[1]], node_level = node_level) |>
    dplyr::rename(
      "{names(by_tree)}" := nodes,
      .rows = children
    )
  dplyr::new_grouped_df(x, groups = groups)
}

#' @export
group_by_inner.PlySummarizedExperiment <- function(x, .by, node_level = NULL) {
  by_tree <- dplyr::select(x, {{ .by }})
  row_data <- SummarizedExperiment::rowData(by_tree)
  nrow_sel <- length(row_data)
  if (nrow_sel != 0) {
    stopifnot(
      "`.by` should only select one rowData variable" = nrow_sel == 1L,
      "`.by` should be a 'tree_vctr'" = inherits(row_data[[1]], "tree_vctr")
    )
    row_data <- generate_inner_slice(row_data[[1]], node_level = node_level) |>
      dplyr::rename(
        "{names(row_data)}" := nodes,
        .indices = children
      ) |>
      dplyr::mutate(.indices_group_id = seq_len(dplyr::n())) |>
      methods::as("DFrame")
  } else {
    row_data <- NULL
  }
  col_data <- SummarizedExperiment::colData(by_tree)
  ncol_sel <- length(col_data)
  if (ncol_sel != 0) {
    stopifnot(
      "`.by` should only select one colData variable" = ncol_sel == 1L,
      "`.by` should be a 'tree_vctr'" = inherits(col_data[[1]], "tree_vctr")
    )
    col_data <- generate_inner_slice(col_data[[1]], node_level = node_level) |>
      dplyr::rename(
        "{names(row_data)}" := nodes,
        .indices = children
      ) |>
      dplyr::mutate(.indices_group_id = seq_len(dplyr::n())) |>
      methods::as("DFrame")
  } else {
    col_data <- NULL
  }
  groups <- structure(
    list(
      row_groups = row_data,
      col_groups = col_data
    ),
    class = "plyxp_groups"
  )
  se_obj <- plyxp::se(x)
  plyxp:::group_data_se_impl(se_obj) <- groups
  plyxp::se(x) <- se_obj
  x
}

#' @export
summarize_by_inner <- function(x, ..., .by) {
  x <- group_by_inner(x, {{ .by }})
  dplyr::summarise(x, ...)
}
