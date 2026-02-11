#' @export
group_by_tree_nodes <- function(x, .by, ..., only_inner = TRUE) {
  UseMethod("group_by_tree_nodes")
}

#' @export
group_by_tree_nodes.data.frame <- function(
  x, .by, ..., only_inner = TRUE,
  node_level = NULL
) {
  by_tree <- dplyr::select(x, {{ .by }})
  stopifnot(
    "`.by` should only select one variable" = length(by_tree) == 1L,
    "`.by` should be a 'tree_vctr'" = inherits(by_tree[[1]], "tree_vctr")
  )
  node_level <- rlang::eval_tidy({{ node_level }}, data = x)
  groups <- generate_inner_slice(by_tree[[1]],
    node_level = node_level,
    only_inner = only_inner
  ) |>
    dplyr::rename(
      "{names(by_tree)}" := nodes,
      .rows = children
    )
  dplyr::new_grouped_df(x, groups = groups)
}

#' @export
group_by_nodes <- function(x, ...) {
  UseMethod("group_by_nodes")
}

#' @export
group_by_nodes.data.frame <- function(x, ...) {
  quos <- rlang::enquos(..., .named = TRUE)
  result <- lapply(quos, rlang::eval_tidy, data = x)
  groups <- nodes_grouped(result, data = x)
  dplyr::new_grouped_df(x, groups = groups)
}

#' @export
group_by_nodes.PlySummarizedExperiment <- function(x, ...) {
  quos <- plyxp:::plyxp_quos(..., .ctx = c("assays", "rows", "cols"))
  mask <- plyxp:::new_plyxp_manager(plyxp::se(x))
  ctxs <- vapply(quos, attr, FUN.VALUE = "", which = "plyxp:::ctx")
  nms <- names(quos)
  mask <- plyxp:::plyxp_evaluate(mask, quos, ctxs, nms, rlang::caller_env())
  results <- mask$results()
  groups <- structure(
    list(
      row_groups = nodes_grouped_DF(
        results$rows, SummarizedExperiment::rowData(x)
      ),
      col_groups = nodes_grouped_DF(
        results$cols, SummarizedExperiment::colData(x)
      )
    ),
    class = "plyxp_groups"
  )
  se_obj <- plyxp::se(x)
  plyxp:::group_data_se_impl(se_obj) <- groups
  plyxp::se(x) <- se_obj
  x
}

nodes_grouped <- function(result, data) {
  len <- length(result)
  if (is.null(result) || len == 0L) {
    return(NULL)
  }
  if (len > 1L) {
    rlang::warn("`group_by_nodes()` should have one result per context. Using first.")
  }
  result <- result[[1L]]
  ptype <- result[0L]
  if (!inherits(ptype, "tree_vctr")) {
    stop("`group_by_nodes()` should result in a `tree_vctr`")
  }
  pos <- base::Position(
    \(obj) {
      inherits(obj, "tree_vctr") && identical(obj[0], ptype)
    },
    x = data
  )
  if (is.na(pos)) {
    stop("tree_vctr result does not match any data")
  }
  data <- data[pos]
  generate_inner_slice(data[[1]],
    node_level = result,
    only_inner = FALSE
  ) |>
    dplyr::rename(
      "{names(data)}" := nodes,
      .indices = children
    )
}

nodes_grouped_DF <- function(result, data) {
  x <- nodes_grouped(result, data)
  if (is.null(x)) {
    return(NULL)
  }
  x |>
    dplyr::mutate(.indices_group_id = seq_len(dplyr::n())) |>
    methods::as("DFrame")
}

#' @export
group_by_tree_nodes.PlySummarizedExperiment <- function(
  x, .by, ...,
  only_inner = TRUE,
  row_node_level = NULL, col_node_level = NULL
) {
  by_tree <- dplyr::select(x, {{ .by }})
  row_data <- SummarizedExperiment::rowData(by_tree)
  nrow_sel <- length(row_data)
  quos <- plyxp:::plyxp_quos(
    {{ row_node_level }},
    {{ col_node_level }},
    .ctx = c("assays", "rows", "cols")
  )
  mask <- plyxp:::new_plyxp_manager(plyxp::se(x))
  ctxs <- vapply(quos, attr, FUN.VALUE = "", which = "plyxp:::ctx")
  nms <- names(quos)
  mask <- plyxp:::plyxp_evaluate(mask, quos, ctxs, nms, rlang::caller_env())
  results <- mask$results()
  mask <- plyxp:::new_plyxp_manager(plyxp::se(x))
  if (nrow_sel != 0) {
    stopifnot(
      "`.by` should only select one rowData variable" = nrow_sel == 1L,
      "`.by` should be a 'tree_vctr'" = inherits(row_data[[1]], "tree_vctr")
    )
    node_level <- if (length(results$rows)) {
      results$rows[[1]]
    } else {
      NULL
    }
    row_data <- generate_inner_slice(row_data[[1]],
      node_level = node_level,
      only_inner = only_inner
    ) |>
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
    node_level <- if (length(results$cols)) {
      results$cols[[1]]
    } else {
      NULL
    }
    col_data <- generate_inner_slice(col_data[[1]],
      node_level = node_level,
      only_inner = only_inner
    ) |>
      dplyr::rename(
        "{names(col_data)}" := nodes,
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
summarize_by_tree_nodes <- function(x, ..., .by) {
  x <- group_by_tree_nodes(x, {{ .by }})
  dplyr::summarise(x, ...)
}
