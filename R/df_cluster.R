assert_discrete <- function(x, arg_name) {
  if (length(x) != 1L) {
    cli::cli_abort("{.arg {arg_name}} should select one column.")
  }
  x <- x[[1]]
  if (is.factor(x)) {
    return(x)
  }
  if (!is.character(x)) {
    cli::cli_abort("{.arg {arg_name}} should select a discrete column.")
  }
  factor(x, unique(x))
}

vec_encode <- function(x, use_labels = FALSE) {
  unique_loc <- vctrs::vec_unique_loc(x)
  uniques <- vctrs::vec_slice(x, unique_loc)
  matches <- vctrs::vec_match(x, uniques)
  attr(matches, "unique_loc") <- unique_loc
  if (use_labels) {
    labels <- lapply(uniques, as.character)
    attr(matches, "labels") <- do.call("paste", labels)
  }
  matches
}

row_cols_id <- function(
  row_ids, col_ids, labeled = FALSE
) {
  row_ids <- vec_encode(row_ids, labeled)
  col_ids <- vec_encode(col_ids, labeled)
  row_uniq_loc <- attr(row_ids, "unique_loc")
  nrows <- length(row_uniq_loc)
  col_uniq_loc <- attr(col_ids, "unique_loc")
  ncols <- length(col_uniq_loc)
  row_ids_ <- row_ids - 1L
  col_ids_ <- col_ids - 1L
  mat_ids <- (col_ids_ * nrows) + row_ids_ + 1L
  list(
    df_to_mat = match(seq_len(nrows * ncols), mat_ids),
    mat_to_df = mat_ids,
    row_uniq_loc = row_uniq_loc,
    nrows = nrows,
    rol_uniq_loc = col_uniq_loc,
    ncols = ncols,
    row_labels = attr(row_ids, "labels"),
    col_labels = attr(col_ids, "labels")
  )
}

cluster_mat_by_rows <- function(
  data, values, by,
  .distance_method, .linkage_method, .p
) {
  if (length(by)) {
    row_tree <- vctrs::vec_init(new_pure_tree_vctr(), nrow(values))
    trees <- vctrs::vec_chop(
      values,
      indices = by
    ) |> lapply(
      function(mat) {
        hclust(
          dist(mat, method = .distance_method, p = .p),
          method = .linkage_method
        ) |> tree_vctr()
      }
    )
    row_tree <- vctrs::vec_init(
      vctrs::vec_ptype_common(!!!trees),
      nrow(values)
    )
    for (i in seq_along(trees)) {
      vctrs::vec_slice(row_tree, by[[i]]) <- trees[[i]]
    }
  } else {
    row_tree <- hclust(
      dist(
        values,
        method = .distance_method, p = .p
      ),
      method = .linkage_method
    ) |> tree_vctr()
  }
  row_tree
}

#' cluster data frame
#' @examples
#' mat <- matrix(0, 6, 5)
#' mat[1:3, 1:3] <- rnorm(9, 10)
#' mat[1:2, 4:5] <- rnorm(4, 4)
#' mat[4:6, 1:2] <- rnorm(6, 7)
#' mat[4:6, 3:5] <- rnorm(9, 2)
#' mat[3, 4:5] <- rnorm(2, 15)
#' dimnames(mat) <- list(sprintf("rows%i", 1:6), sprintf("cols%i", 1:5))
#' df <- tibble::as_tibble(mat, rownames = "rows") |>
#'   dplyr::mutate(grow = rep(1:2, each = 3)) |>
#'   tidyr::pivot_longer(cols = -c(rows, grow), names_to = "cols") |>
#'   dplyr::mutate(gcol = rep(c(1, 1, 1, 2, 2), 6))
#' @export
df_cluster <- function(
  data,
  id_rows,
  id_cols,
  values_from,
  rows_to = "row_tree",
  cols_to = "col_tree",
  ...,
  .linkage_method = "complete",
  .distance_method = "euclidean",
  .p = 2,
  cluster_row = TRUE,
  cluster_col = TRUE,
  by_rows = NULL,
  by_cols = NULL,
  use_labels = FALSE
) {
  key_map <- row_cols_id(
    row_ids = dplyr::select(data, {{ id_rows }}),
    col_ids = dplyr::select(data, {{ id_cols }}),
    labeled = use_labels
  )
  values <- dplyr::select(data, {{ values_from }})
  stopifnot(length(values) == 1L)

  values <- matrix(
    values[[1L]][key_map$df_to_mat],
    nrow = key_map$nrows,
    ncol = key_map$ncols,
    dimnames = list(key_map$row_labels, key_map$col_labels)
  )
  if (isTRUE(cluster_row)) {
    by_rows <- if (!missing(by_rows)) {
      row_ <- dplyr::select(data, {{ by_rows }}) |>
        vctrs::vec_slice(key_map$row_uniq_loc) |>
        vctrs::vec_group_loc() |>
        _$loc
      row_[lengths(row_) > 1L]
    } else {
      NULL
    }
    row_tree <- cluster_mat_by_rows(
      data,
      values = values,
      by = by_rows,
      .distance_method = .distance_method,
      .linkage_method = .linkage_method,
      .p = .p
    )
    data[[rows_to]] <- vctrs::vec_rep(row_tree, key_map$ncols) |>
      vctrs::vec_slice(key_map$mat_to_df)
  }

  if (isTRUE(cluster_col)) {
    by_cols <- if (!missing(by_cols)) {
      col_ <- dplyr::select(data, {{ by_cols }}) |>
        vctrs::vec_slice(key_map$col_uniq_loc) |>
        vctrs::vec_group_loc() |>
        _$loc
      col_[lengths(col_) > 1L]
    } else {
      NULL
    }
    col_tree <- cluster_mat_by_rows(
      data,
      values = t(values),
      by = by_cols,
      .distance_method = .distance_method,
      .linkage_method = .linkage_method,
      .p = .p
    )
    data[[cols_to]] <- vctrs::vec_rep_each(col_tree, key_map$nrows) |>
      vctrs::vec_slice(key_map$mat_to_df)
  }

  data
}
