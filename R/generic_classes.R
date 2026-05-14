#' Supported Classes
#' @name supported_classes
#' @description Possible classes to dispatch on
NULL

#' @rdname supported_classes
#' @export
GR <- methods::getClass("GenomicRanges", where = "GenomicRanges")

#' @rdname supported_classes
#' @export
GRL <- methods::getClass("GenomicRangesList", where = "GenomicRanges")

#' @rdname supported_classes
#' @export
SE <- methods::getClass("SummarizedExperiment", where = "SummarizedExperiment")

#' @rdname supported_classes
#' @export
RSE <- methods::getClass(
  "RangedSummarizedExperiment",
  where = "SummarizedExperiment"
)


#' @rdname suppoted_classes
#' @export
class_simple_range <- S7::new_S3_class("simple_range")

class_nested_range <- S7::new_S3_class("nested_range")

s4_simple_range <- methods::setOldClass("simple_range")


new_simple_range <- function(
  start = integer(),
  stop = integer(),
  names = NULL
) {
  vctrs::new_rcrd(
    Filter(
      Negate(is.null),
      list(
        start = start,
        stop = stop,
        NAMES = names
      )
    ),
    class = "simple_range"
  )
}

#' a light weight way to constuct a "Range" class
#' when we only care about the start and stop positions
#' since splitting a GRangesList can be fairly expensive
#' @export
simple_range <- function(x) UseMethod("simple_range")

#' @export
simple_range.GRanges <- function(x) {
  out <- new_simple_range(
    start = IRanges::start(x),
    stop = IRanges::end(x),
    names = names(x)
  )
  out
}

#' @export
names.simple_range <- function(x) {
  if ("NAMES" %in% vctrs::fields(x)) {
    vctrs::field(x, "NAMES")
  } else {
    NULL
  }
}

#' @export
`names<-.simple_range` <- function(x, value) {
  if (is.null(value)) {
    new_simple_range(
      vctrs::field(x, "start"),
      vctrs::field(x, "stop")
    )
  } else {
    value <- vctrs::vec_cast(value, character())
    x <- new_simple_range(
      vctrs::field(x, "start"),
      vctrs::field(x, "stop"),
      character(vctrs::vec_size(x))
    )
    vctrs::field(x, "NAMES") <- value
    x
  }
}

#' @export
format.simple_range <- function(x, ...) {
  sprintf("%i-%i", vctrs::field(x, "start"), vctrs::field(x, "stop"))
}

#' @export
simple_range.GRangesList <- function(x) {
  range <- unlist(x)
  data <- simple_range(range)
  out <- vctrs::new_vctr(
    vctrs::vec_chop(data, sizes = lengths(x)),
    class = "nested_range"
  )
  names(out) <- names(x)
  out
}

#' @export
format.nested_range <- function(x, ...) {
  sprintf(
    "%s",
    seq_along(x),
    vapply(
      vctrs::vec_data(x),
      \(x) paste(format(utils::head(x)), collapse = ", "),
      character(1)
    )
  )
}

#' @importFrom IRanges start
#' @importFrom IRanges end
NULL

#' @export
methods::setMethod(
  "start",
  s4_simple_range,
  function(x, ...) {
    vctrs::field(x, "start")
  }
)

#' @export
methods::setMethod(
  "end",
  s4_simple_range,
  function(x, ...) {
    vctrs::field(x, "stop")
  }
)

S7::method(
  list_unchop,
  list(methods::getClass("CompressedList"), NULL)
) <- function(x, ptype, ..., indices = NULL) {
  unlistData <- unlist(x)
  data_class <- class(unlistData)
  slice_ <- methods::selectMethod(
    "[",
    c(x = data_class, i = "integer", j = "missing", drop = "missing")
  )

  if (!is.null(indices)) {
    lens_x <- lengths(x)
    lens_i <- lengths(indices)
    if (!all(res <- lens_x == lens_i | lens_x == 1)) {
      failure <- which(!res)[1]
      stop(
        sprintf(
          "Can't recycle `x[[%i]]` (size %i) to size %i",
          failure,
          lens_x[failure],
          lens_i[failure]
        )
      )
    }
    unlistData <- slice_(unlistData, vctrs::list_unchop(indices))
  }
  unlistData
}


S7::method(
  list_unchop,
  list(S7::class_list, methods::getClass("CompressedList"))
) <- function(x, ptype, ..., indices = NULL) {
  unlistData <- do.call("c", lapply(x, unlist))
  widths <- lapply(x, lengths)
  ends <- cumsum(vctrs::list_unchop(widths))
  if (!is.null(indices)) {
    part <- IRanges::PartitioningByEnd(ends)
    start <- IRanges::start(part)
    end <- IRanges::end(part)
    ends <- vctrs::list_unchop(widths, indices = indices) |>
      cumsum()
    ord <- order(vctrs::list_unchop(indices))
    slice_i <- vctrs::list_unchop(Map(\(s, e) s:e, s = start, e = end)[ord])
    unlistData <- unlistData[slice_i]
  }

  relist(unlistData, skeleton = IRanges::PartitioningByEnd(ends))
}
