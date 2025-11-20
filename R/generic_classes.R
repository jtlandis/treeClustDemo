#' Supported Classes
#' @name supported_classes
#' Possible classes to dispatch on
NULL

#' @rdname supported_classes
#' @export
GR <- methods::getClass("GenomicRanges")

#' @rdname supported_classes
#' @export
GRL <- methods::getClass("GenomicRangesList")

#' @rdname supported_classes
#' @export
SE <- methods::getClass("SummarizedExperiment")

#' @rdname supported_classes
#' @export
RSE <- methods::getClass("RangedSummarizedExperiment")


#' @rdname suppoted_classes
#' @export
class_simple_range <- S7::new_S3_class("simple_range")

class_nested_range <- S7::new_S3_class("nested_range")

s4_simple_range <- methods::setOldClass("simple_range")


new_simple_range <- function(start = integer(), stop = integer(), names = NULL) {
  vctrs::new_rcrd(
    Filter(Negate(is.null), list(
      start = start,
      stop = stop,
      NAMES = names
    )),
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
