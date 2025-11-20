#' @export
intersect_bp <- function(start, stop, ref_start, ref_stop) {
  # browser()
  # ref:              [----]
  # others:   [----]              # not compared
  #               [-----]
  #                   [----]
  #                      [-----]
  to_compare <- (ref_start <= stop & ref_stop >= start)
  #| (ref_start > start & ref_stop > stop)
  #  only look at positions that will overlap
  # res <- numeric(length(start))
  # uni -> union
  # int -> intersection
  to_compare <- which(to_compare)
  uni_start <- int_start <- start[to_compare]
  # ref:              [----]
  # int_start:
  #               [-->[
  #                   [
  #                      [
  int_start[int_start < ref_start] <- ref_start
  uni_stop <- int_stop <- stop[to_compare]
  # ref:              [----]
  # int_stop
  #                     ]
  #                        ]
  #                        ]<--]
  int_stop[int_stop > ref_stop] <- ref_stop
  # ref:              [----]
  # int_ranges:
  #                   [-]
  #                   [----]
  #                      [-]

  # ref:              [----]
  # uni_start:
  #               [
  #                   [
  #                   [<-[
  uni_start[uni_start > ref_start] <- ref_start
  # ref:              [----]
  # uni_stop:
  #                     ]->]
  #                        ]
  #                            ]
  uni_stop[uni_stop < ref_stop] <- ref_stop

  # ref:              [----]
  # uni_ranges:
  #               [--------]
  #                   [----]
  #                   [-------]

  # ref:              [----]
  # res:
  #   comp1:
  #                   [-]
  #               [--------]
  #   comp2:
  #                   [----]
  #                   [----]
  #   comp3:
  #                      [-]
  #                   [-------]

  # res[to_compare] <- (int_stop - int_start + 1L) / (uni_stop - uni_start + 1L)
  # # if(any(res<0))
  # #   browser()
  # # res
  # res
  list(
    index = to_compare,
    intersection = int_stop - int_start + 1L,
    union = uni_stop - uni_start + 1L
  )
}


#' @param start numeric vector of start positions
#' @param stop numeric vector of end positions
#' @param trans some function to transform results off
#' intersections and unions
#' @return a sparseMatrix
#' @export
bp_similarity_ <- function(start, stop, trans) {
  n <- length(start)
  # mat <- matrix(nrow = n, ncol = n)
  #
  is <- vector("list", n)
  js <- vector("list", n)
  xs <- vector("list", n)
  # browser()
  for (i in seq_len(n)) {
    ref_start <- start[i]
    ref_stop <- stop[i]
    seq_i <- i:n
    data <- intersect_bp(
      start[seq_i], stop[seq_i],
      ref_start = ref_start, ref_stop = ref_stop
    )
    vec <- trans(data)
    indx <- data$index
    js[[i]] <- rep_len(i, length(indx))
    xs[[i]] <- vec
    is[[i]] <- seq_i[indx]
  }
  Matrix::sparseMatrix(
    i = unlist(is),
    j = unlist(js),
    x = unlist(xs),
    dims = c(n, n),
    symmetric = TRUE
  )
}

#' @export
bp_simi_percent_ <- function(start, stop) {
  bp_similarity_(
    start, stop,
    function(.data) {
      .data$intersection / .data$union
    }
  )
}

#' @export
bp_simi_width_ <- function(start, stop) {
  bp_similarity_(
    start, stop,
    function(.data) {
      .data$intersection
    }
  )
}

compute_bp_simi_width <- function(object) {
  bp_simi_width_(IRanges::start(object), IRanges::end(object))
}

compute_bp_simi_percent <- function(object) {
  bp_simi_percent_(IRanges::start(object), IRanges::end(object))
}
