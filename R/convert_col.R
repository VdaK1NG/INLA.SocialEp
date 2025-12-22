#' Generate stable color labels
#'
#' This function generates stable color labels to ease comparison when representing several figures that
#' we would want to have the same palette.
#'
#' @param data Vector of numeric values to be convert into colors
#' @param breaks Vector of breaks to use for cutting
#' @param pal_fun Vector of colors for the palette
#' @param right Right option for the cut function
#' @param include.lowest Include.lowest for the cut function
#' @param na.col Colour to use in case of NAs
#' @return Returns a list including the colors generated (fill_by), the tags for the intervals (tags) and the colors corresponding for each interval (colors)
#' @export

convert_col <- function(data,
                        breaks,
                        pal_fun,
                        include.lowest=TRUE,
                        right = TRUE,
                        na.col = NA) {

  if (!is.numeric(breaks) || length(breaks) < 2)
    stop("`breaks` must be a numeric vector with length >= 2.")

  nbins <- length(breaks) - 1
  cols  <- pal_fun(nbins)
  btxt   <- format(breaks, trim = TRUE, scientific = FALSE)
  left   <- head(btxt, -1)
  rightb <- tail(btxt, -1)

  labs <- if (right) {
    paste0("(", left, ",", rightb, "]")
  } else {
    paste0("[", left, ",", rightb, ")")
  }
  labs <- str_replace_all(labs, "(?<=\\d)\\.0(?=[,\\]])", "")

  data_lab <- if (is.numeric(data)) {
    as.character(base::cut(data, breaks = breaks, right = right, include.lowest = include.lowest))
  } else {
    trimws(as.character(data))
  }

  out <- cols[match(data_lab, labs)]
  out[is.na(data_lab)] <- na.col
  data <- list()
  data[[1]] <- out
  data[[2]] <- labs
  data[[3]] <- cols
  names(data) <- c("fill_by", "tags", "colors")
  data
}

