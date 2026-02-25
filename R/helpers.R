#' Check if an argument is NULL or NA
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_null_na <- function(x, name) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) {
    stop(
      sprintf(
        "Given '%s' is NA or NULL, a value must be provided...", name
      )
    )
  }
}

#' Check if argument is a non-missing logical scalar
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_logical <- function(x, name) {
  check_null_na(x, name)
  if (!is.logical(x) || length(x) != 1) {
    stop(
      sprintf(
        "Given '%s' must be a single logical TRUE or FALSE...", name
      )
    )
  }
}

#' Check if argument is a non-missing numeric scalar
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_numeric <- function(x, name) {
  check_null_na(x, name)
  if (!is.numeric(x) || length(x) != 1) {
    stop(
      sprintf(
        "Given '%s' must be a single numeric value...", name
      )
    )
  }
}

#' Check if argument is a non-empty character vector (no NA)
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @keywords internal
check_character_vector <- function(x, name) {
  check_null_na(x, name)
  if (!is.character(x) || length(x) == 0 || any(is.na(x))) {
    stop(
      sprintf(
        "Given '%s' must be a non-empty character vector without NA...", name
      )
    )
  }
}

#' Check if argument is a numeric vector of exact length
#'
#' @param x Argument value
#' @param name Argument name for error reporting
#' @param len Expected length of numeric vector
#' @keywords internal
check_numeric_range <- function(x, name, len) {
  check_null_na(x, name)
  if (!is.numeric(x) || length(x) != len || any(is.na(x))) {
    stop(
      sprintf(
        "Given '%s' must be a numeric vector of length %d without NA...",
        name,
        len
      )
    )
  }
}

#' Summary log helper for get_gbif filtering
#'
#' @param log Argument value
#' @param step_name Filtering step name
#' @param before Pre-filtering number of observations
#' @param after Post-filtering number of observations
#' @keywords internal
log_step <- function(log, step_name, before, after) {
  rbind(
    log,
    data.frame(
      step = step_name,
      removed = before - after,
      remaining = after,
      stringsAsFactors = FALSE
    )
  )
}

#' Original legend bar (from Dr. Philipp Brun)
#'
#' @param colors Parameter 1
#' @param crds Parameter 2
#' @param horiz Parameter 3
#' @param zrng Parameter 5
#' @param at Parameter 6
#' @param labs Parameter 7
#' @param tickle Parameter 8
#' @param title Parameter 9
#' @param lablag Parameter 10
#' @param titlag Parameter 11
#' @param box Parameter 12
#' @param breaks Parameter 13
#' @param cx Parameter 14
#' @param tria Parameter 15
#' @importFrom graphics polygon lines text
#' @keywords internal
cscl <- function (colors, crds, horiz = FALSE, zrng = c(0, 100), at = 10 *
    0:10, labs = NA, tickle = 0.2, title = 1, lablag = 1, titlag = 2,
    box = TRUE, breaks, cx = 0.8, tria = "n")
{
    rel.lab <- (at - zrng[1])/(zrng[2] - zrng[1])
    if (length(labs) == 1 && is.na(labs)) {
        labs <- at
    }
    if (horiz == TRUE) {
        if (missing(breaks)) {
            brks <- seq(
              from = crds[1],
              to = crds[2],
              length.out = length(colors) + 1
            )
        }
        else if (length(breaks) != (length(colors) + 1)) {
            stop("breaks should be a vector of length colors +1")
        }
        else {
            brks <- sort(breaks) * (crds[2] - crds[1])/(max(breaks) -
                min(breaks)) + crds[1]
        }
        for (i in 1:(length(colors))) {
            polygon(c(brks[i], brks[i + 1], brks[i + 1], brks[i],
                brks[i]), c(crds[3], crds[3], crds[4], crds[4],
                crds[3]), col = colors[i], border = F)
        }
        if (tria %in% c("b", "l")) {
            tipy <- crds[1] - 0.05 * abs(crds[2] - crds[1])
            polygon(c(crds[c(1, 1)], tipy, crds[1]), c(crds[3],
                crds[4], mean(crds[3:4]), crds[3]), col = colors[1])
        }
        if (tria %in% c("b", "u")) {
            tipy <- crds[2] + 0.05 * abs(crds[2] - crds[1])
            polygon(c(crds[c(2, 2)], tipy, crds[2]), c(crds[3],
                crds[4], mean(crds[3:4]), crds[3]), col = rev(colors)[1])
        }
        x.lab <- rel.lab * (crds[2] - crds[1]) + crds[1]
        for (i in 1:length(x.lab)) {
            lines(c(x.lab[i], x.lab[i]), c(crds[3], crds[3] -
                tickle * abs(crds[4] - crds[3])))
        }
        text(labels = labs, x = x.lab, y = rep(crds[3] - lablag *
            abs(crds[4] - crds[3]), 10), cex = cx)
        text(labels = title, x = mean(crds[1:2]), y = crds[3] -
            titlag * abs(crds[4] - crds[3]), cex = cx)
    }
    else if (horiz == FALSE) {
        if (is.na(tickle)) {
            tickle <- 0.2 * (crds[2] - crds[1])
        }
        if (missing(breaks)) {
            brks <- seq(
              from = crds[4],
              to = crds[3],
              length.out = length(colors) + 1
            )
        }
        else if (length(breaks) != (length(colors) + 1)) {
            stop("breaks should be a vector of length colors +1")
        }
        else {
            brks <- sort(breaks) * (crds[4] - crds[3])/(max(breaks) -
                min(breaks)) + crds[3]
        }
        for (i in 1:(length(colors))) {
            polygon(c(crds[1], crds[1], crds[2], crds[2], crds[1]),
                c(brks[i], brks[i + 1], brks[i + 1], brks[i],
                  brks[i]), col = rev(colors)[i], border = F)
        }
        if (tria %in% c("b", "l")) {
            tipy <- crds[3] - 0.05 * abs(crds[4] - crds[3])
            polygon(c(crds[1], crds[2], mean(crds[1:2]), crds[1]),
                c(crds[c(3, 3)], tipy, crds[3]), col = colors[1])
        }
        if (tria %in% c("b", "t")) {
            tipy <- crds[4] + 0.05 * abs(crds[4] - crds[3])
            polygon(c(crds[1], crds[2], mean(crds[1:2]), crds[1]),
                c(crds[c(4, 4)], tipy, crds[4]), col = rev(colors)[1])
        }
        x.lab <- rev(rel.lab) * (crds[4] - crds[3]) + crds[3]
        for (i in 1:length(x.lab)) {
            lines(c(crds[1], crds[1] - tickle * abs(crds[2] -
                crds[1])), c(x.lab[i], x.lab[i]))
        }
        text(labels = rev(labs), y = x.lab, x = rep(crds[1] -
            lablag * abs(crds[2] - crds[1]), 10), cex = cx, pos = 2)
        text(labels = title, y = mean(crds[3:4]), x = crds[1] -
            titlag * abs(crds[2] - crds[1]), cex = cx, srt = 90)
    }
    if (box) {
        polygon(c(crds[1], crds[2], crds[2], crds[1], crds[1]),
            c(crds[3], crds[3], crds[4], crds[4], crds[3]))
    }
}

#' Custom figure label plot
#'
#' @param text Parameter 1
#' @param region Parameter 2
#' @param pos Parameter 3
#' @param cex Parameter 5
#' @param margin Parameter 6
#' @param ... Parameter 7
#' @importFrom grDevices dev.size
#' @importFrom graphics grconvertX grconvertY strheight strwidth
#' @keywords internal
fig_label <- function(
    text, region="figure",
    pos="topleft", cex=NULL, margin=0.02, ...) {
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    if(region == "figure") {
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  x_margin <- (x[2] - x[1]) * margin
  y_margin <- (y[2] - y[1]) * margin
  x1 <- switch(pos,
    topleft     = x[1] + sw + x_margin, 
    left        = x[1] + sw + x_margin,
    bottomleft  = x[1] + sw + x_margin,
    top         = (x[1] + x[2]) / 2,
    center      = (x[1] + x[2]) / 2,
    bottom      = (x[1] + x[2]) / 2,
    topright    = x[2] - sw - x_margin,
    right       = x[2] - sw - x_margin,
    bottomright = x[2] - sw - x_margin)
  y1 <- switch(pos,
    topleft     = y[2] - sh - y_margin,
    top         = y[2] - sh - y_margin,
    topright    = y[2] - sh - y_margin,
    left        = (y[1] + y[2]) / 2,
    center      = (y[1] + y[2]) / 2,
    right       = (y[1] + y[2]) / 2,
    bottomleft  = y[1] + sh + y_margin,
    bottom      = y[1] + sh + y_margin,
    bottomright = y[1] + sh + y_margin)
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x, y)))
}
