#' Check Whether an Argument Is Missing
#'
#' @param x Argument value.
#' @param name Argument name used in error messages.
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

#' Check Whether an Argument Is a Logical Scalar
#'
#' @param x Argument value.
#' @param name Argument name used in error messages.
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

#' Check Whether an Argument Is a Numeric Scalar
#'
#' @param x Argument value.
#' @param name Argument name used in error messages.
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

#' Check Whether an Argument Is a Character Vector
#'
#' @param x Argument value.
#' @param name Argument name used in error messages.
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

#' Check Whether an Argument Is a Numeric Vector of Fixed Length
#'
#' @param x Argument value.
#' @param name Argument name used in error messages.
#' @param len Expected length of the numeric vector.
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

#' Append a Filtering Step to a Log Table
#'
#' @param log Existing log data frame.
#' @param step_name Name of the filtering step.
#' @param before Number of records before filtering.
#' @param after Number of records after filtering.
#' @return An updated log data frame.
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

#' Draw a Custom Legend Bar
#'
#' Internal plotting helper adapted from Philipp Brun.
#'
#' @param colors Vector of fill colors.
#' @param crds Numeric vector of length four giving the plotting coordinates of
#' the legend box.
#' @param horiz Logical. Should the legend be drawn horizontally?
#' @param zrng Numeric vector of length two giving the value range represented
#' by the legend.
#' @param at Numeric values at which tick marks should be drawn.
#' @param labs Optional labels for \code{at}. If \code{NA}, \code{at} is used.
#' @param tickle Numeric tick-mark length scaling factor.
#' @param title Legend title.
#' @param lablag Numeric offset multiplier for tick labels.
#' @param titlag Numeric offset multiplier for the title.
#' @param box Logical. Should a border be drawn around the legend?
#' @param breaks Optional custom break points. Must have length
#' \code{length(colors) + 1} if supplied.
#' @param cx Numeric text expansion factor.
#' @param tria Character flag controlling triangular tips at the legend ends.
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

#' Draw a Figure Label at a Standard Position
#'
#' @param text Label text.
#' @param region Character string specifying whether the label should be placed
#' relative to the \code{"figure"}, \code{"plot"}, or \code{"device"}.
#' @param pos Character string giving the label position.
#' @param cex Optional text expansion factor.
#' @param margin Numeric margin used to offset the label from the selected
#' boundary.
#' @param ... Additional arguments passed to \code{graphics::text()}.
#' @return Invisibly returns the plotting coordinates used to position the
#' label.
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
