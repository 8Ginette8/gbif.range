### ==================================================================
### cv_range
### ==================================================================
#' Evaluate a Range Map by Cross-Validation
#' 
#' Rebuild a \code{get_range()} model repeatedly from subsets of the occurrence
#' data stored in a \code{getRange} object and evaluate each rebuild against
#' held-out observations.
#' 
#' @param range_object Object of class \code{getRange}, typically returned by
#' \code{get_range()}.
#' @param cv Character string specifying the cross-validation strategy:
#' \code{"random-cv"} or \code{"block-cv"}.
#' @param nfolds Numeric number of folds.
#' @param nblocks Numeric multiplier used when \code{cv = "block-cv"} to define
#' the total number of spatial blocks as \code{nfolds * nblocks}.
#' @param backpoints Numeric number of regularly spaced background points used
#' as pseudo-absences. Default is \code{10000}.
#' @details The function rebuilds the range map \code{nfolds} times. In each
#' iteration, one fold is reserved for evaluation and the remaining folds are
#' used for training.
#'
#' Two strategies are available: random cross-validation and spatial block
#' cross-validation. The latter reduces the influence of spatial autocorrelation
#' by grouping nearby observations before splitting them across folds.
#'
#' Because true absences are generally unavailable, the evaluation uses a
#' regular grid of background points as pseudo-absences and reports precision,
#' sensitivity, specificity, and TSS.
#' @return A data frame with one row per fold plus a \code{Mean} row, and the
#' columns \code{TP}, \code{FA}, \code{TA}, \code{FP}, \code{Precision},
#' \code{Sensitivity}, \code{Specificity}, and \code{TSS}.
#' @references
#' Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera‐
#' Arroita, G., ... & Dormann, C. F. (2017). Cross‐validation strategies
#' for data with temporal, spatial, hierarchical, or phylogenetic structure.
#' Ecography, 40(8), 913-929.
#' 
#' Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P.,
#' & Thuiller, W. (2021). Novel methods to correct for observer and sampling
#' bias in presence‐only species distribution models. Global Ecology and
#' Biogeography, 30(11), 2312-2325.
#' @example inst/examples/cv_range_help.R
#' @importFrom terra ext extract
#' @export
cv_range <- function(range_object = NULL,
                     cv = 'random-cv',
                     nfolds = 5,
                     nblocks = 2,
                     backpoints = 1e4){

  ######################################################
  ### Stop messages
  ######################################################


  # Main
  if (!class(range_object)%in%"getRange") {
    stop("Given 'range_object' must be a 'getRange' object...")
  }
  if (!cv %in% c("random-cv","block-cv")) {
    stop("Given 'cv' must be 'random-cv' or 'block-cv'...")
  }
  check_numeric(nfolds, "nfolds")
  check_numeric(nblocks, "nblocks")
  check_numeric(backpoints, "backpoints")


  ######################################################
  ### Code
  ######################################################


  # First remove observations considered outliers in get_range()
  # (outside range extent)
  xy.df <- range_object$init.args$occ_coord
  r.ext <- terra::ext(range_object$rangeOutput)
  Xrm.cond <- xy.df$decimalLongitude >= r.ext[1] &
                  xy.df$decimalLongitude <= r.ext[2]
  Yrm.cond <-  xy.df$decimalLatitude >= r.ext[3] &
                  xy.df$decimalLatitude <= r.ext[4]
  xy.df <- xy.df[Xrm.cond & Yrm.cond, ]

  # Samples n regular background points over the original range extent
  x.interv <- (r.ext[2] - r.ext[1]) / (sqrt(backpoints) - 1)
  y.interv <- (r.ext[4] - r.ext[3]) / (sqrt(backpoints) - 1)
  lx <- seq(r.ext[1], r.ext[2], x.interv)
  ly <- seq(r.ext[3], r.ext[4], y.interv)
  bp.xy <- expand.grid(decimalLongitude = lx,
                       decimalLatitude = ly)

  # Combine observations with background
  obs.xy <- xy.df[, c("decimalLongitude", "decimalLatitude")]
  all.xy <- rbind(obs.xy, bp.xy)
  all.xy$Pres <- 0
  all.xy[seq_len(nrow(obs.xy)), "Pres"] <- 1

  # Create cv strata
  if (cv %in% "block-cv") {

    # Normal strata vector
    xy.pres <- all.xy$Pres

    # Inputs for bcv strata
    if (nrow(obs.xy) > 1e4) {
      xy.samp <- sample(x = seq_len(nrow(obs.xy)), size = 1e4)
      xy.pres[seq_len(nrow(obs.xy))][-xy.samp] <- 0
    }

    # Run blocking function
    cv.strat <- make_blocks(
      nfolds = nfolds,
      df = all.xy[,c("decimalLongitude","decimalLatitude")],
      nblocks = nfolds * nblocks,
      pres = xy.pres
    )

  # For random-cv
  } else {
    cv.strat <- make_blocks(
      nfolds = nfolds,
      npoints = nrow(all.xy)
    )
  }

  # Prepare the evaluation df
  cv.df <- as.data.frame(matrix(NA, nfolds + 1, 8))
  row.names(cv.df) <- c(sprintf("CV%d", 1:nfolds),"Mean")
  colnames(cv.df) <- c("TP", "FA", "TA", "FP",
                        "Precision", "Sensitivity", "Specificity", "TSS")

  # Run nfolds time the get_range function + evaluation
  for (i in 1:nfolds)
  {
    cat("...fold", i, sep="")

    # Extract all but %nfolds
    xy.fit <- all.xy[!cv.strat %in% i,]
    xy.eval <- all.xy[cv.strat %in% i,]

    # Run get_range
    range.cv <- get_range(
      occ_coord = xy.fit[xy.fit$Pres %in% 1,],
      ecoreg = range_object$init.args$ecoreg,
      ecoreg_name = range_object$init.args$ecoreg_name,
      degrees_outlier = range_object$init.args$degrees_outlier,
      clust_pts_outlier = range_object$init.args$clust_pts_outlier,
      buff_width_point = range_object$init.args$buff_width_point,
      buff_incrmt_pts_line = range_object$init.args$buff_incrmt_pts_line,
      buff_width_polygon = range_object$init.args$buff_width_polygon,
      dir_temp = range_object$init.args$dir_temp,
      raster = TRUE,
      res = range_object$init.args$res,
      verbose = FALSE
    )
    names(range.cv$rangeOutput) <- "layer"

    # Create data to calculate TP,FA,TA & FP
    Pred <- terra::extract(
      x = range.cv$rangeOutput,
      y = as.data.frame(xy.eval[ ,c("decimalLongitude", "decimalLatitude")])
    )
    xy.eval$Pred <- Pred$layer
    xy.eval[is.na(xy.eval$Pred), "Pred"] <- 0

    # Generate evaluation metrics
    cv.df[i, "TP"] <-
      as.integer(length(which(xy.eval$Pres %in% 1 & xy.eval$Pred %in% 1)))
    cv.df[i, "FA"] <-
      as.integer(length(which(xy.eval$Pres %in% 1 & xy.eval$Pred %in% 0)))
    cv.df[i, "TA"] <-
      as.integer(length(which(xy.eval$Pres %in% 0 & xy.eval$Pred %in% 0)))
    cv.df[i, "FP"] <-
      as.integer(length(which(xy.eval$Pres %in% 0 & xy.eval$Pred %in% 1)))
    cv.df[i, "Precision"] <-
      cv.df[i, "TP"] / (cv.df[i, "TP"] + cv.df[i, "FP"])
    cv.df[i, "Sensitivity"] <-
      cv.df[i, "TP"] / (cv.df[i, "TP"] + cv.df[i, "FA"])
    cv.df[i, "Specificity"] <-
      cv.df[i, "TA"] / (cv.df[i, "TA"] + cv.df[i, "FP"])
    cv.df[i, "TSS"] <-
      cv.df[i, "Sensitivity"] + cv.df[i, "Specificity"] - 1
  }

  cat("","\n")
  
  # Finalize average
  cv.df[nfolds+1, ] <- apply(cv.df[1:nfolds, ], 2, mean, na.rm = TRUE)
  return(cv.df)
}
