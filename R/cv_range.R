### ==================================================================
### cv_range
### ==================================================================
#' Evaluate the performance of a range map with cross-validation
#' 
#' Assesses the accuracy of a species range map by applying cross-validation using the observations
#' (restricted to the extent) and function arguments of a get_range() object. By using these initial parameters,
#' the function iteratively re-generates the range map n times, each time using a defined percentage of the extent
#' observations for training, while evaluating the quality of the range map using the remaining records
#' (by default, nfolds = 5 -> calibration = 80%, evaluation = 20%). Two cross-validation methods are available:
#' random and spatial block cross-validation. In random cross-validation, a random subset of the observations
#' is chosen for training in each fold, with the generated map evaluated on the remaining data. In spatial block
#' cross-validation, the observations are spatially divided into blocks, and each fold uses a different set
#' of blocks for training and testing, ensuring that spatial dependencies are properly considered in the
#' evaluation process. Evaluation metrics here employed are Precision, Sensitivity, Specificity and TSS.
#' It is important to note that since no absences are available for evaluation, a uniform random layer of
#' background points is first generated over the study area extent and used as pseudo-absences proxy.
#' 
#' @param range_object A get_range() object.
#' @param cv Character. Should the range map be evaluated with random ('random-cv') or spatial block
#' cross validation ('block-cv').
#' @param nfolds Numeric. Number of chosen folds for cross-validation.
#' @param nblocks Numeric. Only applies if 'block-cv' is employed. Defined the number of blocks per fold.
#' @param backpoints Numeric (optional). Number of regular background points that should be sampled. Defaut is 10,000.
#' @return A data.frame with 'nfolds' rows and 8 evaluation columns.
#' @references
#' Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera‐Arroita, G., ... & Dormann, C. F.
#' (2017). Cross‐validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure.
#' Ecography, 40(8), 913-929.
#' 
#' Chauvier, Y., Zimmermann, N. E., Poggiato, G., Bystrova, D., Brun, P., & Thuiller, W. (2021).
#' Novel methods to correct for observer and sampling bias in presence‐only species distribution models.
#' Global Ecology and Biogeography, 30(11), 2312-2325.
#' @example inst/examples/cv_range_help.R
#' @importFrom terra ext extract
#' @export
cv_range <- function(range_object = NULL,
                     cv = 'random-cv',
                     nfolds = 5,
                     nblocks = 2,
                     backpoints = 1e4){

  # Firt remove observations considered as outliers in get_range (outside range extent)
  xy.df <- range_object[[1]]$occ_coord
  r.ext <- terra::ext(range_object[[2]])
  Xrm.cond <- xy.df$decimalLongitude >= r.ext[1] & xy.df$decimalLongitude <= r.ext[2]
  Yrm.cond <-  xy.df$decimalLatitude >= r.ext[3] & xy.df$decimalLatitude <= r.ext[4]
  xy.df <- xy.df[Xrm.cond & Yrm.cond, ]

  # Samples n regular background points over the original range extent
  x.interv <- (r.ext[2] - r.ext[1]) / (sqrt(backpoints) - 1)
  y.interv <- (r.ext[4] - r.ext[3]) / (sqrt(backpoints) - 1)
  lx <- seq(r.ext[1], r.ext[2], x.interv)
  ly <- seq(r.ext[3], r.ext[4], y.interv)
  bp.xy <- expand.grid(decimalLongitude = lx, decimalLatitude = ly)

  # Combine observations with background
  obs.xy <- xy.df[, c("decimalLongitude", "decimalLatitude")]
  all.xy <- rbind(obs.xy, bp.xy)
  all.xy$Pres <- 0
  all.xy[1:nrow(obs.xy), "Pres"] <- 1

  # Create cv strata
  if (cv %in% "block-cv") {

    # Normal strata vector
    xy.pres <- all.xy$Pres

    # Inputs for bcv strata
    if (nrow(obs.xy) > 1e4) {
      xy.samp <- sample(x = 1:nrow(obs.xy), size = 1e4)
      xy.pres[1:nrow(obs.xy)][-xy.samp] <- 0
    }

    # Run blocking function
    cv.strat <- make_blocks(nfolds = nfolds,
                           df = all.xy[,c("decimalLongitude","decimalLatitude")],
                           nblocks = nfolds * nblocks,
                           pres = xy.pres)

  # For random-cv
  } else {
    cv.strat <- make_blocks(nfolds = nfolds,
                           npoints = nrow(all.xy))
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
    range.cv <- get_range(occ_coord = xy.fit[xy.fit$Pres %in% 1,],
                  bioreg = range_object[[1]]$bioreg,
                  bioreg_name = range_object[[1]]$bioreg_name,
                  degrees_outlier = range_object[[1]]$degrees_outlier,
                  clustered_points_outlier = range_object[[1]]$clustered_points_outlier,
                  buffer_width_point = range_object[[1]]$buffer_width_point,
                  buffer_increment_point_line = range_object[[1]]$buffer_increment_point_line,
                  buffer_width_polygon = range_object[[1]]$buffer_width_polygon,
                  dir_temp = range_object[[1]]$dir_temp,
                  raster = TRUE,
                  res = range_object[[1]]$res,
                  verbose = FALSE)

    # Create data to calculate TP,FA,TA & FP
    Pred <- terra::extract(range.cv[[2]], xy.eval[ ,c("decimalLongitude", "decimalLatitude")])
    xy.eval$Pred <- Pred$layer
    xy.eval[is.na(xy.eval$Pred), "Pred"] <- 0

    # Generate evaluation metrics
    cv.df[i, "TP"] <- as.integer(length(which(xy.eval$Pres %in% 1 & xy.eval$Pred %in% 1)))
    cv.df[i, "FA"] <- as.integer(length(which(xy.eval$Pres %in% 1 & xy.eval$Pred %in% 0)))
    cv.df[i, "TA"] <- as.integer(length(which(xy.eval$Pres %in% 0 & xy.eval$Pred %in% 0)))
    cv.df[i, "FP"] <- as.integer(length(which(xy.eval$Pres %in% 0 & xy.eval$Pred %in% 1)))
    cv.df[i, "Precision"] <- cv.df[i, "TP"] / (cv.df[i, "TP"] + cv.df[i, "FP"])
    cv.df[i, "Sensitivity"] <- cv.df[i, "TP"] / (cv.df[i, "TP"] + cv.df[i, "FA"])
    cv.df[i, "Specificity"] <- cv.df[i, "TA"] / (cv.df[i, "TA"] + cv.df[i, "FP"])
    cv.df[i, "TSS"] <- cv.df[i, "Sensitivity"] + cv.df[i, "Specificity"] - 1
  }

  cat("","\n")
  
  # Finalize average
  cv.df[nfolds+1, ] <- apply(cv.df[1:nfolds, ], 2, mean, na.rm = TRUE)
  return(cv.df)
}
