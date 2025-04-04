### =========================================================================
### cv blocks
### =========================================================================
#' Block-wise split data into training and testing
#'
#' Creates a stratum vector based on a data.frame with n columns. If the data.frame
#' has one column, strata vector is created based on clusters separated by quantiles. If
#' the data.frame has two or more columns, strata vector is created based on 'Clustering
#' Large Applications' (function 'clara' from package cluster). Also, instead of a data.frame
#' the argument 'npoints' can be provided, which create groups by random sampling. An
#' opitimization algorithm (function 'gridSearch' from package NMOF) optimizes for equal stratum sizes.
#'
#' @param nfolds Numeric. Number of approximately equal-sized classes (folds) to separate groups in block-cross
#' validation.
#' @param df Object of class 'data.frame' with n columns containing critera for cluster building.
#' Not necessary if argument npoints is supplied.
#' @param nblocks Number of clusters (blocks) based on the number of folds that should be built.
#' Minimum is the same number as 'nFolds'. Maximum is nrow(df)/10.
#' @param pres Binary vector. Optional argument. If 'df' is supplied, this argument can be used to
#' save processing time. '1' stands for the points on which CLARA is appplied, and '0' stands for
#' the points on which K-nearest neighbors is applied relative to the '1'. If 'df' is not supplied,
#' for which points should random sampling be made? 
#' @param npoints Optional argument if 'df' is not supplied. For how many points should random sampling be made?
#' @return Object of class 'vector' of length nrow(df) or 'npoints', with integers defining
#' different folds.
#' @references
#' Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O., Wang, Z., & Zimmermann, N. E. (2020).
#' Model complexity affects species distribution projections under climate change.
#' Journal of Biogeography, 47(1), 130-142.
#' @example inst/examples/make_blocks_help.R
#' @importFrom cluster clara
#' @importFrom class knn
#' @importFrom NMOF gridSearch
#' @importFrom stats quantile
#' @export
make_blocks <- function(nfolds = 5,
                      df = data.frame(),
                      nblocks = nfolds*2,
                      npoints = NA,
                      pres = numeric()){

  ### ------------------------
  ### check input data
  ### ------------------------

  if (nrow(df) == 0 & is.na(npoints)){
    stop("Please supply number of points if no data.frame is supplied")
  }

  ### ------------------------
  ### generate clusters
  ### ------------------------

  if (nrow(df) == 0){

    ### do ordinary sampling if no strata are supplied
    out.strat <- sample(rep(1:nfolds, ceiling(npoints / nfolds)), size = npoints)

  } else {

    # check for reasonable number of boxes
    if (nrow(df) < 4*nblocks){
      stop("Too many boxes required!")
    }

    if (ncol(df) == 1){

      ### do quantile-based clustering if df contains only one column
      rngi <- stats::quantile(df[, 1], probs = 0:(nblocks) / (nblocks))
      rngi[1] <- rngi[1] - 1
      rngi[length(rngi)] <- rngi[length(rngi)] + 1
      clist <- as.numeric(cut(df[, 1], breaks = rngi, right = TRUE))

    } else {

      # Scale input data
      scd <- apply(df, 2, scale)

      if (length(pres) == 0){

        # do kmedoid clustering for 2 or more columns in df
        kmed <- cluster::clara(scd, k = nblocks, metric = "euclidean")

        # get clusters
        clist <- kmed$clustering

      } else {

        # do kmedoid clustering for 2 or more columns in df
        kmed <- cluster::clara(scd[which(pres == 1), ], k = nblocks, metric = "manhattan")
        knnab <- class::knn(train = scd[which(pres == 1), ], test = scd[which(pres == 0), ], cl = kmed$clustering)
        clist <- kmed$clustering
        cliful <- rep(NA,nrow(scd))
        cliful[which(pres == 1)] <- kmed$clustering
        cliful[which(pres == 0)] <- as.numeric(knnab)
      }
    }

    # sort obtained clusters
    tbl <- sort(table(clist), decreasing = TRUE)

    ### ------------------------
    ### regularly assign clusters to strata
    ### ------------------------

    if (nblocks != nfolds){

      # prepare strata layers
      grps <- rep(list(numeric()), nfolds)

      # for the clusters with many observations
      # distribute them regularly among strata but keep
      # last six clusters for estimating most
      # regular distribution

      if (length(tbl) > 6){

        for (i in 1:(length(tbl) - 6)){

          # determine to which stratum the cluster should be
          # added
          fl <- (floor((i - 1) / nfolds))
          if (fl %% 2 == 0){

            j <- round(1 + nfolds * ((i - 1) / nfolds - (floor((i - 1) / nfolds))))

          } else {

            j <- (nfolds + 1) - round(1 + nfolds * ((i - 1) / nfolds - (floor((i - 1) / nfolds))))
          }
          # add cluster
          grps[[j]] <- append(grps[[j]], tbl[i])
        }
      }

      # prepare for optimal distribution of last 6 clusters
      vlis <- factor(1:nfolds, levels = 1:nfolds)
      prs <- rep(list(vlis), min(length(tbl), 6))
      sstab <- tbl[max(1, (length(tbl) - 5)):length(tbl)]

      # Run brute-forcing gridSearch obtimization
      srch <- NMOF::gridSearch(levels = prs, fun = optme, nms = as.vector(sstab), grps = grps, tot = sum(tbl))

      # pull out results
      wi <- as.numeric(as.character(srch$minlevels))

      # combine results with predistributed clusters
      for (i in 1:length(grps)){
        grps[[i]] <- append(grps[[i]], sstab[wi == i])
      }

      # define vector with output strata
      out.strat <- rep(NA, nrow(df))
      for (i in 1:length(grps)){

        if (length(pres) == 0){

          out.strat[which(as.character(clist) %in% names(grps[[i]]))] = i

        } else {

          out.strat[which(as.character(cliful) %in% names(grps[[i]]))] = i
        }
      }

    } else {

      # if as many strata as clusters are required, simply return clusters
      if (length(pres) == 0){

        out.strat <- clist
      } else{

        out.strat <- cliful
      }
    }
  }
  # return result
  return(out.strat)
}
