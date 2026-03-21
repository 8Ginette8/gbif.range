### =========================================================================
### cv blocks
### =========================================================================
#' Split Data into Approximately Balanced Folds
#'
#' Create a fold-assignment vector for random or spatially structured
#' cross-validation.
#'
#' @param nfolds Numeric number of folds to create.
#' @param df Optional \code{data.frame} whose columns define the structure used
#' to build folds. When omitted, random folds are created from
#' \code{npoints}.
#' @param nblocks Numeric number of initial clusters used to build the folds.
#' Must be at least \code{nfolds}.
#' @param pres Optional binary vector used to restrict CLARA clustering to a
#' subset of rows and assign the remainder by nearest neighbours.
#' @param npoints Optional number of observations to split when \code{df} is
#' not supplied.
#' @details If \code{df} has one column, folds are based on quantile bins. If
#' \code{df} has two or more columns, folds are based on CLARA clustering.
#' Remaining clusters are assigned to folds with an optimization step that
#' tries to balance fold sizes as evenly as possible.
#' @return An integer vector of length \code{nrow(df)} or \code{npoints},
#' giving the fold assignment for each observation.
#' @references
#' Brun, P., Thuiller, W., Chauvier, Y., Pellissier, L., Wüest, R. O.,
#' Wang, Z., & Zimmermann, N. E. (2020). Model complexity affects species
#' distribution projections under climate change. Journal of Biogeography,
#' 47(1), 130-142.
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

  ###########################################
  ### Check input data
  ###########################################


  # N points
  if (nrow(df) == 0 & is.na(npoints)){
    stop("Please supply number of points if no data.frame is supplied")
  }

  # General
  check_numeric(nfolds, "nfolds")
  check_numeric(nblocks, "nblocks")


  ###########################################
  ### Generate clusters
  ###########################################


  if (nrow(df) == 0){

    ### do ordinary sampling if no strata are supplied
    out.strat <- sample(
      rep(seq_len(nfolds),
        ceiling(npoints / nfolds)
      ),
      size = npoints
    )

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
        kmed <- cluster::clara(
          x = scd[which(pres == 1), ],
          k = nblocks,
          metric = "manhattan"
        )
        knnab <- class::knn(
          train = scd[which(pres == 1), ],
          test = scd[which(pres == 0), ],
          cl = kmed$clustering
        )
        clist <- kmed$clustering
        cliful <- rep(NA,nrow(scd))
        cliful[which(pres == 1)] <- kmed$clustering
        cliful[which(pres == 0)] <- as.numeric(knnab)
      }
    }

    # sort obtained clusters
    tbl <- sort(table(clist), decreasing = TRUE)


    ###########################################
    ### Regularly assign clusters to strata
    ###########################################


    if (nblocks != nfolds){

      # prepare strata layers
      grps <- rep(list(numeric()), nfolds)

      # for the clusters with many observations
      # distribute them regularly among strata but keep
      # last six clusters for estimating most
      # regular distribution

      if (length(tbl) > 6){

        for (i in seq_len(length(tbl) - 6)) {

          # determine to which stratum the cluster should be
          # added
          fl <- (floor((i - 1) / nfolds))
          if (fl %% 2 == 0){

            j <- round(
              1 + nfolds * ((i - 1) / nfolds - (floor((i - 1) / nfolds)))
            )

          } else {

            j <- (nfolds + 1) - round(
              1 + nfolds * ((i - 1) / nfolds - (floor((i - 1) / nfolds)))
            )
          }
          # add cluster
          grps[[j]] <- append(grps[[j]], tbl[i])
        }
      }

      # prepare for optimal distribution of last 6 clusters
      vlis <- factor(seq_len(nfolds), levels = seq_len(nfolds))
      prs <- rep(list(vlis), min(length(tbl), 6))
      sstab <- tbl[max(1, (length(tbl) - 5)):length(tbl)]

      # Run brute-forcing gridSearch obtimization
      srch <- NMOF::gridSearch(
        levels = prs,
        fun = optme,
        nms = as.vector(sstab),
        grps = grps,
        tot = sum(tbl)
      )

      # pull out results
      wi <- as.numeric(as.character(srch$minlevels))

      # combine results with predistributed clusters
      for (i in seq_along(grps)){
        grps[[i]] <- append(grps[[i]], sstab[wi == i])
      }

      # define vector with output strata
      out.strat <- rep(NA, nrow(df))
      for (i in seq_along(grps)){

        if (length(pres) == 0){

          out.strat[which(as.character(clist) %in% names(grps[[i]]))] <- i

        } else {

          out.strat[which(as.character(cliful) %in% names(grps[[i]]))] <- i
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
