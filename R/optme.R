### =========================================================================
### optimization function for cluster distribution (meta)
### =========================================================================
#' Objective Function for \code{make_blocks()}
#'
#' Internal helper used by \code{make_blocks()} to distribute clusters across
#' folds as evenly as possible.
#'
#' @param x Candidate allocation of the remaining clusters to folds.
#' @param nms Named vector of cluster sizes.
#' @param grps Current fold assignments built before optimization.
#' @param tot Total number of observations across all clusters.
#' @keywords internal
#' @importFrom stats aggregate
optme <- function(x, nms, grps, tot){

  # determine number of observations in each groups from initial step
  grp <- vapply(grps, sum, FUN.VALUE = numeric(1))

  # aggregate remaining observations by suggested cluster grouping
  x <- as.numeric(as.character(x))
  agg.vals <- stats::aggregate(nms, by = list(x), FUN = "sum")

  for (i in seq_along(grp)){

    if (i %in% agg.vals$Group.1){
      
      grp[i] <- agg.vals$x[which(agg.vals$Group.1 == i)] + grp[i]
    }
  }

  # Calculate difference from equal distribution
  pen <- (grp - tot / length(grp))^2

  return(sum(pen))

}
