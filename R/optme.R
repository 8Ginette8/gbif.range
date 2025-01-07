### =========================================================================
### optimization function for cluster distribution (meta)
### =========================================================================
#' Optimization function to create equal-sized strata in the 'make_blocks' function
#'
#' Not to be called directly by the user.
#' @author Philipp Brun
#' @param x Input 1
#' @param nms Input 2
#' @param grps Input 3
#' @param tot Input 4
#' @keywords internal
#' @export
#' @importFrom stats aggregate
optme=function(x,nms,grps,tot){

  # determine number of observations in each groups from initial step
  grp = sapply(grps,"sum")

  # aggregate remaining observations by suggested cluster grouping
  x = as.numeric(as.character(x))
  agg.vals = stats::aggregate(nms,by=list(x),FUN="sum")

  for (i in 1:length(grp)){

    if (i%in%agg.vals$Group.1){
      
      grp[i] = agg.vals$x[which(agg.vals$Group.1==i)]+grp[i]
    }
  }

  # Calculate difference from equal distribution
  pen = (grp-tot/length(grp))^2

  return(sum(pen))

}
