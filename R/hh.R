#' @title Show first rows and columns
#' @description Displays the first mydims rows of the first mydims columns
#' @param d object to display
#' @param mydims how many rows of how many columns to show
#' @return Matrix of the first mydims rows of the first mydims columns


#' @export



hh = function ( d, mydims=5 ) {
  # 29.1.15 data.table included
  if("data.table" %in% class(d)) { 
    d[1 : min(dim(d)[1],mydims), names(d)[ 1 : min(dim(d)[2],mydims)], with = F]
  } else   d [ 1 : min(dim(d)[1],mydims) , 1 : min(dim(d)[2],mydims) , drop =F]
}