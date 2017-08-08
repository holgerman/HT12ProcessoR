#' @title Show first and last rows 
#' @description Displays the first and last rows of a table
#' @param d object to display
#' @param myrows how many of the first and last rows to show
#' @return Matrix of the first and last myrow rows 


#' @export



ht = function ( d, myrows=10 ) 
{ ## updated 11.3. to show all if dim 1 smaller than myrows*2
  rows2show = min(dim(d)[1],myrows)
  if(dim(d)[1] <= 2*rows2show) return(d) 
  rbind ( head ( d ,  rows2show ), tail ( d ,  rows2show ))
}