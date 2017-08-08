#' plotNormTransform
#'
#' This function plots an eset  befor and another eset after normalisation 
#'

#' @param eset_before  TODO Erklaerung
#' @param eset_after  TODO Erklaerung 
#' @param png_fn TODO Erklaerung
#' @return a png plot
#' @import data.table
#' @export

## debug




plotNormTransform = function(eset_before, eset_after, titeldet = "", png_fn = NULL) {
  
  ## plotten der transformation normalisation
  
  
  main_before = paste0("BEFORE ", titeldet, "(Total ", dim(eset_before)[1], " probes, \nvalid ", dim(eset_before)[2], " samples)")
  main_after = paste0("AFTER ", titeldet, "(Total ", dim(eset_after)[1], " probes, \nvalid ", dim(eset_after)[2], " samples)")
  # versuch = try(vsn::meanSdPlot(eset_before, main = main_before), silent = T, plot = F)
  # if (class(versuch) == "try-error") {
  p1 = vsn::meanSdPlot(eset_before, plot = F)
  p2 = vsn::meanSdPlot(eset_before, ranks = F, plot = F)
  p3 = vsn::meanSdPlot(eset_after, plot = F)
  p4 = vsn::meanSdPlot(eset_after, ranks = F, plot = F)
  
  multiplot(p1$gg + ggplot2::ggtitle(main_before), p2$gg + ggplot2::ggtitle(main_before), p3$gg + ggplot2::ggtitle(main_after), p4$gg + ggplot2::ggtitle(main_after), cols = 2)
  # } else {
  #   
  #   # par(mfrow = c(2, 2))
  #   # vsn::meanSdPlot(eset_before, main = main_before)
  #   # vsn::meanSdPlot(eset_before, main = main_before, ranks = F)
  #   # vsn::meanSdPlot(eset_after, main = main_after)
  # vsn::meanSdPlot(eset_after, main = main_after, ranks = F)
  
  
  # }
  ## speichern
  if(is.null(png_fn)==F){
    if (class(versuch) == "try-error") {
      png(png_fn, width = 700, height = 600)
      multiplot(p1$gg + ggplot2::ggtitle(main_before), p2$gg + ggplot2::ggtitle(main_before), p3$gg + ggplot2::ggtitle(main_after), p4$gg + ggplot2::ggtitle(main_after), cols = 2)
      dev.off()
    } else {
      png(png_fn, width = 600, height = 600)
      par(mfrow = c(2, 2))
      vsn::meanSdPlot(eset_before, main = main_before)
      vsn::meanSdPlot(eset_before, main = main_before, ranks = F)
      vsn::meanSdPlot(eset_after, main = main_after)
      vsn::meanSdPlot(eset_after, main = main_after, ranks = F)
      dev.off()
    }
  }
  
}