#' @title Initial PCA of  attributes from sample files originatting from Illumina GenomeStudio
#' @description Numerical Attributes from the slot `$chipsamples` from an HT12prepro object are dimension-reduced by PCR
#' @param ht12object A list object of class HT12prepro
#' @return a list object including slots `pca_object` which is a copy of the $chipsamples  with additional columns representing the principal components. Additionally, `pc.data` is the PCA-object created by function princomp()

#' @import data.table
#' @export

calcInitialPCA = function(ht12object) {
  ## debug
  # ht12object =  prepro_ht12

  ### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

  stopifnot( "HT12prepro" %in% class(ht12object))

  sample_overview_l3 =   ht12object$chipsamples


  numcols <- names(sample_overview_l3)[sapply(sample_overview_l3, is.numeric)]
  numcols
  numcols <- setdiff(numcols, c("Index", "Sentrix.Barcode"))
  data.pc <- sample_overview_l3[ is.na(sample_overview_l3$old_ID) == F & sample_overview_l3$in_study == T, numcols]
  names(data.pc)
  data.pc <- sapply(data.pc, scale)
  class(data.pc)
  rownames(data.pc) <- sample_overview_l3[is.na(sample_overview_l3$old_ID) == F & sample_overview_l3$in_study == T, 'new_ID']
  hh(data.pc)
  # throw out columns without data
  nas <- showNA(data.frame(data.pc))
  badcolna <- names(nas[nas == dim(data.pc)[1]])
  data.pc <- data.pc[ ,colnames(data.pc) %nin% badcolna]
  data.pc <- data.pc[ ,showNA(data.frame(data.pc))$NAs == 0]
  # pc.data <- princomp(data.pc, cor = T, scores = T)
  pc.data <- prcomp(data.pc)
  plot(pc.data, main="Eigenvalues")
  biplot(pc.data)

  pc123 = data.frame(new_ID = rownames(pc.data$x), pc.data$x)
  ht(pc123)
  pc123 = merge(pc123, sample_overview_l3, by = "new_ID")

  ht(pc123,1)

  ## ----save----------------------------------------------------------------
  res = c()
  res$pca_object = pc.data
  res$scores = pc123
  return(res)

}
