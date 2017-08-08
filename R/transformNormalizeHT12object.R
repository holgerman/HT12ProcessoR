#' @title Normalize and transform expression sets included in a HT12prepro object
#'
#' @description Normalizes and transforms the expression sets included in a HT12prepro object and create additional slots saving the resulting data.  See vignette for an example. 
#'

#' @param ht12object A list object of class HT12prepro created with function createExpressionSet()
#' @param paramfile Path to the file specifying parameters
#' @param normalisation_method Method used for normalisation. Either 'quantile' (quantile normalisation) or 'rsn' (robust spline normalisation). If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @param fn_plot Location where the plot showing normalization and transformation can be saved. If NULL (default), the plot will not be saved
#' @return A list object of class HT12prepro including  additional slot with an expression set with normalized and transformed data excluding control probe information named `$all_nobkgd_eset_ql`  and an additional slot with an expression set with normalized and transformed data including control probe information named `$total_nobkgd_eset_ql`    + the slot with the history of the commands named `$history`` is updated
#' @import data.table
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# normalisation_method = "from_paramfile"
# ht12object =  prepro_ht12
# fn_plot = "from_paramfile"

transformNormalizeHT12object = function(ht12object,paramfile = NULL, normalisation_method = "from_paramfile", fn_plot = NULL ) {
### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)
  
  myparameters = match.call()
  showVennplots = F
  
  
  # status checken
  historie =  ht12object$history$calls
  if(any(grepl("createExpressionSet", historie))==F) stop("Function 'createExpressionSet()' has to be run before!")
  
  #laden parameter
  if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)
  
  
  
sample_overview_l6 <- ht12object$chipsamples
mytable(sample_overview_l6$in_study)
sample_overview_l6instudy <- sample_overview_l6[sample_overview_l6$in_study, ]
dim(sample_overview_l6)
dim(sample_overview_l6instudy)
table(table(sample_overview_l6instudy$new_ID)) 
if(length(table(table(sample_overview_l6instudy$new_ID))) != 1)
  stop("IDs (column new_ID) must be unique....stopping...")


# laden expressionsets


total_nobkgd_eset =  ht12object$total_nobkgd_eset


## ----process6------------------------------------------------------------
goodind  <- sample_overview_l6instudy$new_ID
# str(goodind)
total_nobkgd_eset
message("Transforming and Normalizing expression sets of ", length(goodind), " samples still good according to column 'in_study'....")
total_nobkgd_eset_goodind <- total_nobkgd_eset[, goodind]
# total_nobkgd_eset_goodind <- total_nobkgd_eset_goodind[,goodind[1:50]]  # debug
total_nobkgd_eset_goodind

## ----normandtrans--------------------------------------------------------
if(normalisation_method == "from_paramfile")  methodtransform <- getParam2("normalisation_method", myparam = param) else  methodtransform = normalisation_method  
methodtransform

total_nobkgd_eset_ql <- transformNormalize(total_nobkgd_eset_goodind,
                                           methodtransform = methodtransform,
                                           dolog2 = T)
                                          
total_nobkgd_eset_ql



plotNormTransform(total_nobkgd_eset_goodind, total_nobkgd_eset_ql, png_fn = fn_plot)

## ----speichern---------------------------------------------------------


ht12object$total_nobkgd_eset_ql = total_nobkgd_eset_ql
ht12object$dokuobjects_transformNormalizeHT12object$normmethod = methodtransform
ht12object$all_nobkgd_eset_ql = total_nobkgd_eset_ql[Biobase::featureNames(ht12object$all_nobkgd_eset),]
ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))

ht12object
}