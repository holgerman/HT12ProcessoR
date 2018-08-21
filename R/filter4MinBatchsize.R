#' @title Filter samples of a HT12prepro object for min batch size n>1
#'
#' @description This functions filters for min Batch size n>1  as this is a limitation for following removal of batch effects in the following step applying function removeBatchEffects()
#'
#'

#' @param ht12object A list object of class HT12prepro created with function filterTechnicallyFailed()
#' @param paramfile Path to the file specifying parameters
#' @param second_combat_withvar A second round of removeBatchEffects() can be intended. In this case, batches from the first round checking Sentrix.Barcode batches as well as batches from this second round specified as a column in ht12object$chipsamples named according to this variable must be checked for minimum batch size n>1.If "", no 2nd combat will be done. If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#'
#' @return A list object of class HT12prepro where the  slot with  sample-related attributes of the current processing-stage named `$chipsamples` is update. Excluded individual are characterized by column in_study ==F and reason4exclusion
#' @import data.table
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# ht12object =  prepro_ht12
# second_combat_withvar= "from_paramfile"

filter4MinBatchsize = function(ht12object,paramfile = NULL, second_combat_withvar= "from_paramfile") {

### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)
myparameters = match.call()
showVennplots = F

# status checken
historie =  ht12object$history$calls
if(any(grepl("filterTechnicallyFailed", historie))==F) stop("Function 'filterTechnicallyFailed()' has to be run before!")

#laden parameter
if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)

ht12object$chipsamples = as.character(ht12object$chipsamples) ## added 22.8.18 , da int64 nicht  wie eine numerische Variable mit einem character verglichen werden kann
sample_overview_l7pre <-ht12object$chipsamples
mytable(sample_overview_l7pre$in_study)
sample_overview_l7preinstudy <- sample_overview_l7pre[ sample_overview_l7pre$in_study, ]
dim(sample_overview_l7pre)
dim(sample_overview_l7preinstudy)
table(table(sample_overview_l7preinstudy$new_ID))
if(length(table(table(sample_overview_l7preinstudy$new_ID))) != 1)
  stop("IDs (column new_ID) must be unique....stopping...")


#laden expressionsets
total_nobkgd_eset_ql = ht12object$total_nobkgd_eset_ql

## ----goodind-------------------------------------------------------------
# auszahlen aktualisiere barcode
table(table(sample_overview_l7pre$Sentrix.Barcode))
# table(sample_overview_l7pre[grep("mis", sample_overview_l7pre$Sentrix.Barcode), 5])
Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber <- sample_overview_l7pre[match_hk(Biobase::pData(total_nobkgd_eset_ql)$sampleID,
                                                                                                      sample_overview_l7pre$new_ID), "Sentrix.Barcode" ]
check <- showNA(Biobase::pData(total_nobkgd_eset_ql))$NAs
if(sum(check) != 0)
  stop("not all smples appear to ahave a Sentrix ID!")
ht(Biobase::pData(total_nobkgd_eset_ql), 1)

# gute individuen
goodind  <- sample_overview_l7preinstudy$new_ID
# str(goodind)
total_nobkgd_eset_ql
total_nobkgd_eset_ql <- total_nobkgd_eset_ql[,goodind]
total_nobkgd_eset_ql

## ----criteria2-----------------------------------------------------------
message("Filtering for batch size > 1 for Sentrix IDs...")



tabled <- table(Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber)
table(tabled)
singlbarcoders <- names(tabled[tabled == 1])
singlbarcoders
singlbarcoders_ind <- sample_overview_l7preinstudy[sample_overview_l7preinstudy$Sentrix.Barcode %in% singlbarcoders, "new_ID"]
singlbarcoders_ind
# str(goodind)
goodind <- setdiff(goodind, singlbarcoders_ind)
# str(goodind)
total_nobkgd_eset_ql <- total_nobkgd_eset_ql[,goodind]
total_nobkgd_eset_ql




## ----combat8b------------------------------------------------------------
head(Biobase::pData(total_nobkgd_eset_ql))
Biobase::pData(total_nobkgd_eset_ql)$subgroup <- sample_overview_l7pre[ match_hk(Biobase::pData(total_nobkgd_eset_ql)$sampleID,
                                                                                 sample_overview_l7pre$new_ID), "subgroup"]
showNA(Biobase::pData(total_nobkgd_eset_ql))
singularcheck <- data.table(Biobase::pData(total_nobkgd_eset_ql))
singularcheck2 <- singularcheck[,.N, by = list(hybridisierungchipserialnumber, subgroup) ]
singularcheck3 <- singularcheck2[allDuplicatedEntries(hybridisierungchipserialnumber)]
singularcheck3



if(second_combat_withvar== "from_paramfile") second_combat_withvar_used <- getParam2("second_combat_withvar", myparam = param) else second_combat_withvar_used = second_combat_withvar
second_combat_withvar_used

if(second_combat_withvar_used != "") {
  message("Filtering for batch size > 1 as adjusting genexpression data in a second combat round for `",
          second_combat_withvar_used, "` is planned...")

  # Add data
  Biobase::pData(total_nobkgd_eset_ql_combat)[,second_combat_withvar_used] <-  sample_overview_l7pre[match_hk(Biobase::pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                                                              sample_overview_l7pre$new_ID),
                                                                                                     second_combat_withvar_used]
  check <- showNA(Biobase::pData(total_nobkgd_eset_ql_combat))
  if(sum(check) != 0)
    stop("Not all individuals have a Sentrix ID!")
  ht(Biobase::pData(total_nobkgd_eset_ql_combat), 1)

  # Exclude single barcoders
  tabled <- table(Biobase::pData(total_nobkgd_eset_ql_combat)[, second_combat_withvar_used])
  table(tabled)
  newsinglbarcoders <- names(tabled[tabled == 1])
  newsinglbarcoders
  newsinglbarcoders_ind <- sample_overview_l7preinstudy[sample_overview_l7preinstudy$Sentrix.Barcode %in% newsinglbarcoders, "new_ID"]
  newsinglbarcoders_ind
  # str(goodind)
  goodind <- setdiff(goodind, newsinglbarcoders_ind)
  # str(goodind)
}

# individuenattrib
sample_overview_l7pre[sample_overview_l7pre$new_ID %in% singlbarcoders_ind,"in_study"] = F
sample_overview_l7pre[sample_overview_l7pre$new_ID %in% singlbarcoders_ind,"reason4exclusion"] <- "Filtering as Batch/effect correction via ComBat not possible - only 1 Ind. within a batch  found"
mytable(sample_overview_l7pre$reason4exclusion)
mytable(sample_overview_l7pre$in_study)
ht12object$chipsamples =  sample_overview_l7pre

all_removed_singlbarcoders = unique(na.omit(c(singlbarcoders_ind, singlbarcoders)))
message("Removed ", length(all_removed_singlbarcoders), " samples where only a single sample was found in at least one batch. (Removed IDs: ", paste(all_removed_singlbarcoders, collapse = ", "), ")")

fordoku =c(
  "singlbarcoders_ind",
  "singlbarcoders",
  "sample_overview_l7pre"
)




stopifnot(sum(duplicated(fordoku))==0)


ht12object$dokuobjects_filter4MinBatchsize = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_filter4MinBatchsize) = fordoku


ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object$history
ht12object


}
