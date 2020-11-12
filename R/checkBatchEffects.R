#' @title Identify remaining batch effects in expression data
#'
#' @description Identification of transcripts where association with batches is still stronger than expected by chance after batch-adjustment. Stronger than chance (i.e. overinflation) is defined as association stronger p-value after Bonferroni correction. All defined Batches will be tested, i.e. `Sentrix.Barcode`,  `processingbatch`, and `strangebatch` in the table from slot `$chipsamples` if more than one category of the batch is found. Probes are not filtered but annotated when found to be overinflated. This function makes use of package MatrixEQTL in order to allow fast ANOVA on matrices. Therefore, some messages from these functions report doing an "eQTL analysis" allthough in fact an ANOVA is done on batch effects

#' @param ht12object A list object of class HT12prepro created with function removeBatchEffects()
#' @param paramfile Path to the file specifying parameters
#' @param round4ANOVAcheck parameter for control how identical Results of MATRIXEQTL and standard ANova have to be. Provided number rounds the F statistic value. Applied in function runMAtrixEQTLAnova called by calcAnovaSentrixRunSpecialbatchViaMatrixEQTL2() . Note that slightly different results may occur  when e.g. using two completely nested covariates
#' @return A list object of class HT12prepro where  the slot with  probe-related attributes of the current processing-stage named `$genesdetail` is updated as well as the  slot with the history of the commands named `$history`.  QQ plots of the association are shown.
#' @import data.table
#' @export

## debug
# paramfile = myparamfile
# ht12object =  copy(prepro_sorbv2)
# showPlots = T



checkBatchEffects = function(ht12object,paramfile = NULL, showPlots=T,round4ANOVAcheck = 5) {

### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

  myparameters = match.call()
  showVennplots = F

  # status checken
  historie =  ht12object$history$calls
  if(any(grepl("removeBatchEffects", historie))==F) stop("Function 'removeBatchEffects()' has to be run before!")

  #laden parameter
  if(is.null(paramfile)==F) param <- data.frame(data.table::fread(paramfile))



sample_overview_l8 <-  ht12object$chipsamples
sample_overview_l8instudy <- sample_overview_l8[sample_overview_l8$in_study, ]
dim(sample_overview_l8)
dim(sample_overview_l8instudy)
table(table(sample_overview_l8instudy$new_ID))
if(length(table(table(sample_overview_l8instudy$new_ID))) != 1)
  stop("IDs (column new_ID) must be unique....stopping...")



# laden annotation probes
genesdetail <- ht12object$genesdetail

ht(genesdetail, 2)
mytable(stringr::str_sub(genesdetail$ilmn, 1, 4))
subgroups <- unique(sample_overview_l8instudy$subgroup)



#laden expressionsets nach combat
total_nobkgd_eset_ql_combat = ht12object$total_nobkgd_eset_ql_combat
total_nobkgd_eset_ql_combat
#laden expressionsets VOR combat
total_nobkgd_eset_ql = ht12object$total_nobkgd_eset_ql
total_nobkgd_eset_ql

# laden specialbatch
anova_with_special_batch <- unique(sample_overview_l8instudy$strangebatch)
anova_with_special_batch
anova_with_special_batch <- anova_with_special_batch %in% c("", 0) == F
anova_with_special_batch
stopifnot(is.logical(anova_with_special_batch))

## ----goodind-------------------------------------------------------------
# gute individuen post combat schon gefiltert
total_nobkgd_eset_ql_combat
total_nobkgd_eset_ql

# gute individuen
goodind <- sample_overview_l8instudy$new_ID
# str(goodind)
check1 <- venn3(goodind,
                pData(total_nobkgd_eset_ql)$sampleID,
                pData(total_nobkgd_eset_ql_combat)$sampleID, plotte = showVennplots)
# str(check1)
# removed check1$q6 from this, which is the overlap of samples pre and post combat - which does not have to be 0
stopifnot(length(c(check1$q2, check1$q4, check1$q5,  check1$q7, check1$q3)) == 0)

# Its no longer valid to compare sample size before/after combat due to removal of IDs... check on other ways to successfully check this
# venn2(ht12object$chipsamples[ht12object$chipsamples$in_study==F, "new_ID"],
# check1$q6)

# i might have to remove check1$q4 in case samples with NAs in the covariates get removed during ComBat

## ----AnnotationCombatCheck-----------------------------------------------
ht(pData(total_nobkgd_eset_ql_combat), 2)
total_nobkgd_eset_ql_combat
pData(total_nobkgd_eset_ql)$Sentrix.Barcode <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql)$sampleID,
                                                                           sample_overview_l8$new_ID),
                                                                  "Sentrix.Barcode"]
pData(total_nobkgd_eset_ql)$bigbatch <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql)$sampleID,
                                                                    sample_overview_l8$new_ID),
                                                           "fileset_id"]
pData(total_nobkgd_eset_ql)$strangebatch <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql)$sampleID,
                                                                        sample_overview_l8$new_ID),
                                                               "strangebatch"]
pData(total_nobkgd_eset_ql_combat)$Sentrix.Barcode <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                                  sample_overview_l8$new_ID),
                                                                         "Sentrix.Barcode"]
pData(total_nobkgd_eset_ql_combat)$bigbatch <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                               sample_overview_l8$new_ID),
                                                                      "fileset_id"]
pData(total_nobkgd_eset_ql_combat)$strangebatch <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                               sample_overview_l8$new_ID),
                                                                      "strangebatch"]

pData(total_nobkgd_eset_ql)$subgroup <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql)$sampleID,
                                                                        sample_overview_l8$new_ID),
                                                               "subgroup"]

pData(total_nobkgd_eset_ql_combat)$subgroup <- sample_overview_l8[match_hk(pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                               sample_overview_l8$new_ID),
                                                                      "subgroup"]


check <- sum(showNA(pData(total_nobkgd_eset_ql))$NAs)
if(check != 0)
  stop("check phenotypes for missing items")
check <- sum(showNA(pData(total_nobkgd_eset_ql_combat))$NAs)
if(check != 0)
  stop("check phenotypes for missing items")

## ----overinflated--------------------------------------------------------
subgroups
anova_with_special_batch
n_subgroups <- length(unique(pData(total_nobkgd_eset_ql_combat)$subgroup))
n_subgroups

showNA(genesdetail)
ht(pData(total_nobkgd_eset_ql),1)
anovres <- calcAnovaSentrixRunSpecialbatchViaMatrixEQTL2(total_nobkgd_eset_ql,
                                                         total_nobkgd_eset_ql_combat,
                                                         genesdetail,
                                                         subgroups,
                                                         anova_with_special_batch,
                                                         strangebatch,
                                                         adjust4subgroup = (n_subgroups>1),
                                                         round4ANOVAcheck = round4ANOVAcheck
                                                         )
genesdetail <- anovres$genesdetail
ht(genesdetail, 2)
bonf_pwert <- anovres$bonf_pwert
bonf_pwert
showNA(genesdetail)
mytable(genesdetail$all_pval_after_sentrix_bf)
mytable(genesdetail$all_pval_before_sentrix_bf)
mytable(genesdetail$all_pval_after_bigbatch_bf)
mytable(genesdetail$all_pval_before_bigbatch_bf)
xtabs_hk(~genesdetail$all_pval_after_sentrix_bf + genesdetail$all_pval_after_bigbatch_bf)
xtabs_hk(~genesdetail$all_pval_after_sentrix_bf + genesdetail$all_pval_after_strangebatch_bf)
xtabs_hk(~genesdetail$all_pval_after_bigbatch_bf + genesdetail$all_pval_after_strangebatch_bf)
xtabs_hk(~genesdetail$all_pval_before_sentrix_bf + genesdetail$all_pval_before_bigbatch_bf)
xtabs_hk(~genesdetail$all_pval_before_sentrix_bf + genesdetail$all_pval_before_strangebatch_bf)
xtabs_hk(~genesdetail$all_pval_before_bigbatch_bf + genesdetail$all_pval_before_strangebatch_bf)

if(showPlots==T) {
  par(mfrow =c(2,2))
  n_sentrix = length(unique(sample_overview_l8$Sentrix.Barcode[ sample_overview_l8$in_study]))
  n_fileset_id = length(unique(sample_overview_l8$fileset_id[ sample_overview_l8$in_study]))

  # str(genesdetail)
  ggd.qqplot(genesdetail$all_pval_before_sentrix, main= "ANOVA Sentrix before Combat")
  mtext(paste(n_sentrix, " batch(es) present", collapse = " "))
  ggd.qqplot(genesdetail$all_pval_after_sentrix, main= "ANOVA Sentrix after Combat")
  mtext(paste(n_sentrix, " batch(es) present", collapse = " "))
  ggd.qqplot(genesdetail$all_pval_before_bigbatch, main= "ANOVA fileset_id before Combat")
  mtext(paste(n_fileset_id, " batch(es) present", collapse = " "))
  ggd.qqplot(genesdetail$all_pval_after_bigbatch, main= "ANOVA fileset_id after Combat")
  mtext(paste(n_fileset_id, " batch(es) present", collapse = " "))


}

## ----save----------------------------------------------------------------

ht12object$genesdetail = genesdetail

# for doku
dim_total_nobkgd_eset_ql_combat <- dim(total_nobkgd_eset_ql_combat)
dim_total_nobkgd_eset_ql_combat
dim_total_nobkgd_eset_ql <- dim(total_nobkgd_eset_ql)
dim_total_nobkgd_eset_ql

fordoku =c("sample_overview_l8",
             "genesdetail",
             "bonf_pwert",
             "dim_total_nobkgd_eset_ql_combat",
           'anova_with_special_batch')




stopifnot(sum(duplicated(fordoku))==0)


ht12object$dokuobjects_checkBatchEffect = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_checkBatchEffect) = fordoku


ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object$history
ht12object

}
