#' @title Filter samples with atypical expressoin data
#'
#' @description Filter samples for atypical gene expression levels. This is defined as atypical large Euclidian distance of all expressed, not batch-associated expression probes to an artificial individual. This artificial individual was defined as the average of samples after removing 10% samples farthest away from the average of all samples done separately for each subgroup (implemented in the R / Bioconductor package lumi [PMID:18467348]).

#' @param ht12object A list object of class HT12prepro created with function checkBatchEffects()
#' @param paramfile Path to the file specifying parameters
#' @param filter2ind_atypischEuklid Filter for extreme combination of values from expressed transcripts summarized as Euclidian distance of expression values of all expression probes previously classified as expressed and non-overinflated regarding ANOVA on batch-effects.  Valid is ln('Euclidian distance') < median(ln('Euclidian Distance' )) + [value] *  IQR(ln('Euclidian Distance' )). If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @return A list object of class HT12prepro where the slot with  sample-related attributes of the current processing-stage named `$chipsamples` is updated.  Excluded individual are characterized by column in_study ==F and reason4exclusion. The slot with the history of the commands named `$history`` is updated.
#' @import data.table
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# ht12object =  prepro_ht12
# filter2ind_atypischEuklid =  "from_paramfile"

filterAtypicalExpressed = function(ht12object,paramfile = NULL,filter2ind_atypischEuklid = "from_paramfile") {
  
### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)
  
  myparameters = match.call()
  showVennplots = F
  
  # status checken
  historie =  ht12object$history$calls
  if(any(grepl("checkBatchEffects", historie))==F) stop("Function 'checkBatchEffects()' has to be run before!")
  
  #laden parameter
  if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)
  
  
  
  
  sample_overview_l9 <- ht12object$chipsamples
  mytable(sample_overview_l9$in_study)
  sample_overview_l9instudy <- sample_overview_l9[sample_overview_l9$in_study,]
  dim(sample_overview_l9)
  dim(sample_overview_l9instudy)
  table(table(sample_overview_l9instudy$new_ID)) 
  if(length(table(table(sample_overview_l9instudy$new_ID))) != 1)
    stop("IDs (column new_ID) must be unique....stopping...")
  
  # laden annotation probes
  subgroups <- unique(sample_overview_l9instudy$subgroup)
  subgroups
  
  genesdetail <- ht12object$genesdetail
  
  goodexpressedprobesnum <- plyr::ddply(data.frame(subgroups), "subgroups",
                                        function(x) {
                                          expressed_var <- paste0("goodexpressedprobe_", reformate_subgroup(x$subgroups[1]))
                                          res <- mytable(genesdetail[expressed_var])[, "observed", drop = F]
                                          res["TRUE", , drop = F]
                                        })
  goodexpressedprobesnum
  
  #laden expressionsets
  
  total_nobkgd_eset_ql_combat = ht12object$total_nobkgd_eset_ql_combat
  
  ## ----Euclidean-----------------------------------------------------------
  if(filter2ind_atypischEuklid== "from_paramfile") outlierCriteriumEuklid <- as.numeric(getParam2("filter2ind_atypischEuklid", myparam = param)) else outlierCriteriumEuklid = filter2ind_atypischEuklid
  outlierCriteriumEuklid
  
  
  subgroups
  
  # filtern auf exprimierte und nicht-kontrollen fuer euklidische dist der expressionswerte
  getEuklidianOutlyer <- function (subgroup, genesdetail, sample_overview_l9instudy, total_nobkgd_eset_ql_combat) {
    expressed_subgroup <- paste0("expressed_", subgroup)
    expressed_subgroup
    
    # kontrollen herauswerfen
    goodprobesNocons <- genesdetail [genesdetail$is_purecontrol == F & genesdetail[,expressed_subgroup], "nuid"]
    ind_subgroup <- sample_overview_l9instudy[sample_overview_l9instudy$subgroup == i, "new_ID"]
    all_nobkgd_expressed_eset_ql <- total_nobkgd_eset_ql_combat[goodprobesNocons, ind_subgroup]
    all_nobkgd_expressed_eset_ql
    
    # euklidische Distance berechnen bzgl. des Zentrums was allen Individuen ohne die
    # 10% extremsten Individuen (bzgl. euklidischer Distanz) entspricht
    time33 = Sys.time()
    message("Starting calculation of Euclidian distance")
    eukliddist <- lumi::detectOutlier(all_nobkgd_expressed_eset_ql,
                                      metric = "euclidean",
                                      standardize = TRUE,
                                      Th = 2,
                                      ifPlot = F)
    message("Finished calculation of Euclidian distance within ", formateTimediff(Sys.time()-time33))
    # markus mag den Th2 threshold nicht, da er unberechenbar ist, wieviel er hinauswirft,
    # dennoch ist die funktion schoen, um ueberhaupt ersteinmal die numerische eukliddist zu berechnen
    # summary(eukliddist) ##irrelevant, da ich da 3
    #als dataframe organisieren mit dem michinteressierenden Abstand
    hh(attr(eukliddist, "sampleDistance"))
    euklid <- data.frame(euklid=attr(eukliddist, "sampleDistance")[, "Center"], platzhalter = NA)
    euklid <- euklid[row.names(euklid) != "Center",]
    ht(euklid, 2)
    euklid$platzhalter <- NULL
    euklid$euklid <- log(euklid$euklid)
    
    # Define outliers
    high_cutoff <- outlierCriteriumEuklid*IQR(euklid$euklid) + median(euklid$euklid) # ohne zentrum
    high_cutoff
    bad_euklid <- row.names(euklid)[euklid$euklid > high_cutoff]
    bad_euklid
    length(bad_euklid)
    length(bad_euklid)/(dim(euklid)[1])
    res <- c()
    res$high_cutoff <- high_cutoff
    res$median <- median(euklid$euklid)
    res$bad_euklid <- bad_euklid
    res$euklid <- euklid
    res
  }
  suppl_info_euklid <- vector(mode = "list", length = length(subgroups))
  names(suppl_info_euklid) = subgroups
  suppl_info_euklid
  for(i in subgroups){
    # message("processing subgroup ", i)
    euklidinfo = getEuklidianOutlyer(subgroup = i, genesdetail,sample_overview_l9instudy, total_nobkgd_eset_ql_combat)
    suppl_info_euklid[[i]]["high_cutoff"] = euklidinfo$high_cutoff
    suppl_info_euklid[[i]]["median"] = euklidinfo$median
    suppl_info_euklid[[i]]["bad_euklid"] = list(euklidinfo$bad_euklid)
    suppl_info_euklid[[i]]["euklid"] = list(euklidinfo$euklid)
  }
  # str(suppl_info_euklid)
  num_bad_euklid <- plyr::ldply(suppl_info_euklid, function(x) length(x[["bad_euklid"]]))
  num_bad_euklid
  
  ## ----attach--------------------------------------------------------------
  sample_overview_l9$euklid <- NA
  for(i in subgroups){
    indfilter <- sample_overview_l9$subgroup == i
    euklid <- suppl_info_euklid[[i]]$euklid
    sample_overview_l9[indfilter, 'euklid'] <- euklid[match_hk(sample_overview_l9[indfilter, 'new_ID'],
                                                               row.names(euklid)), "euklid"]
    table(is.na(sample_overview_l9$euklid))
    bad_euklid <- suppl_info_euklid[[i]]$bad_euklid
    
    message("Removed ", length(bad_euklid), " Individuals in Subgroup ", i, ' due to atypical Euclidian distance of gene expression data')
    sample_overview_l9[sample_overview_l9$new_ID %in% bad_euklid, "in_study"] = F
    sample_overview_l9[sample_overview_l9$new_ID %in% bad_euklid, "reason4exclusion"] = "extreme Euclidean distance of expression values"
    xtabs_hk(~sample_overview_l9$reason4exclusion + sample_overview_l9$in_study)
    
    # plotten
    # plotten histogram der euklidischen Distanz ( stets Bezug auf Zentrum bzgl. 90% nicht extremen Individuen)
    high_cutoff <- suppl_info_euklid[[i]]$high_cutoff
    mymedian <- suppl_info_euklid[[i]]$median
    mytitle <- paste0("Filter extreme Euclidean Distance of Expression values subgroup ", i)
    # par(mfrow=c(2,1))
    
    # plotbereich = c(0, 1.2*max(c(sample_overview_l9$euklid,high_cutoff ),na.rm = T))
    plotbereich = c(min(c(sample_overview_l9$euklid,high_cutoff ), na.rm = T),
                    1.2*max(c(sample_overview_l9$euklid, high_cutoff ),na.rm = T))
    hist(euklid$euklid,
         breaks = 40,
         col = "grey",
         freq = F,
         main = mytitle,
         cex.main = 0.8,
         xlim = plotbereich);
    lines(density(euklid$euklid), col = "orange")
    abline(v = high_cutoff, col = "red")
    abline(v = mymedian, col = "blue")
    
    # #plotten boxplot
    # boxplot(euklid$euklid, main = mytitle, ylim = plotbereich,cex.main = 0.8)
    # abline(h = high_cutoff, col = "red")    
    # abline(h = mymedian, col="blue")
  }
  
  ## ----releuclid-----------------------------------------------------------
  #setDT(sample_overview_l9)
    sample_overview_l9 = data.table::data.table(data.frame(sample_overview_l9))
  
  sample_overview_l9[, mymedian := median(euklid, na.rm = T), by = subgroup]
  sample_overview_l9[, iqr:= IQR(euklid,na.rm = T ), by = subgroup]
  sample_overview_l9[, euklid_relIQR := (euklid - mymedian)/iqr, by = subgroup]
  sample_overview_l9[, range(euklid_relIQR, na.rm = T), by = subgroup]
  sample_overview_l9[in_study == T, stopifnot(max((euklid_relIQR)) <= outlierCriteriumEuklid)]
  sample_overview_l9[, mymedian := NULL]
  sample_overview_l9[, iqr := NULL]
  setDF(sample_overview_l9)
  
  ## ----detailplot----------------------------------------------------------
  detplot <- ggplot2::ggplot(sample_overview_l9[sample_overview_l9$in_study, ],
                             aes(x = euklid, fill = factor(Sentrix.Barcode))) +
    geom_density(alpha = 0.5) + facet_grid( subgroup ~ ., scale = "free_y")
  detplot 
  detplot2 <- ggplot2::ggplot(sample_overview_l9[sample_overview_l9$in_study, ],
                              aes(x = euklid, fill = subgroup)) +
    geom_density(alpha = 0.5) + facet_grid( Sentrix.Barcode ~ ., scale = "free_y")
  detplot2
  detplot3 <- ggplot2::ggplot(sample_overview_l9[sample_overview_l9$in_study, ],
                              aes(x = euklid, fill = fileset_id)) +
    geom_density(alpha = 0.5) +
    facet_grid(subgroup ~ ., scale = "free_y")
  detplot3
  detplot4 <- ggplot2::ggplot(sample_overview_l9[sample_overview_l9$in_study, ],
                              aes(x = euklid, fill = subgroup)) +
    geom_density(alpha = 0.5) +
    facet_grid( fileset_id ~ ., scale = "free_y")
  detplot4
  
  ## ----Save----------------------------------------------------------------
  # samples as used 
  ht12object$chipsamples = sample_overview_l9
  
  # status fuer docku
  dim_total_nobkgd_eset_ql_combat <- dim(total_nobkgd_eset_ql_combat)
  plots <- grep("detplot", ls(), value = T)
  fordoku <-  c("sample_overview_l9",
                "goodexpressedprobesnum",
                "num_bad_euklid",
                "suppl_info_euklid",
                "dim_total_nobkgd_eset_ql_combat",
                "outlierCriteriumEuklid",
                plots)
  
  stopifnot(sum(duplicated(fordoku))==0)
  
  
  ht12object$dokuobjects_filterAtypicalExpressed = lapply(fordoku, function(x) get(x))
  
  
  names(ht12object$dokuobjects_filterAtypicalExpressed) = fordoku
  
  
  ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
  ht12object$history
  ht12object
  
}