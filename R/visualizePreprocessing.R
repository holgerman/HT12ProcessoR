#' @title Visualize expression data before and after preprocessing
#'
#' @description Visualisation of control features and expression data before and after preprocessing in several plots. Thereby, raw data before preprocessing is log2 transformed for a more meaningful comparison

#' @param ht12object A list object of class HT12prepro created with function filterAtypicalExpressed()
#' @param paramfile Path to the file specifying parameters
#' @param show_qcplots_individSamplelevel Show QC plots comparing expression values before / after preprocessing either binned on SentrixID ('FALSE') or binned on Sample ID ('TRUE')? If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @param showPlots Show plots. If FALSE, plots are only stored in slot $dokuobjects_visualizePreprocessing
#' @return A list object of class HT12prepro including additional slots  $dokuobjects_visualizePreprocessing  
#' @import data.table
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# ht12object =  prepro_ht12
# show_qcplots_individSamplelevel =  "from_paramfile"
# showPlots = T

visualizePreprocessing = function(ht12object,paramfile = NULL,show_qcplots_individSamplelevel = "from_paramfile", showPlots=T) {
### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)
  
  myparameters = match.call()
  showVennplots = F
  
  # status checken
  historie =  ht12object$history$calls
  if(any(grepl("filterAtypicalExpressed", historie))==F) stop("Function 'filterAtypicalExpressed()' has to be run before!")
  
  #laden parameter
  if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)
  
  
  
  
  sample_overview_l9 <- ht12object$chipsamples
  
  
  
  ## ----functions-----------------------------------------------------------
  annotateEset <- function (e1, annotcon, attrib = sample_overview_l10instudy, showPlots = T) {
    # newnames = attrib[match(names(e1), attrib$new_ID), "old_ID"]
    # grep(" ", attrib$old_ID, value = T)
    # newnames =  stringr::str_replace_all(newnames, " ", "")
    # names(e1) = newnames
    
    e1$ilmn <- suppressWarnings(lumi::nuID2IlluminaID(row.names(e1)))[, 'Probe_Id']
    e1[grep("ILMN", row.names(e1)), "ilmn"] <- row.names(e1[grep("ILMN", row.names(e1)),])
    # print(table(is.na(e1$ilmn)))
    e1 <- moveColFront(e1, "ilmn")
    e1 <- merge(e1, annotcon[, c("ilmn", "Reporter_Group_Name", "Reporter_Group_id")], by = "ilmn", all.x = T)
    e1 <- moveColFront(e1, "Reporter_Group_Name")
    e1 <- moveColFront(e1, "Reporter_Group_id")
    e1
  }
  
  makePrintDF <- function (e2, var_of_interest, attrib = sample_overview_l10instudy) {
    # debug 
    # var_of_interest =  "housekeeping"
    # attrib = sample_overview_l10instudy
    d <- e2[(e2$Reporter_Group_id == var_of_interest & is.na(e2$Reporter_Group_id) == F)|
              (e2$Reporter_Group_Name == var_of_interest & is.na(e2$Reporter_Group_Name) == F), ]
    # print(hh(d))
    d <- t(d[, 4:(dim(d)[2])])
    myid <- row.names(d)
    d <- data.frame(d)
    d$id <- myid
    d$Sentrix.Barcode <- attrib[match_hk(d$id, attrib$new_ID), "Sentrix.Barcode"]
    # print(table(d$Sentrix.Barcode, useNA = "always"))
    d$Sentrix.Barcode <- attrib[match_hk(d$id, attrib$new_ID), "Sentrix.Barcode"]
    # print(table(d$Sentrix.Barcode, useNA = "always"))
    d$subgroup <- attrib[match_hk(d$id, attrib$new_ID), "subgroup"]
    # print(mytable(d$subgroup))
    dm <- melt(d, id.vars = c("id", "Sentrix.Barcode", "subgroup"))
    dm$value <- as.numeric(dm$value)
    names(dm)[names(dm) == "value"] <- var_of_interest
    dm$fileset_id <- attrib[match_hk(d$id, attrib$new_ID), "fileset_id"]
    dm$strange_batch <- attrib[match_hk(d$id, attrib$new_ID), "strangebatch"]
    names(dm) <- stringr::str_replace_all(names(dm), ":", "_")
    # print(ht(dm, 2))
    dm
  }
  
  plotGenes <- function (var_of_interest = "housekeeping", e1m) {  
    pd <- e1m[, c(var_of_interest, 'Sentrix.Barcode', 'fileset_id', 'strange_batch')]
    pd <- pd[order(pd$fileset_id), ]
    pd$Sentrix.Barcode <- factor(pd$Sentrix.Barcode, levels  = unique(pd$Sentrix.Barcode))
    pd$strange_batch <- as.numeric(factor(pd$strange_batch))
    # print(ht(pd))
    p1 <- ggplot2::ggplot(pd, ggplot2::aes_string(x = 'Sentrix.Barcode',
                                                  y = var_of_interest,
                                                  group = 'Sentrix.Barcode',
                                                  fill = 'fileset_id',
                                                  col = 'fileset_id',
                                                  alpha = 'strange_batch')) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(col = FALSE, alpha = F, fill = guide_legend(ncol = 4)) +
      ggplot2::theme(axis.text.x = element_blank(), legend.text = element_text(size = 8)) +
      ggplot2::theme(legend.position="top")
    p1
  }
  plotGenes2 <- function (var_of_interest = "housekeeping", e1m) {  
    pd <- e1m[, c(var_of_interest, 'Sentrix.Barcode', 'fileset_id', 'strange_batch', "subgroup")]
    pd <- pd[order(pd$fileset_id), ]
    pd$Sentrix.Barcode <- factor(pd$Sentrix.Barcode, levels = unique(pd$Sentrix.Barcode))
    pd$strange_batch <- as.numeric(factor(pd$strange_batch))
    # print(ht(pd))
    p1 <- ggplot2::ggplot(pd, ggplot2::aes_string(x = 'Sentrix.Barcode',
                                                  y = var_of_interest,
                                                  group = 'Sentrix.Barcode',
                                                  fill = 'fileset_id',
                                                  col = 'fileset_id',
                                                  alpha = 'strange_batch')) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(col = FALSE, alpha = F, fill = guide_legend(ncol = 4)) +
      ggplot2::theme(axis.text.x = element_blank(), legend.text = element_text(size = 8)) +
      ggplot2::facet_grid(. ~subgroup, scales = "free_x") +
      ggplot2::theme(legend.position="top")
    p1
  }
  plotGenes2log2 <- function (var_of_interest =  "housekeeping", e1m) {  
    pd <- e1m[, c(var_of_interest, 'Sentrix.Barcode', 'fileset_id', 'strange_batch', "subgroup")]
    pd <- pd[order(pd$fileset_id), ]
    pd$Sentrix.Barcode <- factor(pd$Sentrix.Barcode, levels = unique(pd$Sentrix.Barcode))
    pd$strange_batch <- as.numeric(factor(pd$strange_batch))
    pd[, var_of_interest] <- log2( pd[, var_of_interest])
    # print(ht(pd))
    p1 <- ggplot2::ggplot(pd, ggplot2::aes_string(x = 'Sentrix.Barcode',
                                                  y = var_of_interest,
                                                  group = 'Sentrix.Barcode',
                                                  fill = 'fileset_id',
                                                  col = 'fileset_id',
                                                  alpha = 'strange_batch')) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(col = FALSE, alpha = F, fill = guide_legend(ncol = 4)) +
      ggplot2::theme(axis.text.x = element_blank(), legend.text = element_text(size = 8)) +
      ggplot2::facet_grid(. ~subgroup, scales = "free_x") +
      ylab(paste("log2", var_of_interest))+
      ggplot2::theme(legend.position="top")
    p1
  }
  plotGenes3 <- function (var_of_interest = "housekeeping", e1m) {  
    pd <- e1m[, c(var_of_interest, 'Sentrix.Barcode', 'fileset_id', 'strange_batch', "subgroup", "id")]
    pd$Sentrix.Barcode <- factor(pd$Sentrix.Barcode, levels = unique(pd$Sentrix.Barcode))
    pd$strange_batch <- as.numeric(factor(pd$strange_batch))
    pd <- pd[order(pd$fileset_id, pd$Sentrix.Barcode), ]
    pd$id <- factor(pd$id, levels = unique(pd$id))    
    # print(ht(pd))
    p1 <- ggplot2::ggplot(pd, ggplot2::aes_string(x = 'id',
                                                  y = var_of_interest,
                                                  group = 'id',
                                                  fill = 'Sentrix.Barcode',
                                                  col = 'Sentrix.Barcode',
                                                  alpha = 'strange_batch')) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(col = FALSE, alpha = F, fill = guide_legend(ncol = 4)) +
      ggplot2::theme(axis.text.x = element_blank(), legend.text = element_text(size = 8)) +
      ggplot2::facet_grid(. ~subgroup, scales = "free_x")+
      ggplot2::theme(legend.position="top")
    p1
  }
  plotGenes3log2 <- function (var_of_interest = "housekeeping", e1m) {  
    pd <- e1m[, c(var_of_interest, 'Sentrix.Barcode', 'fileset_id', 'strange_batch', "subgroup", "id")]
    pd$Sentrix.Barcode <- factor(pd$Sentrix.Barcode, levels = unique(pd$Sentrix.Barcode))
    pd$strange_batch <- as.numeric(factor(pd$strange_batch))
    pd[, var_of_interest] <- log2(pd[,var_of_interest])
    pd <- pd[order(pd$fileset_id, pd$Sentrix.Barcode), ]
    pd$id <- factor(pd$id, levels = unique(pd$id))
    # print(ht(pd))
    p1 <- ggplot2::ggplot(pd, ggplot2::aes_string(x = 'id',
                                                  y = var_of_interest,
                                                  group = 'id',
                                                  fill = 'Sentrix.Barcode',
                                                  col = 'Sentrix.Barcode' ,
                                                  alpha = 'strange_batch')) +
      ggplot2::geom_boxplot() +
      ggplot2::guides(col = FALSE, alpha = F,  fill = guide_legend(ncol = 4)) +
      ggplot2::theme(axis.text.x = element_blank(), legend.text = element_text(size=8)) +
      ggplot2::facet_grid(. ~subgroup, scales = "free_x") + ylab(paste("log2", var_of_interest))+
      ggplot2::theme(legend.position="top")
    p1
  }
  
  sample_overview_l10 <- ht12object$chipsamples
  mytable(sample_overview_l10$in_study)
  mytable(sample_overview_l10$reasons4exclusion)
  
  ## ----load92--------------------------------------------------------------
  sample_overview_l10instudy <- sample_overview_l10[ sample_overview_l10$in_study, ]
  dim(sample_overview_l10)
  dim(sample_overview_l10instudy)
  table(table(sample_overview_l10instudy$new_ID)) 
  if(length(table(table(sample_overview_l10instudy$new_ID))) != 1)
    stop("modify code, script expects ids in column new_ID!")
  
  # laden annotation probes
  genesdetail <- ht12object$genesdetail
  ht(genesdetail, 2)
  mytable(stringr::str_sub(genesdetail$ilmn, 1, 4))
  
  #laden expressionsets nach combat
  total_nobkgd_eset_ql_combat = ht12object$total_nobkgd_eset_ql_combat
  total_nobkgd_eset_ql_combat
  
  total_nobkgd_eset= ht12object$total_nobkgd_eset
  total_nobkgd_eset
  
  # details kontrollinfo
  head(annotcon)
  
  if("Probe_Id" %in% names(annotcon)) annotcon <- reshape::rename(annotcon, c(Probe_Id = "ilmn"))
  head(annotcon)
  
  unique(annotcon$Reporter_Group_Name)
  
  annotcon$nuid <- suppressWarnings(lumi::IlluminaID2nuID(IlluminaID = annotcon$ilmn))[, "nuID"]
  annotcon[is.na(annotcon$nuid), "nuid"] <- annotcon[is.na(annotcon$nuid), "ilmn"]
  ht(annotcon, 2)
  
  
  ## ----add-----------------------------------------------------------------
  # details erccc probes
  head(ercc_det)
  matcher <- annotcon[annotcon$Reporter_Group_Name %in% ercc_det$ERCC.ID, ]
  ercc_det$ilmn <- matcher[match_hk(ercc_det$ERCC.ID, matcher$Reporter_Group_Name), "ilmn"]
  ercc_det$ercc_gruppe1 <- cut(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.),
                               breaks = c(-3,0,1,2,3,4,5), include.lowest = T )
  mytable(ercc_det$ercc_gruppe1)
  ercc_det$ercc_gruppe1short <- paste0("ercc_", stringr::str_sub(sapply(stringr::str_split(ercc_det$ercc_gruppe1, ","), "[",2),1,1))
  mytable(ercc_det$ercc_gruppe1short)
  ht(ercc_det, 2)
  ht(annotcon, 2)
  toadd <- annotcon[0, ][1:dim(ercc_det)[1], ]
  toadd$ilmn <- ercc_det$ilmn
  matcher <- unique(annotcon[, c('nuid', "ilmn")])
  toadd$nuid <- matcher[match_hk(toadd$ilmn, matcher$ilmn), "nuid"]
  toadd$Reporter_Group_Name <- "ercc_det"
  toadd$Reporter_Group_id <- ercc_det$ercc_gruppe1short
  ht(toadd, 3)
  annotcon <- rbind(annotcon, toadd)
  
  ## ----good----------------------------------------------------------------
  xtabs_hk(~ sample_overview_l10$reason4exclusion + sample_overview_l10$in_study)
  goodind <- sample_overview_l10instudy$new_ID
  # str(goodind)
  total_nobkgd_eset
  total_nobkgd_eset <- total_nobkgd_eset[,goodind]
  total_nobkgd_eset
  total_nobkgd_eset_ql_combat  # das rausnehmen nach etablierung der 2. Runde normal combat nach euklid
  total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql_combat[,goodind]
  total_nobkgd_eset_ql_combat
  qlist33 <- venn2(rownames(exprs(total_nobkgd_eset_ql_combat)),rownames(exprs(total_nobkgd_eset)), plotte = showVennplots)
  # stopifnot(identical(dim(total_nobkgd_eset_ql_combat), dim(total_nobkgd_eset)))
  
  ## ----looking, fig.width=20, fig.height=12--------------------------------
  e1 <- data.frame(exprs(total_nobkgd_eset))
  e1 <- annotateEset(e1, annotcon)
  hh(e1)
  dim(e1)
  
  # das gnze fuer e2m
  e2 <- data.frame(exprs(total_nobkgd_eset_ql_combat))
  e2 <- annotateEset(e2, annotcon)
  hh(e2)
  dim(e2)
  
  # ERCC anpassen um zu plotten
  grep("ERCC", e1$Reporter_Group_id, value = T)
  grep("ERCC", e1$Reporter_Group_Name, value = T)
  e1[grep("ERCC", e1$Reporter_Group_id), c("Reporter_Group_id", "Reporter_Group_Name")] = "ERCC"
  e2[grep("ERCC", e2$Reporter_Group_id), c("Reporter_Group_id", "Reporter_Group_Name")] = "ERCC"
  
  if(show_qcplots_individSamplelevel== "from_paramfile") show_qcplots_individSamplelevel_used <- getParam2("show_qcplots_individSamplelevel", myparam = param) else show_qcplots_individSamplelevel_used = show_qcplots_individSamplelevel
  
  
  stopifnot(show_qcplots_individSamplelevel_used %in% c(T, F))
  show_qcplots_individSamplelevel_used
  
  ## ----erccstatus----------------------------------------------------------
  ercc <- genesdetail[ genesdetail$is_ercc, "nuid"]
  ercc
  ercc_vorher_da <- any(ercc %in% rownames(exprs(total_nobkgd_eset)))
  ercc_vorher_da
  ercc_nachher_da <- any(ercc %in% rownames(exprs(total_nobkgd_eset_ql_combat)))
  ercc_nachher_da
  
  ## ----indlev--------------------------------------------------------------
  if(show_qcplots_individSamplelevel_used) {
    plot_housekeeping_vor <- plotGenes3log2("housekeeping", e1m = makePrintDF(e1, "housekeeping")) +
      ggtitle("vor Praeprozessierung")
    plot_housekeeping_vor
    plot_housekeeping_nach <- plotGenes3("housekeeping", e1m = makePrintDF(e2, "housekeeping")) +
      ggtitle("NACH Praeprozessierung")
    plot_housekeeping_nach
    plot_negative_vor <- plotGenes3log2("negative", e1m = makePrintDF(e1, "negative")) +
      ggtitle("vor Praeprozessierung")
    plot_negative_vor
    plot_negative_nach <- plotGenes3("negative", e1m = makePrintDF(e2, "negative")) +
      ggtitle("NACH Praeprozessierung")
    plot_negative_nach
    plot_cy3_hyb_vor <- plotGenes3log2("cy3_hyb", e1m = makePrintDF(e1, "cy3_hyb")) +
      ggtitle("vor Praeprozessierung")
    plot_cy3_hyb_vor
    plot_cy3_hyb_nach <- plotGenes3("cy3_hyb", e1m = makePrintDF(e2, "cy3_hyb")) +
      ggtitle("NACH Praeprozessierung")
    plot_cy3_hyb_nach
    plot_low_stringency_hyb_vor <- plotGenes3log2("low_stringency_hyb", e1m = makePrintDF(e1, "low_stringency_hyb")) +
      ggtitle("vor Praeprozessierung")
    plot_low_stringency_hyb_vor
    plot_low_stringency_hyb_nach <- plotGenes3("low_stringency_hyb", e1m = makePrintDF(e2, "low_stringency_hyb")) +
      ggtitle("NACH Praeprozessierung")
    plot_low_stringency_hyb_nach
    plot_biotin_vor <- plotGenes3log2("biotin", e1m = makePrintDF(e1, "biotin")) +
      ggtitle("vor Praeprozessierung")
    plot_biotin_vor
    plot_biotin_nach <- plotGenes3("biotin", e1m = makePrintDF(e2, "biotin")) +
      ggtitle("NACH Praeprozessierung")
    plot_biotin_nach
    plot_phage_lambda_genome_high_vor <- plotGenes3log2("phage_lambda_genome_high", e1m = makePrintDF(e1, "phage_lambda_genome:high")) + 
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_high_vor
    plot_phage_lambda_genome_high_nach <- plotGenes3("phage_lambda_genome_high", e1m = makePrintDF(e2, "phage_lambda_genome:high")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_high_nach
    plot_phage_lambda_genome_med_vor <- plotGenes3log2("phage_lambda_genome_med", e1m = makePrintDF(e1, "phage_lambda_genome:med")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_med_vor
    plot_phage_lambda_genome_med_nach <- plotGenes3("phage_lambda_genome_med", e1m = makePrintDF(e2, "phage_lambda_genome:med")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_med_nach
    plot_phage_lambda_genome_low_vor <- plotGenes3log2("phage_lambda_genome_low", e1m = makePrintDF(e1, "phage_lambda_genome:low")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_low_vor
    plot_phage_lambda_genome_low_nach <- plotGenes3("phage_lambda_genome_low", e1m = makePrintDF(e2, "phage_lambda_genome:low")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_low_nach
    plot_labeling_vor <- plotGenes3log2("labeling", e1m = makePrintDF(e1, "labeling")) +
      ggtitle("vor Praeprozessierung")
    plot_labeling_vor
    plot_labeling_nach <- plotGenes3("labeling", e1m = makePrintDF(e2, "labeling")) +
      ggtitle("NACH Praeprozessierung")
    plot_labeling_nach
    plot_phage_lambda_genome_pm_vor <- plotGenes3log2("phage_lambda_genome_pm", e1m = makePrintDF(e1, "phage_lambda_genome:pm")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_pm_vor
    plot_phage_lambda_genome_pm_nach <- plotGenes3("phage_lambda_genome_pm", e1m = makePrintDF(e2, "phage_lambda_genome:pm")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_pm_nach
    plot_phage_lambda_genome_mm2_vor <- plotGenes3log2("phage_lambda_genome_mm2", e1m = makePrintDF(e1, "phage_lambda_genome:mm2")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_mm2_vor
    plot_phage_lambda_genome_mm2_nach <- plotGenes3("phage_lambda_genome_mm2", e1m = makePrintDF(e2, "phage_lambda_genome:mm2")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_mm2_nach
    if(ercc_vorher_da) plot_ercc_vor <- plotGenes3log2("ERCC", e1m = makePrintDF(e1, "ERCC")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_vor
    if(ercc_nachher_da)plot_ercc_nach <- plotGenes3("ERCC", e1m = makePrintDF(e2, "ERCC")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_nach
    if(ercc_vorher_da) plot_ercc_0_vor <- plotGenes3log2("ercc_0", e1m = makePrintDF(e1, "ercc_0")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_0_vor
    if(ercc_nachher_da)  plot_ercc_0_nach <- plotGenes3("ercc_0", e1m = makePrintDF(e2, "ercc_0")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da)  plot_ercc_0_nach
    if(ercc_vorher_da)  plot_ercc_1_vor <- plotGenes3log2("ercc_1", e1m = makePrintDF(e1, "ercc_1")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da)  plot_ercc_1_vor
    if(ercc_nachher_da) plot_ercc_1_nach <- plotGenes3("ercc_1", e1m = makePrintDF(e2, "ercc_1")) + 
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_1_nach
    if(ercc_vorher_da) plot_ercc_2_vor <- plotGenes3log2("ercc_2", e1m = makePrintDF(e1, "ercc_2")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da)  plot_ercc_2_vor
    if(ercc_nachher_da) plot_ercc_2_nach <- plotGenes3("ercc_2", e1m = makePrintDF(e2, "ercc_2")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_2_nach
    if(ercc_vorher_da) plot_ercc_3_vor <- plotGenes3log2("ercc_3", e1m = makePrintDF(e1, "ercc_3")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_3_vor
    if(ercc_nachher_da) plot_ercc_3_nach <- plotGenes3("ercc_3", e1m = makePrintDF(e2, "ercc_3")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_3_nach
    if(ercc_vorher_da) plot_ercc_4_vor <- plotGenes3log2("ercc_4", e1m = makePrintDF(e1, "ercc_4")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_4_vor
    if(ercc_nachher_da) plot_ercc_4_nach <- plotGenes3("ercc_4", e1m = makePrintDF(e2, "ercc_4")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_4_nach
    if(ercc_vorher_da) plot_ercc_5_vor <- plotGenes3log2("ercc_5", e1m = makePrintDF(e1, "ercc_5")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_5_vor
    if(ercc_nachher_da) plot_ercc_5_nach <- plotGenes3("ercc_5", e1m = makePrintDF(e2, "ercc_5")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_5_nach
  }
  
  ## ----plotsentrix---------------------------------------------------------
  if(show_qcplots_individSamplelevel_used==F) {
    
    plot_housekeeping_vor <- plotGenes2log2("housekeeping", e1m = makePrintDF(e1, "housekeeping")) +
      ggtitle("vor Praeprozessierung")
    plot_housekeeping_vor
    plot_housekeeping_nach <- plotGenes2("housekeeping", e1m = makePrintDF(e2, "housekeeping")) +
      ggtitle("NACH Praeprozessierung")
    plot_housekeeping_nach
    plot_negative_vor <- plotGenes2log2("negative", e1m = makePrintDF(e1, "negative")) +
      ggtitle("vor Praeprozessierung")
    plot_negative_vor
    plot_negative_nach <- plotGenes2("negative", e1m = makePrintDF(e2, "negative")) +
      ggtitle("NACH Praeprozessierung")
    plot_negative_nach
    plot_cy3_hyb_vor <- plotGenes2log2("cy3_hyb", e1m = makePrintDF(e1, "cy3_hyb")) +
      ggtitle("vor Praeprozessierung")
    plot_cy3_hyb_vor
    plot_cy3_hyb_nach <- plotGenes2("cy3_hyb", e1m = makePrintDF(e2, "cy3_hyb")) +
      ggtitle("NACH Praeprozessierung")
    plot_cy3_hyb_nach
    plot_low_stringency_hyb_vor <- plotGenes2log2("low_stringency_hyb", e1m = makePrintDF(e1, "low_stringency_hyb")) +
      ggtitle("vor Praeprozessierung")
    plot_low_stringency_hyb_vor
    plot_low_stringency_hyb_nach <- plotGenes2("low_stringency_hyb", e1m = makePrintDF(e2, "low_stringency_hyb")) +
      ggtitle("NACH Praeprozessierung")
    plot_low_stringency_hyb_nach
    plot_biotin_vor <- plotGenes2log2("biotin", e1m = makePrintDF(e1, "biotin")) +
      ggtitle("vor Praeprozessierung")
    plot_biotin_vor
    plot_biotin_nach <- plotGenes2("biotin", e1m = makePrintDF(e2, "biotin")) +
      ggtitle("NACH Praeprozessierung")
    plot_biotin_nach
    plot_phage_lambda_genome_high_vor <- plotGenes2log2("phage_lambda_genome_high", e1m = makePrintDF(e1, "phage_lambda_genome:high")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_high_vor
    plot_phage_lambda_genome_high_nach <- plotGenes2("phage_lambda_genome_high", e1m = makePrintDF(e2, "phage_lambda_genome:high")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_high_nach
    plot_phage_lambda_genome_med_vor <- plotGenes2log2("phage_lambda_genome_med", e1m = makePrintDF(e1, "phage_lambda_genome:med")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_med_vor
    plot_phage_lambda_genome_med_nach <- plotGenes2("phage_lambda_genome_med", e1m = makePrintDF(e2, "phage_lambda_genome:med")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_med_nach
    plot_phage_lambda_genome_low_vor <- plotGenes2log2("phage_lambda_genome_low", e1m = makePrintDF(e1, "phage_lambda_genome:low")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_low_vor
    plot_phage_lambda_genome_low_nach <- plotGenes2("phage_lambda_genome_low", e1m = makePrintDF(e2, "phage_lambda_genome:low")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_low_nach
    plot_labeling_vor <- plotGenes2log2("labeling", e1m = makePrintDF(e1, "labeling")) + ggtitle("vor Praeprozessierung")
    plot_labeling_vor
    plot_labeling_nach <- plotGenes2("labeling", e1m = makePrintDF(e2, "labeling")) + ggtitle("NACH Praeprozessierung")
    plot_labeling_nach
    plot_phage_lambda_genome_pm_vor <- plotGenes2log2("phage_lambda_genome_pm", e1m = makePrintDF(e1, "phage_lambda_genome:pm")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_pm_vor
    plot_phage_lambda_genome_pm_nach <- plotGenes2("phage_lambda_genome_pm", e1m = makePrintDF(e2, "phage_lambda_genome:pm")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_pm_nach
    plot_phage_lambda_genome_mm2_vor <- plotGenes2log2("phage_lambda_genome_mm2", e1m = makePrintDF(e1, "phage_lambda_genome:mm2")) +
      ggtitle("vor Praeprozessierung")
    plot_phage_lambda_genome_mm2_vor
    plot_phage_lambda_genome_mm2_nach <- plotGenes2("phage_lambda_genome_mm2", e1m = makePrintDF(e2, "phage_lambda_genome:mm2")) +
      ggtitle("NACH Praeprozessierung")
    plot_phage_lambda_genome_mm2_nach
    if(ercc_vorher_da) plot_ercc_vor <- plotGenes2log2("ERCC", e1m = makePrintDF(e1, "ERCC")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_vor
    if(ercc_nachher_da)  plot_ercc_nach <- plotGenes2("ERCC", e1m = makePrintDF(e2, "ERCC")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da)  plot_ercc_nach
    if(ercc_vorher_da) plot_ercc_0_vor <- plotGenes2log2("ercc_0", e1m = makePrintDF(e1, "ercc_0")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_0_vor
    if(ercc_nachher_da) plot_ercc_0_nach <- plotGenes2("ercc_0", e1m = makePrintDF(e2, "ercc_0")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_0_nach
    if(ercc_vorher_da)plot_ercc_1_vor <- plotGenes2log2("ercc_1", e1m = makePrintDF(e1, "ercc_1")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_1_vor
    if(ercc_nachher_da) plot_ercc_1_nach <- plotGenes2("ercc_1", e1m = makePrintDF(e2, "ercc_1")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_1_nach
    if(ercc_vorher_da) plot_ercc_2_vor <- plotGenes2log2("ercc_2", e1m = makePrintDF(e1, "ercc_2")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_2_vor
    if(ercc_nachher_da) plot_ercc_2_nach <- plotGenes2("ercc_2", e1m = makePrintDF(e2, "ercc_2")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_2_nach
    if(ercc_vorher_da) plot_ercc_3_vor <- plotGenes2log2("ercc_3", e1m = makePrintDF(e1, "ercc_3")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da)  plot_ercc_3_vor
    if(ercc_nachher_da)  plot_ercc_3_nach <- plotGenes2("ercc_3", e1m = makePrintDF(e2, "ercc_3")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_3_nach
    if(ercc_vorher_da) plot_ercc_4_vor <- plotGenes2log2("ercc_4", e1m = makePrintDF(e1, "ercc_4")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_4_vor
    if(ercc_nachher_da) plot_ercc_4_nach <- plotGenes2("ercc_4", e1m = makePrintDF(e2, "ercc_4")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_4_nach
    if(ercc_vorher_da) plot_ercc_5_vor <- plotGenes2log2("ercc_5", e1m = makePrintDF(e1, "ercc_5")) +
      ggtitle("vor Praeprozessierung")
    if(ercc_vorher_da) plot_ercc_5_vor
    if(ercc_nachher_da) plot_ercc_5_nach <- plotGenes2("ercc_5", e1m = makePrintDF(e2, "ercc_5")) +
      ggtitle("NACH Praeprozessierung")
    if(ercc_nachher_da) plot_ercc_5_nach
  }
  
  ## ----repamong------------------------------------------------------------
  # categspalte bauen
  hh(e1)
  makePrintAllDF <- function (e1, attrib = sample_overview_l10instudy) {
    
    # unique(e1[order(e1$Reporter_Group_Name),1:2])
    e1$plotcateg <- e1$Reporter_Group_id
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome"] <- "biotin"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "thrB"] <- "labeling"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "permuted_negative"] <- "negative"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome:high"] <- "hybrid_high"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome:med"] <- "hybrid_med"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome:low"] <- "hybrid_low"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome:mm2"] <- "string_mm"
    e1$plotcateg[is.na(e1$plotcateg) == F & e1$plotcateg == "phage_lambda_genome:pm"] <- "string_pm"
    plotcategories <- unique(na.omit(e1$plotcateg))
    plotcategories
    e1c <- data.table::data.table(e1[is.na(e1$plotcateg) == F, names(e1) %nin% c('Reporter_Group_id',
                                                                                 'Reporter_Group_Name',
                                                                                 "subgroup",
                                                                                 'ilmn')])
    myplot <- melt(e1c,id.vars = "plotcateg")
    myplot[, subgroup := attrib[match_hk(variable, attrib$new_ID), "subgroup"]]
    myplot[, mytable(subgroup, doprint = F)]
    myplot
  }
  e1a <- makePrintAllDF(e1)
  e1a
  plot_allcateg_vor <- ggplot2::ggplot(data.frame(e1a),
                                       ggplot2::aes(plotcateg,
                                                    log2(value),
                                                    fill = plotcateg) ) +
    ggplot2::geom_violin(scale = 'width') +
    ggplot2::theme(axis.text.x = element_text(angle = -30, hjust = 0, size = 12), axis.text.y = element_text(size = 12)) +
    ggplot2::guides(fill = F)  +
    ggtitle("log2 expression values BEFORE preprocess") +
    ggplot2::facet_grid( subgroup ~., scales = "free_x")
  plot_allcateg_vor
  e2a <- makePrintAllDF(e2)
  e2a
  plot_allcateg_nach <- ggplot2::ggplot(data.frame(e2a),
                                        ggplot2::aes(plotcateg,
                                                     value,
                                                     fill = plotcateg) ) +
    ggplot2::geom_violin(scale = 'width') +
    ggplot2::theme(axis.text.x = element_text(angle = -30, hjust = 0, size = 12), axis.text.y = element_text(size = 12)) +
    ggplot2::guides(fill = F) +
    ggtitle("expression values AFTER preprocess") +
    ggplot2::facet_grid( subgroup ~., scales = "free_x")
  plot_allcateg_nach
  
  # einzeln plotten housekeeping
  head(ilmnAnnot014allgInfos)
  
  makeIndPlot <- function(e1, gruppe = "housekeeping", attrib = sample_overview_l10instudy) {
    e1_house <- data.table::data.table(e1[is.na(e1$Reporter_Group_id) == F & e1$Reporter_Group_id %in% gruppe, ])
    messvar <- names(e1_house)
    messvar <- messvar[messvar %in% attrib$new_ID]
    e1_house <- melt(e1_house,id.vars = c( 'Reporter_Group_id', 'Reporter_Group_Name',"ilmn"), measure.vars = messvar)
    e1_house$Reporter_Group_id[e1_house$Reporter_Group_id == "phage_lambda_genome:mm2"] = "string_mm"
    e1_house$Reporter_Group_id[e1_house$Reporter_Group_id == "phage_lambda_genome:pm"] = "string_pm"
    
    # ersetzen ilmn durch gene falls bekannt
    e1_house[, ilmn := ifelse(ilmn %in% paste0("ILMN_", ilmnAnnot014allgInfos$ilmn),
                              ilmnAnnot014allgInfos[match_hk(e1_house$ilmn,
                                                             paste0("ILMN_",ilmnAnnot014allgInfos$ilmn)),
                                                    "SymbolReannotated_orgHsEg"], ilmn)]
    ordered_ilmn <- e1_house[order(Reporter_Group_id), unique(ilmn)]
    e1_house[, ilmn := factor(ilmn, levels = ordered_ilmn)]
    e1_house[, subgroup := attrib[match_hk(variable, attrib$new_ID), "subgroup"]]
    e1_house[, mytable(subgroup, doprint = F)]
    e1_house
  }
  e1_detail <- makeIndPlot(e1, gruppe = c("housekeeping", "phage_lambda_genome:mm2", "phage_lambda_genome:pm"))
  plot_details_vor <- ggplot2::ggplot(e1_detail,
                                      ggplot2::aes(ilmn, log2(value),
                                                   fill = Reporter_Group_id)) +
    ggplot2::geom_violin(scale = 'width') +
    ggplot2::theme(axis.text.x = element_text(angle = -30, hjust = 0, size = 12), axis.text.y = element_text(size = 12)) +
    ggtitle("expression values BEFORE preprocess") +
    ggplot2::facet_grid( subgroup ~., scales = "free_x")
  plot_details_vor
  e2_detail <- makeIndPlot(e2, gruppe = c("housekeeping", "phage_lambda_genome:mm2", "phage_lambda_genome:pm")) 
  plot_details_nach <- ggplot2::ggplot(e2_detail,
                                       ggplot2::aes(ilmn, value, fill = Reporter_Group_id  )) +
    ggplot2::geom_violin(scale = 'width') +
    ggplot2::theme(axis.text.x = element_text(angle = -30, hjust = 0, size = 12), axis.text.y = element_text(size = 12)) +
    ggtitle("expression values AFTER preprocess") +
    ggplot2::facet_grid( subgroup ~., scales = "free_x")
  plot_details_nach
  
  ## ----pdfplot-------------------------------------------------------------
  
  vglPlotte = function(plot_biotin_vor, plot_biotin_nach, mywidth = c(1,1)) {
    ylimit = c(min(c(plot_biotin_vor$data[,1], plot_biotin_nach$data[,1]), na.rm = T), max(c(plot_biotin_vor$data[,1], plot_biotin_nach$data[,1]), na.rm = T))
    suppressMessages(gridExtra::grid.arrange(gridExtra::arrangeGrob((plot_biotin_vor + ggplot2::theme(legend.position="top")+ggplot2::guides(fill=guide_legend(ncol=8)) + ggplot2::ggtitle('BEFORE preprocessing')+ ggplot2::ylim(ylimit)), (plot_biotin_nach + ggplot2::scale_y_continuous(labels = function(x) signif(x,5))+ ggplot2::ggtitle('AFTER preprocessing') + ggplot2::theme(legend.position="top")+ggplot2::guides(fill=guide_legend(ncol=8))+ ggplot2::ylim(ylimit)), ncol=2, widths=mywidth)))}
  
  
  
  vglPlotte2 = function(plot_biotin_vor, plot_biotin_nach, mywidth = c(1,2)) {
    suppressMessages(gridExtra::grid.arrange(gridExtra::arrangeGrob((plot_biotin_vor  + ggplot2::ggtitle('BEFORE preprocessing')+ggplot2::guides(fill=guide_legend(ncol=6))+ ggplot2::theme(legend.position="top")), (plot_biotin_nach+ ggplot2::ggtitle('AFTER preprocessing')+ggplot2::guides(fill=guide_legend(ncol=6))+ ggplot2::theme(legend.position="top")), ncol=2, widths=mywidth)))}
  
  
  
  if(showPlots==T) {
    
    
    
    
    
    vglPlotte(plot_biotin_vor, plot_biotin_nach)
    vglPlotte(plot_cy3_hyb_vor, plot_cy3_hyb_nach)
    vglPlotte(plot_housekeeping_vor, plot_housekeeping_nach)
    vglPlotte(plot_labeling_vor, plot_labeling_nach)
    vglPlotte(plot_negative_vor, plot_negative_nach)
    
    vglPlotte(plot_low_stringency_hyb_vor, plot_low_stringency_hyb_nach)
    vglPlotte(plot_phage_lambda_genome_mm2_vor, plot_phage_lambda_genome_mm2_nach)
    vglPlotte(plot_phage_lambda_genome_pm_vor, plot_phage_lambda_genome_pm_nach)
    vglPlotte(plot_phage_lambda_genome_low_vor, plot_phage_lambda_genome_low_nach)
    vglPlotte(plot_phage_lambda_genome_med_vor, plot_phage_lambda_genome_med_nach)
    vglPlotte(plot_phage_lambda_genome_high_vor, plot_phage_lambda_genome_high_nach)
    
    
    if(ercc_vorher_da & ercc_nachher_da){
      vglPlotte(plot_ercc_vor, plot_ercc_nach)
      vglPlotte(plot_ercc_0_vor, plot_ercc_0_nach)
      vglPlotte(plot_ercc_1_vor, plot_ercc_1_nach)
      vglPlotte(plot_ercc_2_vor, plot_ercc_2_nach)
      vglPlotte(plot_ercc_3_vor, plot_ercc_3_nach)
      vglPlotte(plot_ercc_4_vor, plot_ercc_4_nach)
      vglPlotte(plot_ercc_5_vor, plot_ercc_5_nach)
    }
    
    
    ylimits = c(5,17)
    ybreaks = seq(5,17,1)
    vglPlotte2(plot_allcateg_vor + ggplot2::scale_y_continuous(lim=ylimits, breaks = ybreaks), plot_allcateg_nach  + scale_y_continuous(lim=ylimits, breaks =ybreaks), mywidth = c(1,1))
    
    
    vglPlotte2(plot_details_vor +ggplot2::scale_y_continuous(lim=ylimits, breaks = ybreaks), plot_details_nach  + scale_y_continuous(lim=ylimits, breaks =ybreaks), mywidth = c(1,1))
  }
  
  
  
  # 
  # plot_housekeeping_vor
  # plot_housekeeping_nach
  # plot_negative_vor
  # plot_negative_nach
  # plot_cy3_hyb_vor
  # plot_cy3_hyb_nach
  # plot_low_stringency_hyb_vor
  # plot_low_stringency_hyb_nach
  # plot_biotin_vor
  # plot_biotin_nach
  # plot_phage_lambda_genome_high_vor
  # plot_phage_lambda_genome_high_nach
  # plot_phage_lambda_genome_med_vor
  # plot_phage_lambda_genome_med_nach
  # plot_phage_lambda_genome_low_vor
  # plot_phage_lambda_genome_low_nach
  # plot_labeling_vor
  # plot_labeling_nach
  # plot_phage_lambda_genome_pm_vor
  # plot_phage_lambda_genome_pm_nach
  # plot_phage_lambda_genome_mm2_vor
  # plot_phage_lambda_genome_mm2_nach
  # if(ercc_vorher_da) plot_ercc_vor
  # if(ercc_nachher_da)plot_ercc_nach
  # if(ercc_vorher_da) plot_ercc_0_vor
  # if(ercc_nachher_da) plot_ercc_0_nach
  # if(ercc_vorher_da) plot_ercc_1_vor
  # if(ercc_nachher_da) plot_ercc_1_nach
  # if(ercc_vorher_da) plot_ercc_2_vor
  # if(ercc_nachher_da) plot_ercc_2_nach
  # if(ercc_vorher_da) plot_ercc_3_vor
  # if(ercc_nachher_da) plot_ercc_3_nach
  # if(ercc_vorher_da) plot_ercc_4_vor
  # if(ercc_nachher_da) plot_ercc_4_nach
  # if(ercc_vorher_da)  plot_ercc_5_vor
  # if(ercc_nachher_da) plot_ercc_5_nach
  
  # plot_goodexpression_nach
  # plot_allcateg_vor 
  # plot_allcateg_nach
  # plot_details_vor 
  # plot_details_nach
  
  
  plots <- grep('^plot_', ls(), value = T)
  
  
  fordoku <-   c(plots, "ercc_vorher_da", "ercc_nachher_da")
  
  stopifnot(sum(duplicated(fordoku))==0)
  
  
  ht12object$dokuobjects_visualizePreprocessing = lapply(fordoku, function(x) get(x))
  
  
  names(ht12object$dokuobjects_visualizePreprocessing) = fordoku
  
  
  ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
  ht12object$history
  ht12object
  
}