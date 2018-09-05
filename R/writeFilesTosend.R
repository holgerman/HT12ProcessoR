#' @title Writes File used for follow up analyses in a folder
#'
#' @description Save the preprocessed data in a folder for further using filenames specified as parameters. Parameter `file_renaming_samples_tosend` can specify a filename where a column named `oldname` and `newname` is provided in order to rename sample. Column oldname must match the names of the samples in column new_ID  in the `$chipsamples` data.frame included in the HT12 object. Can be empty if no renaming is required

#' @param ht12object A list object of class HT12prepro created with function checkExtractChipsamples
#' @param paramfile Path to the file specifying parameters
#' @param datafolder_results directory of the resulting preprocessed data
#'
#' @param file_sampleannot_final.txt default is	sampleannot_HT12v4.txt Filename for sample-related attributes after finished preprocessing
#' @param file_probeannot_final.txt Filename for expression probe-related attributes after finished preprocessing
#' @param file_annot_final.xlsx Filename for an excel file combining 'file_sampleannot_final.txt' and 'file_probeannot_final.txt'.
#' @param file_final_expression_set Filename of expression set  after finished preprocessing
#' @param file_final_expression_matrix Filename of expression matrix as R object after finished preprocessing
#' @param file_all_transcripts_good_incl_remapping_ok Filename where probes are listed that are classified as expressed, not over-inflated for batch effects and classified at least 'good' according to remapping on the human genome as described in Nucleic Acids Res. Januar 2010;38(3):e17. in all user-provided subgroups
#' @param file_final_expression_matrix_allProbes.txt Filename of expression matrix as tab delimited text file after finished preprocessing
#' @param file_renaming_samples_tosend Filename where a column named `oldname` and `newname` is provided in order to rename samples, column oldname must match the names of the sample in column new_ID  in the `chipsamples` data.frame. Can be empty if no renaming is required


#' @return Written files defined in the parameters
#' @import data.table
#' @export

## debug
# paramfile = myparamfile
# ht12object =  prepro_ht12
# file_sampleannot_final.txt = "from_paramfile"
# file_probeannot_final.txt = "from_paramfile"
# file_annot_final.xlsx = "from_paramfile"
# file_final_expression_set = "from_paramfile"
# file_final_expression_matrix = "from_paramfile"
# file_all_transcripts_good_incl_remapping_ok = "from_paramfile"
# file_final_expression_matrix_allProbes.txt = "from_paramfile"
# file_renaming_samples_tosend = "from_paramfile"





writeFilesTosend = function(ht12object,paramfile = NULL,file_sampleannot_final.txt = "from_paramfile",file_probeannot_final.txt= "from_paramfile",file_annot_final.xlsx= "from_paramfile",file_final_expression_set= "from_paramfile",file_final_expression_matrix= "from_paramfile",file_all_transcripts_good_incl_remapping_ok= "from_paramfile",file_final_expression_matrix_allProbes.txt= "from_paramfile",file_renaming_samples_tosend= "from_paramfile") {


### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

  myparameters = match.call()
  showVennplots = F


  #laden parameter
  if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)



sample_overview_l10 <- ht12object$chipsamples
mytable(sample_overview_l10$in_study)
sample_overview_l10instudy <- sample_overview_l10[sample_overview_l10$in_study, ]
dim(sample_overview_l10)
dim(sample_overview_l10instudy)
table(table(sample_overview_l10instudy$new_ID))
if(length(table(table(sample_overview_l10instudy$new_ID))) != 1)
  stop("Modify Code! Script assumes unique IDs!")

# laden annotation probes

genesdetail <- ht12object$genesdetail
mytable(genesdetail$expressed)
ht(genesdetail, 2)
mytable(stringr::str_sub(genesdetail$ilmn, 1, 4))
showClassDF(genesdetail)

#laden expressionsets nach combat

total_nobkgd_eset_ql_combat = ht12object$total_nobkgd_eset_ql_combat
total_nobkgd_eset_ql_combat
goodind <- sample_overview_l10instudy$new_ID
# str(goodind)

# laden ilmnAnnot
ht(ilmnAnnot014allgInfos,1)
## ----narrowdown----------------------------------------------------------
qlist577 <- venn2(sampleNames(total_nobkgd_eset_ql_combat), goodind, plotte = showVennplots)
# str(qlist577)
total_nobkgd_eset_ql_combat = total_nobkgd_eset_ql_combat[,goodind]
total_nobkgd_eset_ql_combat
ht(ilmnAnnot014allgInfos, 2)
col2use <- c("ilmn", "nuid", "Species",
            "Source", "Search_Key", "Transcript",
            "ILMN_Gene", "Source_Reference_ID",
            "RefSeq_ID", "Unigene_ID", "Entrez_Gene_ID",
            "GI", "Accession", "Symbol",
            "Protein_Product", "Array_Address_Id",
            "Probe_Type", "Probe_Start", "Probe_Sequence",
            "Chromosome", "Probe_Chr_Orientation",
            "Probe_Coordinates", "Cytoband", "Definition",
            "Ontology_Component", "Ontology_Process", "Ontology_Function",
            "Synonyms", "Obsolete_Probe_Id", "EntrezReannotated",
            "ProbeQuality", "badprobe_dunning", "CodingZone",
            "GenomicLocation", "SecondMatches", "OtherGenomicMatches",
            "RepeatMask", "OverlappingSNP", "RefseqIdReannotated_orgHsEg",
            "SymbolReannotated_orgHsEg", "EnsemblReannotated_orgHsEg", "tr_chr",
            "tr_strand", "tr_start", "tr_ende", "mappings")
doku_excluded_allgPRobeannot <- c("is_control", "ArrayAddress", "ProbeSequence",
                                 "ReporterGroupName", "ReporterGroupID", "SymbolReannotated",
                                 "GenomicLocation_chr", "GenomicLocation_strand", "GenomicLocation_start",
                                 "GenomicLocation_ende","snp_interferes", "validsnp_interferes")
qlist33 <- venn3(col2use, doku_excluded_allgPRobeannot, names(ilmnAnnot014allgInfos),plotte= showVennplots )
ilmnAnnotallg <- ilmnAnnot014allgInfos[, col2use]

## ----narrowdowninf-------------------------------------------------------
ht(genesdetail, 1)
subsetcols <- grep(paste(c("present",
                           "meane",
                           "meanall",
                           "expressed",
                           "MeanRatioExpressed",
                           "MaxRatioExpressed",
                           "goodexpressedprobe"),
                         collapse = "|"),
                   names(genesdetail), value = T)
subsetcols
anovacols <- grep("all_pval", names(genesdetail), value = T)
anovacols
col2use2 <- c("nuid", "is_purecontrol", "is_ercc", subsetcols, "spearmancor_combat", anovacols)
excluded2 <- c("ilmn")
qlist55 <- venn3(excluded2, col2use2, names(genesdetail), plotte =  showVennplots)
# str(qlist55)
genesdetail2 <- genesdetail[, col2use2]
qlist2 <- venn2(names(ilmnAnnotallg), names(genesdetail2), plotte =  showVennplots)

## ----zusammen------------------------------------------------------------
probeannot <- merge(genesdetail2, ilmnAnnotallg, by = "nuid", all = T, sort = F)
probeannot <- moveColFront(probeannot, "ilmn")
probeannot$ilmn <- paste0("ILMN_", probeannot$ilmn)
probeannot$ilmn[probeannot$ilmn == "ILMN_NA"] <- probeannot$nuid[probeannot$ilmn == "ILMN_NA"]
check2 <- table(stringr::str_length(probeannot$ilmn))
check2
stopifnot(identical(names(check2), "12"))
ht(data.frame(t(ht(probeannot, 1))))

## ----probeclassdrann-----------------------------------------------------
# controllinformationen auslesen und anhuebschen
head(all_con)
qlist3 <- venn2(all_con$Probe_Id, probeannot$ilmn, plotte = showVennplots)
#data.table::setDT(all_con)
 all_con = data.table::data.table(data.frame(all_con))
# stopifnot(length(qlist3$q2)==0)
all_con[,Probe_Class := apply(all_con[, c("Reporter_Group_Name","Reporter_Group_id"), with = F], 1,
                             function(x) paste(unique(x), collapse = " | "))]

# ercc infos noch hochaufgeloester drann
ht(ercc_det, 1)
ercc_det$ercc_gruppe1 <- cut(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.),
                             breaks = c(-3,0,1,2,3,4,5), include.lowest = T )
mytable(ercc_det$ercc_gruppe1)
# boxplot(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.) ~ercc_det$ercc_gruppe1)
ercc_det$ercc_gruppe1short <- paste0("(ercc_", stringr::str_sub(sapply(stringr::str_split(ercc_det$ercc_gruppe1, ","),
                                                              "[",2),1,1), ")")
xtabs(~ ercc_det$ercc_gruppe1short + ercc_det$ercc_gruppe1)
ht(ercc_det, 2)
all_con$ercc_det <- ercc_det[ match_hk(all_con$Reporter_Group_Name,
                                       ercc_det$ERCC.ID),
                              'ercc_gruppe1short']
all_con$Probe_Class <- ifelse(is.na(all_con$ercc_det),
                              all_con$Probe_Class,
                              all_con$ercc_det)
mytable(all_con$Probe_Class)
coninfo2 <- plyr::ddply(all_con, "Probe_Id", function(df){
  # df = all_con[ all_con$Probe_Id==all_con$Probe_Id[1],]
  data.frame(Probe_Class = paste(unique(df$Probe_Class), collapse = " | "))
})
mytable(coninfo2$Probe_Class)
coninfo2$Probe_Class[coninfo2$Probe_Class == "biotin | phage_lambda_genome"] <- "(biotin)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "cy3_hyb | phage_lambda_genome:high | low_stringency_hyb | phage_lambda_genome:pm"] <- "(high conc. & perfect match)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "cy3_hyb | phage_lambda_genome:low"] <- "(low conc.)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "cy3_hyb | phage_lambda_genome:med | low_stringency_hyb | phage_lambda_genome:pm"] <- "(medium conc. & perfect match)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "housekeeping"] <- "(housekeeping)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "labeling | thrB"] <- "(labeling)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "low_stringency_hyb | phage_lambda_genome:mm2"] <- "(mismatch)"
coninfo2$Probe_Class[coninfo2$Probe_Class == "negative | permuted_negative"] <- "(negative)"
mytable(coninfo2$Probe_Class)
probeannot$Probe_Class <- ifelse(probeannot$nuid %in% ilmnAnnotallg$nuid, "transcript ", "control ")
probeannot$Probe_Class2add <- coninfo2[match_hk(probeannot$ilmn, coninfo2$Probe_Id), 'Probe_Class']
probeannot$Probe_Class2add[is.na(probeannot$Probe_Class2add)] <- ""
mytable(probeannot$Probe_Class2add)
probeannot$Probe_Class[grep('ercc', probeannot$Probe_Class2add)] <- "control "
probeannot$Probe_Class <- paste0(probeannot$Probe_Class, probeannot$Probe_Class2add)
mytable(probeannot$Probe_Class)
probeannot$Probe_Class2add = NULL

## ----samplannot----------------------------------------------------------
ht(sample_overview_l10, 1)
a <- names(sample_overview_l10)
relnames <- c("old_ID", "new_ID",  "subgroup",
             "nobkgd_f", "con_f", "sample_f", "fileset_id",
             "Sample.Group", "Sentrix.Barcode", "Sample.Section",
             "Detected.Genes..0.01.", "Detected.Genes..0.01._relIQR",
             "Detected.Genes..0.05.", "Signal.Average",
             "Signal.P05", "Signal.P25", "Signal.P50",
             "Signal.P75", "Signal.P95", "BIOTIN", "CY3_HYB",
             "HOUSEKEEPING", "LABELING", "LOW_STRINGENCY_HYB",
             "NEGATIVE..background.", "Noise", "in_study",
             "reason4exclusion", "mahal", "mahal_relIQR",
             "euklid","euklid_relIQR")
doku_excluded_samplecols <- c("paste_id", "different_cluster", "hybrid_low_s07",
                              "hybrid_med_s07", "hybrid_high_s07", "string_pm_s07",
                              "string_mm_s07", "biotin_s07", "housekeeping_s07",
                              "labeling_s07", "ercc_s07", "negative_s07",
                              "samplemean_s07", "p05_s07", "p25_s07",
                              "p50_s07", "p75_s07", "p95_s07")
qlist44 <- venn3(relnames, doku_excluded_samplecols, names(sample_overview_l10), plotte = showVennplots)
sampleinfo <- sample_overview_l10[, relnames]
xtabs(~ sampleinfo$reason4exclusion + is.na(sampleinfo$mahal), exclude = NULL, na.action = "na.pass")

## ----tidyup--------------------------------------------------------------
total_nobkgd_eset_ql_combat
qlist3 <- venn2(featureNames(total_nobkgd_eset_ql_combat), probeannot$nuid, plotte = showVennplots)
featureNames(total_nobkgd_eset_ql_combat) <- probeannot[match_hk(featureNames(total_nobkgd_eset_ql_combat), probeannot$nuid), "ilmn"]

## ----gooditrans----------------------------------------------------------
par(mfrow = c(1,1))
existing_objekts = ls()


probeannot$is_neverNA <- probeannot$ilmn %in% genesdetail$ilmn
mytable(probeannot$is_neverNA)
probeannot <- moveColFront(probeannot, c("ilmn", "nuid", "is_neverNA"))
hh(probeannot, 9)
expressedprobecols <- grep(paste(c("^expressed"), collapse = "|"), names(genesdetail), value = T)
expressedprobecols
showme <- data.frame(eval(expr = parse(text = paste0("with(probeannot, xtabs_hk(~is_purecontrol + is_ercc +",
                                                     paste(expressedprobecols, collapse = "+"), "))"))))
showme[order(showme$is_purecontrol, showme$is_ercc,showme$Freq), ]
probeannot$expressedprobe_allsubgroups = apply(probeannot[, expressedprobecols,
                                                          drop = F], 1,
                                               function(x) all(x == T))
mytable(probeannot$expressedprobe_allsubgroups)
goodexpressedprobecols <- grep(paste(c("goodexpressedprobe"), collapse = "|"),
                               names(genesdetail), value = T)
goodexpressedprobecols
showme <- data.frame(eval(expr = parse(text = paste0("with(probeannot, xtabs_hk(~is_purecontrol + is_ercc +",
                                                     paste(goodexpressedprobecols, collapse = "+"), "))"))))
showme[order(showme$is_purecontrol, showme$is_ercc,showme$Freq), ]
probeannot$goodexpressedprobe_allsubgroups <- apply(probeannot[, goodexpressedprobecols, drop = F], 1,
                                                    function(x) all(x == T))
mytable(probeannot$goodexpressedprobe_allsubgroups)
probeannot$goodprobe_remapping <- probeannot$badprobe_dunning == F
goodfilter <- probeannot$goodexpressedprobe_allsubgroups & probeannot$goodprobe_remapping & probeannot$is_neverNA
mytable(goodfilter)
hh(probeannot)
good_ilmn <- unique(probeannot[ goodfilter, "ilmn"])
# str(good_ilmn)
good_entrez <- unique(na.omit(probeannot[goodfilter, "EntrezReannotated"]))
# str(good_entrez)

## ----ordering------------------------------------------------------------
bb <- names(probeannot)
anova_after_bfcols <- grep("all_pval_after[a-zA-Z_].*_bf",bb, value = T)
anovacols_sub = setdiff(anovacols,anova_after_bfcols)
relcolsprobe = c("ilmn", "nuid", "Array_Address_Id", "Probe_Class",
                 "Probe_Type", "Species", "Source", "Search_Key",
                 "Transcript", "ILMN_Gene", "Source_Reference_ID",
                 "RefSeq_ID", "Unigene_ID", "Entrez_Gene_ID", "GI",
                 "Accession", "Symbol", "Protein_Product",
                 "Probe_Start", "Probe_Sequence", "Chromosome",
                 "Probe_Chr_Orientation", "Probe_Coordinates",
                 "Cytoband", "Definition", "Ontology_Component",
                 "Ontology_Process", "Ontology_Function", "Synonyms",
                 "Obsolete_Probe_Id", "EntrezReannotated",
                 "ProbeQuality", "CodingZone", "GenomicLocation",
                 "SecondMatches", "OtherGenomicMatches",
                 "RepeatMask","OverlappingSNP",
                 "RefseqIdReannotated_orgHsEg", "SymbolReannotated_orgHsEg",
                 "EnsemblReannotated_orgHsEg", "tr_chr", "tr_strand",
                 "tr_start", "tr_ende", "mappings",
                 grep("present", bb, value = T),
                 grep("meane", bb, value = T),
                 grep("meanall", bb, value = T),
                 grep("MeanRatioExpressed", bb, value = T),
                 grep("MaxRatioExpressed", bb, value = T),
                 "spearmancor_combat",
                 anovacols_sub,
                 "is_neverNA",
                 "is_ercc","is_purecontrol",
                 grep( "^expressed", bb, value = T),
                 anova_after_bfcols,
                 setdiff(grep("goodexpressedprobe", bb, value = T),
                         "goodexpressedprobe_allsubgroups"),
                 "goodexpressedprobe_allsubgroups", "badprobe_dunning", 'goodprobe_remapping')
sort(relcolsprobe)
check23 <- venn2(relcolsprobe, names(probeannot), plotte =  showVennplots)
check23
stopifnot(length(c(check23$q2, check23$q3))==0)
relcolsprobe

probeannot <- probeannot[,relcolsprobe]
for(mysubgroup in unique(sampleinfo$subgroup )) {
  # mysubgroup = unique(sampleinfo$subgroup )[1]
  newvar <- paste('perfectprobe', mysubgroup, sep = "_")
  probeannot$tempvar <- probeannot[, paste0("goodexpressedprobe_",mysubgroup)] &
    (probeannot$badprobe_dunning==F)
  probeannot <- rename(probeannot,c(tempvar = newvar))
}
probeannot$perfectprobe_allsubgroups <- probeannot$ilmn %in% good_ilmn
mytable(probeannot$perfectprobe_allsubgroups)

## ----speichern2----------------------------------------------------------


## ----annotis-------------------------------------------------------------
## annotation der samples
ht(annot_sample)

qlist_check3 = venn2(annot_sample$Column, names(sampleinfo), plotte =  showVennplots)
qlist_check3
stopifnot(length(qlist_check3$q3)==0)
rownames(annot_sample) <- annot_sample$Column
annot_sample <- annot_sample[names(sampleinfo), ]
rownames(annot_sample) <- NULL

# nun noch die probes
head(annot_probes_ori)
annot_probes = annot_probes_ori
ht(annot_probes)

# read_excel("../../1406_ge_lifea1/sent/159514_markus_charge/CONFIDENTIAL_s12_sampleUNDprobeannot_HT12v4_annotated.xlsx", 4)
# expand annotationfor subgroup specific parameters
givenparam_subgroupspecific <- c("present",
                                 "expressed",
                                 "goodexpressedprobe",
                                 "meane",
                                 "MeanRatioExpressed",
                                 "MaxRatioExpressed",
                                 "meanall",
                                 "perfectprobe")
newparam_subgroupspecific <- data.frame(expand.grid(givenparam_subgroupspecific,
                                                    unique(sampleinfo$subgroup)))
newparam_subgroupspecific$newparam <- paste(newparam_subgroupspecific[, 1],
                                            newparam_subgroupspecific[, 2],
                                            sep = "_")
names(newparam_subgroupspecific) <- c('Column', 'subgroup', 'newparam')
givenannot_subgroupspecific <- annot_probes[annot_probes$Column %in% givenparam_subgroupspecific, ]
toaddannot_subgroupspecific <- merge(newparam_subgroupspecific, givenannot_subgroupspecific, by = "Column")
toaddannot_subgroupspecific$Description <- paste(toaddannot_subgroupspecific$Description,
                                                 paste0("(for subgroup ", toaddannot_subgroupspecific$subgroup, ")"))
toaddannot_subgroupspecific$Column <- paste(toaddannot_subgroupspecific$Column,
                                            toaddannot_subgroupspecific$subgroup,
                                            sep = "_")
names(toaddannot_subgroupspecific)
names(annot_probes)
toaddannot_subgroupspecific <- toaddannot_subgroupspecific[, names(annot_probes)]
annot_probes <- rbind(annot_probes, toaddannot_subgroupspecific)
qlist_check4 <- venn2(annot_probes$Column, names(probeannot) , plotte =  showVennplots)
qlist_check4
stopifnot(length(qlist_check4$q3) == 0)

# reihenfolge anpassen
annot_probes[allDuplicatedEntries(annot_probes$Column), ]
rownames(annot_probes) <- annot_probes$Column
annot_probes <- annot_probes[names(probeannot), ]
rownames(annot_probes) <- NULL

## ----renamefinal---------------------------------------------------------
if(file_renaming_samples_tosend== "from_paramfile") file_renaming_samples_tosend_fn <- getParam2("file_renaming_samples_tosend", myparam = param) else file_renaming_samples_tosend_fn = file_renaming_samples_tosend

if(file_renaming_samples_tosend_fn == "" )
  renaming <- data.table(oldname = sample_overview_l10$new_ID,
                         newname = sample_overview_l10$new_ID) else {
file_renaming_samples_tosend <- fread(file_renaming_samples_tosend_fn)
renaming <- data.table(oldname = file_renaming_samples_tosend[,oldname],
                       newname = file_renaming_samples_tosend[, newname])
}
renaming
showNA(renaming)
if(sum(showNA(renaming)$NAs)!=0) stop("NAs in renamingfile specified in parameter `file_renaming_samples_tosend` not allowed.... stopping")
if(any(sampleNames(total_nobkgd_eset_ql_combat) %nin% renaming$oldname)) {
  print(sampleNames(total_nobkgd_eset_ql_combat)[sampleNames(total_nobkgd_eset_ql_combat) %nin% renaming$oldname])
  print(ht(renaming))
  stop("Not all IDs from column 'new_ID' of the datafram from slot $chipsamples found in the provided renaming file.... stopping")
  }

# umbenennen und fehlende
sampleinfo$finalID <- renaming[match_hk(sampleinfo$new_ID, renaming$oldname), newname]
sampleinfo <- moveColFront(sampleinfo, "finalID")
hh(sampleinfo)
stopifnot(all(is.na(sampleinfo$finalID)) == F)
class(total_nobkgd_eset_ql_combat)
sampleNames(total_nobkgd_eset_ql_combat) = renaming[ match_hk(sampleNames(total_nobkgd_eset_ql_combat), renaming$oldname), newname]
pData(total_nobkgd_eset_ql_combat)$sampleID = renaming[ match_hk(pData(total_nobkgd_eset_ql_combat)$sampleID, renaming$oldname), newname]

ht(pData(total_nobkgd_eset_ql_combat),2)

showNA(pData(total_nobkgd_eset_ql_combat))

## ----buildmissing objects------------------------------------------------
# speichern weiter vorbereiten
# matrices
gx <- exprs(total_nobkgd_eset_ql_combat)
hh(gx)
gx_spss <- gx
class(gx_spss)
rownames_gx_spss <- rownames(gx_spss)
gx_spss <- data.frame(gx_spss, check.names = F)
#setDT(gx_spss)
gx_spss = data.table::data.table(data.frame(gx_spss))

gx_spss <- cbind(rownames_gx_spss, gx_spss)
setnames(gx_spss, 'rownames_gx_spss', "ilmn_id")
hh(gx_spss)

## ----save----------------------------------------------------------------
if(file_sampleannot_final.txt== "from_paramfile") file_sampleannot_final <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_sampleannot_final.txt", myparam = param)) else file_sampleannot_final = file_sampleannot_final.txt

if(file_probeannot_final.txt== "from_paramfile") file_probeannot_final <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_probeannot_final.txt", myparam = param)) else file_probeannot_final.txt = file_probeannot_final

if(file_annot_final.xlsx== "from_paramfile") file_annot_final <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_annot_final.xlsx", myparam = param)) else file_annot_final = file_annot_final.xlsx

if(file_final_expression_set== "from_paramfile") file_final_expression_set <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_final_expression_set", myparam = param)) else file_final_expression_set = file_final_expression_set

if(file_final_expression_matrix== "from_paramfile") file_final_expression_matrix <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_final_expression_matrix", myparam = param)) else file_final_expression_matrix = file_final_expression_matrix

if(file_all_transcripts_good_incl_remapping_ok== "from_paramfile") file_all_transcripts_good_incl_remapping_ok <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_all_transcripts_good_incl_remapping_ok", myparam = param)) else file_all_transcripts_good_incl_remapping_ok = file_all_transcripts_good_incl_remapping_ok

if(file_final_expression_matrix_allProbes.txt== "from_paramfile") file_final_expression_matrix_allProbes.txt <- paste0(getParam2("datafolder_results",myparam = param), "/", getParam2("file_final_expression_matrix_allProbes.txt", myparam = param)) else file_final_expression_matrix_allProbes.txt = file_final_expression_matrix_allProbes.txt


message("Writing annotation files as txt...")
write.delim(good_ilmn, file_all_transcripts_good_incl_remapping_ok,
            writeColnames = F,
            createDir = T)
dim(sampleinfo)
write.delim(sampleinfo, file_sampleannot_final, createDir = T)
dim(probeannot)
write.delim(probeannot, file_probeannot_final, createDir = T)
eset_preproc <- total_nobkgd_eset_ql_combat
eset_preproc
message("Writing expression data files..")
save(eset_preproc, file = file_final_expression_set)
dim(gx)
save(gx, file = file_final_expression_matrix)
dim(gx_spss)
write.delim(gx_spss, file_final_expression_matrix_allProbes.txt, createDir = T)
total_nobkgd_eset_ql_combat[good_ilmn,]

message("Writing annotation files as xlsx...")
WriteXLS::WriteXLS(c('sampleinfo','annot_sample', 'probeannot', 'annot_probes'),
         ExcelFileName=file_annot_final,
         SheetNames= c("sample annot",
                       "README sample annot",
                       "probe annot",
                       "README probe annot"),
         AutoFilter = T,
         BoldHeaderRow = T,
         FreezeRow = 1,
         FreezeCol = 1)

# for doku
fordoku <- grep("doku|good_", ls(), value = T)
fordoku
stopifnot(sum(duplicated(fordoku))==0)


ht12object$dokuobjects_writeFilesTosend = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_writeFilesTosend) = fordoku


ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object$genesdetail = probeannot
ht12object


}


