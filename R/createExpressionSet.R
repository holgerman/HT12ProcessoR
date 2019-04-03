#' Create an expressionset from Illumina HT12v4 data
#'
#' @description This function loads the expression data from several text files created by Illumina GenomeStudio. These files are specified via the slot $chipsamples of an HT12prepro-object from columns "fileset_id", "nobkgd_f", "con_f",and "sample_f". See vignette for an example.
#' @param ht12object A list object of class HT12prepro created with function checkExtractChipsamples()
#' @param paramfile Path to the file specifying parameters
#' @param colseparator Separator splitting numbers into units of 1000. Typically "." or  ",". If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @return A list object of class HT12prepro with a slot with updated  sample-related attributes of the current processing-stage named `$chipsamples`, a slot with the R-object holding expression data used to create the expression set named `$rawdataWOcons_joined`, a slot with the R-object holding control data used to create the expression set named `$rawdataOnlycons_joined`, a slot with the R-object holding control data used to create the expression set named `$rawdataOnlycons_joined`, a slot with expression-set expression data excluding control data named `$all_nobkgd_eset`, a slot with expression-set expression data including control data named `$total_nobkgd_eset`, a slot with the history of the commands named `$history`including expressionset slots
#' @import data.table
#' @export

## debug
# paramfile = myparamfile
# colseparator= "from_paramfile"
# ht12object =  prepro_ht12
# require(data.table)
#
createExpressionSet = function(ht12object,paramfile = NULL, colseparator= "from_paramfile") {
### strings are imported as strings and not as factors
options(stringsAsFactors=FALSE)

myparameters = match.call()
showVennplots = F

# status checken
historie =  ht12object$history$calls
if(any(grepl("checkExtractChipsamples", historie))==F) stop("Function 'checkExtractChipsamples()' has to be run before!")

#laden parameter
if(is.null(paramfile)==F) param <- data.frame(data.table::fread(paramfile))

laengeilmn = dim(allilmn)[1]
if(laengeilmn != 47323) stop("47323 ilmn-IDs aus file 'data/all47323ilmn.rda' wurden nicht richtig geladen")


if("strangebatch" %nin% names(ht12object$chipsamples)) ht12object$chipsamples$strangebatch = "none"

sample_overview_l3_initial = ht12object$chipsamples

dim(sample_overview_l3_initial)
sample_overview_l3 = sample_overview_l3_initial[ sample_overview_l3_initial$in_study,]
dim(sample_overview_l3)
unique(sample_overview_l3$fileset_id)
showClassDF(sample_overview_l3)


## ----createNoBkgdFile----------------------------------------------------

#bauen der fileuebersicht zum laden
filestable <- unique(sample_overview_l3[,c("fileset_id", "nobkgd_f", "con_f", "sample_f")])
showClassDF(filestable)


#nobkgd importieren
id_overview <- c()

#liste aller vorhandenen ilmns einlesen
all_nobkgd <- data.table( allilmn)
setnames(all_nobkgd, "ilmn", "PROBE_ID")

if(colseparator== "from_paramfile") colseparator_used = getParam2("colseparator", myparam = param) else colseparator_used = colseparator


colseparator_used
for(i in 1:length(filestable[,1])){
  #     i = 1
  #alles laden
  message("Import file ",i,  ": " ,filestable$nobkgd_f[i], "...\n")
  nobkgd <- suppressWarnings(data.table::fread(paste(filestable$nobkgd_f[i], sep="/" )))
  names(nobkgd)
  ht(nobkgd,1)

  #expressiosn spalten reduzieren auf ids im table
  #die IDs, von denen ich daten importieren will
  samples_in <- as.character(sample_overview_l3[as.character(sample_overview_l3$nobkgd_f) ==
                                                  as.character(filestable$nobkgd_f[i]), "old_ID"])
  if(length(grep("^X", samples_in)) != 0) stop("irgendwelche sampleIDs fangen  mit 'X' an -> skript modifizieren, weil z.Z. aus bequemlichkeit davon ausgehe, dass nach dataframe import ein beginnendes 'X' vom Import einer spalte herruehrt, die original eine zahl war!!!")
  if(length(samples_in)==0) stop("keine IDs zu importieren - checke skript!")
  check = table(table(unique(samples_in)))
  check
  if(identical(names(check),  "1")==F) stop("mehrdeutige IDs im table!")
  if(length(setdiff(c("PROBE_ID"), names(nobkgd)) != 0)) stop(' "PROBE_ID" nicht gegeben - skript modifizieren')

  #bezugdatenframe zum umbennenen der IDs bauen
  spalten <- data.frame(old_colname=names(nobkgd))
  spalten$originalid_from_colname <- stringr::str_trim(sapply(stringr::str_split(as.character(spalten$old_colname), "\\.AVG_Signal|\\.Detection.Pval"), "[", 1))

  # spalten$id_from_colname <- as.character(stringr::str_replace(spalten$originalid_from_colname, "^X", ""))
  spalten$id_from_colname <- as.character(spalten$originalid_from_colname)

  #zur ueberpruefung modifiziere ich die zu importierenden ids laut uebersichtstabelle so, dass ich die importierten auf vollstaendigkeit ueberpruefen kann
  # samples_in_mod <- stringr::str_replace_all(samples_in, " ", ".")
  # samples_in_mod <- stringr::str_replace_all(samples_in_mod, "-", ".")
  samples_in_mod <- samples_in
  check = table(table(samples_in_mod))
  if(identical(names(check), "1")==F) stop(' Proben ID nicht unique in diesem Run')
  qlistcheck2 = venn2(samples_in_mod,spalten$id_from_colname, plotte = showVennplots)
  # str(qlistcheck2)
  if(length(setdiff(samples_in_mod, na.omit(spalten$id_from_colname))) != 0) stop("nicht alle zu importierende IDs konnten expressionsspaltenIDs zugeordnet werden")
  spalten$sample_id_toimport <- samples_in[match_hk(spalten$id_from_colname, samples_in_mod)]
  subset_all <-  sample_overview_l3[ sample_overview_l3$nobkgd_f == (filestable$nobkgd_f[i]) &
                                       sample_overview_l3$old_ID %in% spalten$sample_id_toimport,]
  head(subset_all)
  dim(subset_all)
  spalten$new_ID <- subset_all[match(spalten$sample_id_toimport, subset_all$old_ID),"new_ID"]
  zugeordnete_new_ID = spalten[is.na(spalten$sample_id_toimport)==F, "new_ID"]
  if(any(is.na(zugeordnete_new_ID))) stop("new_ID zuordnung zur ausgelesenen ID aus spalten$sample_id_toimport failed")
  spalten$new_colname <-  stringi::stri_replace_all_fixed(spalten$old_colname, spalten$originalid_from_colname, spalten$new_ID)
  ht(spalten,7)

  #zu importierende spalten
  importindex <- is.na(spalten$sample_id_toimport)==F
  spalten_good <- as.character(spalten[importindex, "old_colname"])
  new_names <- spalten[importindex,"new_colname"]
  nobkgd <- nobkgd[,c("PROBE_ID",  spalten_good), with = F]
  names(nobkgd) <- c("PROBE_ID",  new_names  )
  nobkgd[,new_names, with = F]
  test = nobkgd[1:2,]
  if(colseparator_used != ".") {
    for (j in new_names) set(nobkgd,j=j,value=as.numeric(gsub(pattern = colseparator_used, replacement = ".", nobkgd[[j]])))
    # showClassDF(nobkgd)
  }
  setkey(all_nobkgd, PROBE_ID)
  setkey(nobkgd, PROBE_ID)
  # message("object all_nobkgd is datatable check: ", is.data.table(all_nobkgd))
  all_nobkgd <- nobkgd[ all_nobkgd]
  if(length(grep("NA", names(all_nobkgd))) != 0) stop("NA im neuen spaltennamen aufgetreten - checke skript")
  #umbenennung dokumentieren
  spalten$file_used <-  filestable$nobkgd_f[i]
  id_overview <- rbind(id_overview, spalten)
}
ht(id_overview,5)

## ----checkallda----------------------------------------------------------
## unvollstaendige Files:
all_nobkgd = data.frame(all_nobkgd)
checkvollstaendig3 <- dim(all_nobkgd[,grep("Detection[ \\.]Pval", names(all_nobkgd), value=T)])
checkvollstaendig3
if(checkvollstaendig3[2] != dim(sample_overview_l3)[1]) {
  status3s01_allelinfosda = ("Problem found: Could not import expression information (Detection.Pval) from all samples\n")
  stop(status3s01_allelinfosda)
} else  {status3s01_allelinfosda = "Succesfully imported expression information (Detection.Pval) from all samples, no problem found.\n";
message(status3s01_allelinfosda)}

## ----checks-----------------------------------------------------------
all_detection = grep("Detection", names(all_nobkgd), value = T)
all_avg = grep("AVG", names(all_nobkgd), value = T)
all_detection = stringr::str_replace(all_detection, ".Detection.Pval", "")
all_avg = stringr::str_replace(all_avg, ".AVG_Signal", "")
qlist_import = venn2(all_avg, all_detection, plotte = showVennplots)
# str(qlist_import)

## ----checkallesignda-----------------------------------------------------
# unvollstaendige Files:
checkvollstaendig3b <- dim(all_nobkgd[ ,grep("AVG_Signal", names(all_nobkgd), value  = T)])
checkvollstaendig3b
if(checkvollstaendig3b[2] != dim(sample_overview_l3)[1]) {
  status3bs01_allelinfosda = ("Problem found: Could not import expression information (AVG_Signal) from all samples\n")
  stop(status3bs01_allelinfosda)
} else  {status3bs01_allelinfosda = "Succesfully imported expression information (AVG_Signal) from all samples, no problem found.\n";
message(status3bs01_allelinfosda)}

## ----detailcheck---------------------------------------------------------
all_detection <- grep("Detection", names(all_nobkgd), value = T)
all_avg <- grep("AVG", names(all_nobkgd), value = T)
all_detection <- stringr::str_replace(all_detection, ".Detection.Pval", "")
all_avg <- stringr::str_replace(all_avg, ".AVG_Signal", "")
qlist_import <- venn2(all_avg, all_detection, plotte = showVennplots)
# str(qlist_import)
stopifnot(length(qlist_import$q2) + length(qlist_import$q3) == 0)

## ----negcheck------------------------------------------------------------
hh(all_nobkgd,8)
dim(all_nobkgd)
min_nobk <- data.frame(colmin = sapply(all_nobkgd[, sapply(all_nobkgd, is.numeric)],
                                       function (x) min(x)), colnam = names(all_nobkgd[, sapply(all_nobkgd, is.numeric)]))
ht(min_nobk)
head(min_nobk[grep("AVG", min_nobk$colnam, ignore.case = T), "colmin"])
check_na <- table( min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA="ifany")
check_na
table(is.na(names(check_na)))
if(any(is.na(names(check_na))) != 0){
  status3s02_NAexprwerte = ("Problem found: 'NA' as expression value not allowed.\n")
  warning(status3s02_NAexprwerte)
} else {status3s02_NAexprwerte = "Did not find any 'NA' value within expression values, no problem identified.\n";
message(status3s02_NAexprwerte)}

## ----understand----------------------------------------------------------

if(grepl("^Problem found", status3s02_NAexprwerte))  {
  sm_checkneg <- min_nobk[grepl("AVG", min_nobk$colnam) & (min_nobk$colmin <0 | is.na(min_nobk$colmin)),]
  sm_checkneg
  sm_checkneg2 <- all_nobkgd[ , c('PROBE_ID', sm_checkneg$colnam)]
  showNA(sm_checkneg2)
  sm_checkNA = sm_checkneg2[ apply(sm_checkneg2, 1, function(x) any(is.na(x))),]
  sm_checkNA  # NA verstanden, das sind die probes, wo eines fehlt
  # rauswerfen der missing probe
  # loaded3 = load("obj/s01_1_fordoku.RData")
  # loaded3
  # rm(list = c(setdiff(loaded3, "missingilmns_det")))
  # missingilmns_det
  hh(all_nobkgd)
  na_ilmns <- all_nobkgd[,grep("AVG", names(all_nobkgd))]
  hh(na_ilmns)
  na_ilmns <- as.matrix(na_ilmns)
  rownames(na_ilmns) <- all_nobkgd$PROBE_ID
  na_ilmns <- t(na_ilmns)
  hh(na_ilmns)
  na_ilmns_countNAs <- showNA(na_ilmns)
  ht(na_ilmns_countNAs)
  mytable(na_ilmns_countNAs$NAs)
  ## CHANGE start3 20.8.18
  # ilmn2remove_NAfound <- rownames(na_ilmns_countNAs[na_ilmns_countNAs$NAs >
  # 0, , drop = F])
  ilmn2remove_NAfound <- na_ilmns_countNAs[na_ilmns_countNAs$NAs > 0, "var"]
  message("removing NA- ILMNs:\n", paste(ilmn2remove_NAfound, collapse = "\n"))
  ## CHANGE end3 20.8.18

  ilmn2remove_NAfound
  all_nobkgd_initial <- all_nobkgd
  all_nobkgd <- all_nobkgd[ all_nobkgd$PROBE_ID %nin% ilmn2remove_NAfound,]
  dim(all_nobkgd_initial)
  dim(all_nobkgd)
  s03_1_rausgeflogen_weilNA = ilmnAnnot014allgInfos[ paste0("ILMN_", ilmnAnnot014allgInfos$ilmn) %in% ilmn2remove_NAfound,]
} else {
  na_ilmns_countNAs <- data.frame(NAs = character(0))
  sm_checkNA <- character(0)
  s03_1_rausgeflogen_weilNA <- ilmnAnnot014allgInfos[0,]
  ilmn2remove_NAfound <- character(0)
}

## ----negatcheck2---------------------------------------------------------
#false 2213
checkneg <- table(min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA = "ifany")
checkneg
if(any(grepl("TRUE", names(checkneg)))){
  status3s02_negativeexprwerte = ("Problem found: negative gene expression values.")
  warning(status3s02_negativeexprwerte)
} else {status3s02_negativeexprwerte = "Did not find any negative expression value, no problem identified.\n";
message(status3s02_negativeexprwerte)}

## ----understand2---------------------------------------------------------
if( grepl("^Problem found", status3s02_negativeexprwerte)) {
  sm_checkneg <- min_nobk[grepl("AVG", min_nobk$colnam) & (min_nobk$colmin <0 | is.na(min_nobk$colmin)),]
  sm_checkneg
  sm_checkneg2 <- all_nobkgd[ , c('PROBE_ID', sm_checkneg$colnam)]
  showNA(sm_checkneg2)
  sm_checkNeg3 <- sm_checkneg2[ apply(sm_checkneg2[,sapply(sm_checkneg2, is.numeric),drop = F], 1, function(x) any(is.na(x) == F & x <0)),]
  data.table(sm_checkNeg3  )
  neg_ind <- unlist(sm_checkNeg3[, grep("AVG", names(sm_checkNeg3)), drop = F])
  neg_ind <- neg_ind[neg_ind<0]
  neg_ind <- unique(sapply(stringr::str_split(names(neg_ind), "\\."), "[", 1))
  neg_ind
  grep("AVG",names(sm_checkNeg3), value=T)
  strange_probe_expr <- unlist(all_nobkgd[ all_nobkgd$PROBE_ID %in% sm_checkNeg3$PROBE_ID,grep("AVG",names(all_nobkgd))])
  strange_probe_expr[strange_probe_expr<0]
  # str(strange_probe_expr)
  # hist(strange_probe_expr, breaks = 111)
  ## wie sieht das individuum insgesamt aus?
  # loaded3 <- load("obj/s02_1_fordoku.RData")
  # loaded3
  # rm(list = c(setdiff(loaded3, "pc123")))

  # farbe <- as.numeric(rownames(pc123) %in% neg_ind)+1
  # try(rgl::plot3d(pc123[1:3], col = farbe, size = 7)) # sess
  # plot(pc123[1:2], col = farbe, cex = farbe, pch = "x")
} else neg_ind = NULL
neg_ind

## ----indout--------------------------------------------------------------
all_minexpressval <- min(apply(all_nobkgd[, grep("Signal", names(all_nobkgd))], 2, function(x) min(x[x>=0])))
all_minexpressval
sample_overview_l4 <- sample_overview_l3
if(length(neg_ind) >0) negvalue_name2change = grep(paste(neg_ind, collapse = "|"),
                                                   names(all_nobkgd),
                                                   value = T) else negvalue_name2change = NULL
negvalue_name2change <- negvalue_name2change[grep("Signal", negvalue_name2change)]
dim(all_nobkgd)
if(length(neg_ind) >0) {
  negprobes <- all_nobkgd$PROBE_ID[unique(as.vector(unlist(apply(all_nobkgd[ ,names(all_nobkgd) %in% negvalue_name2change, drop = F],
                                                                 1, function(x) which(x<0)))))]
  all_nobkgd[ , names(all_nobkgd) %in% negvalue_name2change][all_nobkgd[ ,names(all_nobkgd) %in% negvalue_name2change]<0] = all_minexpressval
} else negprobes = NULL
dim(all_nobkgd)
dim(sample_overview_l4)

## ----whcheckneg----------------------------------------------------------
min_nobk <- data.frame(colmin = sapply(all_nobkgd[, sapply(all_nobkgd, is.numeric)],
                                       function(x) min(x)),colnam = names(all_nobkgd[, sapply(all_nobkgd, is.numeric)]))
ht(min_nobk)
head(min_nobk[grep("AVG", min_nobk$colnam), "colmin"])
checkneg <- table(min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA = "ifany") #false 2213
checkneg
table(is.na(names(checkneg)))
if(any(grepl("TRUE", names(checkneg))) | any(is.na(names(checkneg))) != 0){
  status3s02b_negativeexprwerte = ("Problem found: negative expression values or NA values still exists but are not allowed.\n")
  warning(status3s02b_negativeexprwerte)
} else  {status3s02b_negativeexprwerte = "Negative expression values or NA expression values do not exist anymore, no problem identified anymore.\n"; if(grepl("^Problem found", status3s02_negativeexprwerte)) message(status3s02b_negativeexprwerte)}

## ----checkNoBkgdFile-----------------------------------------------------
dim(filestable)
genenames <- setdiff(names(all_nobkgd), grep("PROBE_ID", names(all_nobkgd), value = T))
# str(genenames)
dim(all_nobkgd)
all_nobkgd <- all_nobkgd[, c("PROBE_ID", genenames)]
dim(all_nobkgd)

## ----save----------------------------------------------------------------


ht12object$rawdataWOcons_joined = all_nobkgd

# ablegen der namen, die im erstellen der expression verwendet werden inkl. der automatischen umgewandelten Namen von R
ht(id_overview[ is.na(id_overview$new_ID) ==F,],2)
ht12object$id_overview = id_overview


## ----lumi----------------------------------------------------------------
#als expressionset-txt file fuer lumi-import, d.h. entfernen einiger nicht-benoetigter spalten
time6 <- Sys.time()
use_names <- names(all_nobkgd)
use_names <- c("PROBE_ID", grep("AVG_Signal",use_names, value=T),grep("Detection.Pval",use_names, value=T) )
stopifnot(all(use_names %in% names(all_nobkgd)))
all_nobkgd_eset_input <- all_nobkgd[,names(all_nobkgd) %in% use_names]
dim(all_nobkgd_eset_input)
table(names(all_nobkgd_eset_input) %in% use_names , useNA="always")
setdiff( names(all_nobkgd_eset_input) ,use_names)
# filename_lumibauen <- tempfile()
# filename_lumibauen
# debugnames <- names(all_nobkgd_eset_input)[c(1:20, (dim(all_nobkgd_eset_input)[2]-4):dim(all_nobkgd_eset_input)[2])]
# debugnames

# write.delim(all_nobkgd_eset_input[, debugnames], filename_lumibauen)
# write.delim(all_nobkgd_eset_input, filename_lumibauen)
time6 <-  Sys.time() -time6
time7 <- Sys.time()
names(all_nobkgd_eset_input)[1:22]

# all_detection = grep("Detection", names(all_nobkgd), value = T)
# all_avg = grep("AVG", names(all_nobkgd), value = T)
all_nobkgd_eset_input[1,]
all_detection <- grep("Detection", names(all_nobkgd_eset_input), value = T)
all_detection
all_avg <- grep("AVG", names(all_nobkgd_eset_input), value = T)
all_avg
all_detection <- stringr::str_replace(all_detection, ".Detection.Pval", "")
all_avg <- stringr::str_replace(all_avg, ".AVG_Signal", "")
qlist_lumicheck = venn2(all_avg, all_detection, plotte = showVennplots )
# str(qlist_lumicheck)
neg_ind
stopifnot(length(qlist_lumicheck$q2) + length(qlist_lumicheck$q3) == 0)
# message("importing data into Expressionset via lumi...")
# all_nobkgd_eset <- lumi::lumiR(filename_lumibauen,
                               # lib.mapping="lumiHumanIDMapping",
                               # columnNameGrepPattern = list(exprs='AVG_Signal|AVG_SIGNAL',
                                                            # detection='Detection|DETECTION'))
# class(all_nobkgd_eset)       # ExpressionSet
# if(class(all_nobkgd_eset) == "ExpressionSet") all_nobkgd_eset
# hh(detection(all_nobkgd_eset))
# hh(exprs(all_nobkgd_eset))
# stopifnot(all(apply(exprs(all_nobkgd_eset), 2, function(x) any(is.na(x))) == F))
# table(apply(detection(all_nobkgd_eset), 2, function(x) any(is.na(x))))
# stopifnot(all(apply(detection(all_nobkgd_eset), 2, function(x) any(is.na(x))) == F))

## ----savi----------------------------------------------------------------
# status fuer docku
# stati <- grep("^status",ls(), value = T)
# stati
# for(i in stati) print(get(i))

# nacheck_fin <- showNA(exprs(all_nobkgd_eset))
# mytable(nacheck_fin$NAs)
# all_nobkgd_eset
time7 <- Sys.time() - time7


# samples as used
sample_overview_l4doku <- sample_overview_l3_initial
mytable(sample_overview_l4doku$reason4exclusion)
mytable(sample_overview_l4doku$in_study)


sample_overview_l4 <- sample_overview_l4doku

######################################################################3
## create a file incl controls
message("\nImporting control expression data....")
sample_overview_l4_initial = sample_overview_l4
sample_overview_l4 <- sample_overview_l4_initial[sample_overview_l4_initial$in_study,]
dim(sample_overview_l4_initial)
dim(sample_overview_l4)

# checken unique proben ids
table(table(sample_overview_l4$new_ID))  #nein
if(length(table(table(sample_overview_l4$new_ID))) != 1) stop("IDs (column new_ID) must be unique....stopping...")

# Check auf vollstaendigkeit der kontrollfiles
check <- dim(sample_overview_l4)
sample_overview_l4 <- sample_overview_l4[ is.na(sample_overview_l4$con_f) == F, ]
if(identical(check, dim(sample_overview_l4)) == F) stop("not all  control filenames given ...")

# annotation der proben
head(annotcon)
table(table(annotcon$Probe_Id))
table(table(annotcon$Array_Address_Id))
#   table(annotcon$Reporter_Group_id)
sum(table(annotcon$Reporter_Group_id))
# str(unique(annotcon$Array_Address_Id))

# controlprobes_fn = paste0(basicpath, "/genstat/01_daten/1104_ht12_leheart_ok/libraries/illumina_all_control_Array_Address_Id_110416hk.txt", sep="")
# all_con <- unique(read.delim(controlprobes_fn, stringsAsFactors=F))     #liste aller vorhandenen ilmns einlesen
# str(all_con)   #883 obs. of  1 variable:
# table(table(all_con))
#
# qlist0 = venn2(all_con$ProbeID, annotcon$Array_Address_Id)
# str(qlist0)

## ----importcontr---------------------------------------------------------
# bauen des files
# bauen der fileuebersicht zum laden
names(sample_overview_l4)
filestable <- unique(sample_overview_l4[,c("fileset_id", "nobkgd_f", "con_f", "sample_f")])


#con importieren
id_overview <- c()
all_con <- data.frame(ProbeID = unique(annotcon$Array_Address_Id))     #liste aller vorhandenen ilmns einlesen
# str(all_con)   #883 obs. of  1 variable:
if(dim(all_con)[1] != 883) stop("nicht alle 883 kontroll-IDs geladen!")

for(i in 1:length(filestable[,1])){
  #alles laden
  message("Import control file ",i,  ": " ,filestable$con_f[i], "...\n")
  con <- suppressWarnings(data.table::fread(filestable$con_f[i]))
  dim(con)
  if(sum((con$TargetID == "LOW_STRINGENCY_HYB" & con$ProbeID %in% c(1110170, 4610291, 2510500, 4010327))==T) != 4) stop("checke, herausgweorfene redundante - ")
  if(names(table(table(con$ProbeID %in% c(1110170, 4610291, 2510500, 4010327)))[1])!="8")  stop("checke, herausgweorfene redundante - 4 muessen doppelt sein")
  if(identical(sort(unique(all_con[,1])), sort(unique(con$ProbeID)))==F) message("Did not find all 883 controls in file ", filestable$con_f[i])
  # print(setdiff(all_con[,1], con$ProbeID))
  # print(setdiff(con$ProbeID, all_con[,1]))
  con <- con[((con$TargetID == "LOW_STRINGENCY_HYB" & con$ProbeID %in% c(1110170, 4610291, 2510500, 4010327))==F),]
  dim(con)
  names(con)
  #expressiosn spalten reduzieren auf ids im table
  samples_in <- as.character(sample_overview_l4[as.character(sample_overview_l4$con_f) == as.character(filestable$con_f[i]), "old_ID"]) #die IDs, von denen ich daten importieren will
  if(length(grep("^X", samples_in)) != 0) stop("irgendwelche sampleIDs fangen  mit 'X' an -> skript modifizieren, weil z.Z. aus bequemlichkeit davon ausgehe, dass nach dataframe import ein beginnendes 'X' vom Import einer spalte herruehrt, die original eine zahl war!!!")
  if(length(samples_in) == 0) stop("keine IDs zu importieren - checke skript!")
  check <- table(table(unique(samples_in)))
  check
  if(identical(names(check), "1") == F) stop("given IDs in table not unique!")
  #if((length(grep("\\.", samples_in)) != 0 )) stop(paste("punkt im skript -> skript modifizieren!!! checke id:", grep("\\.", samples_in, value=T)))
  if(length(setdiff(c("TargetID", "ProbeID"), names(con)) != 0)) stop('mindestens ein sample_overview_l4ut von "TargetID", "ProbeID" fehlt - skript modifizieren')
  #bezugdatenframe zum umbennenen der IDs bauen
  spalten <- data.frame(old_colname = names(con))
  spalten$originalid_from_colname <- stringr::str_trim(sapply(stringr::str_split(as.character(spalten$old_colname),
                                                      "\\.AVG_Signal|\\.Detection.Pval|\\.BEAD_STDERR|\\.Avg_NBEADS"), "[", 1))
  # spalten$id_from_colname <- as.character(str_replace(spalten$originalid_from_colname, "^X", "")) #aus historischen gruenden drinne, stoert nicht
  spalten$id_from_colname <- as.character(spalten$originalid_from_colname)

  samples_in_mod <- samples_in
  check <- table(table(samples_in_mod))
  if(identical(names(check), "1") == F)
    stop('Proben ID nicht unique in diesem Run')
  qlistcheck2 = venn2(samples_in_mod,spalten$id_from_colname, plotte = showVennplots)
  # str(qlistcheck2)
  if(length(setdiff(samples_in_mod, na.omit(spalten$id_from_colname))) != 0)
    stop("nicht alle zu importierende IDs konnten expressionsspaltenIDs zugeordnet werden")
  spalten$sample_id_toimport <- samples_in[match_hk(spalten$id_from_colname, samples_in_mod)]

  # detailed checke to here
  subset_all <-  sample_overview_l4[(sample_overview_l4$con_f) ==
                                      (filestable$con_f[i]) & (sample_overview_l4$old_ID) %in% (spalten$sample_id_toimport),]
  head(subset_all)
  dim(subset_all)
  spalten$new_ID <- subset_all[match(spalten$sample_id_toimport, subset_all$old_ID),"new_ID"]
  if(any(is.na(spalten[is.na(spalten$sample_id_toimport)==F, "new_ID"])))
    stop("new_ID zuordnung aus spalten$sample_id_toimport failed")
  spalten$new_colname <- stringi::stri_replace_all_fixed(spalten$old_colname,
                                                         spalten$originalid_from_colname,
                                                         spalten$new_ID)
  ht(spalten)

  #zu importierende spalten
  spalten_good <- as.character(spalten[is.na(spalten$sample_id_toimport) == F, "old_colname"])
  new_names <- spalten[is.na(spalten$sample_id_toimport) == F, "new_colname"]
  con <- con[ , c("TargetID", "ProbeID", spalten_good), with = F]
  names(con) <- c("TargetID", "ProbeID", new_names)

  ## start CHANGED new chunk --> 20.8.18 check already here for duplicated columns --> do not merge them, but check for inconsistencies
  namescheck = venn2(names(all_con), names(con), plotte = showVennplots)
  ## check overlapping names
  overlapping_name= setdiff(namescheck$q1, "ProbeID")
  stopifnot(length(overlapping_name)<=1)
  if(length(overlapping_name)==1) {
    vgl_frame = data.table(all_con=all_con[,overlapping_name], con = con[match_hk(all_con$ProbeID, con$ProbeID), get(overlapping_name)])
    ws_check = vgl_frame[all_con != con]
    if(nrow(ws_check)>0) {
      message("I found contradictions in column ", overlapping_name, " for index i = ", i, " file ", filestable[, 1][i], " with previously merged data")
      print(ws_check)
      stop("Stopping. Please resolve!")
    } else con[,(overlapping_name):=NULL]
  }
  ## end  CHANGED new chunk --> 20.8.18
  all_con <- merge(all_con, con, by="ProbeID", all.x = T, sort = F)
  if(length(grep("NA", names(all_con))) != 0)
    stop("NA im neuen spaltennamen aufgetreten - checke skript")
  #umbenennung dokumentieren

  spalten$file_used <-  filestable$con_f[i]
  id_overview <- data.frame(rbindlist(list(id_overview, spalten)))
}
ht(id_overview,2)
ht(all_con)
names(all_con) = stringr::str_replace_all(names(all_con), " ", ".")

## ----exclusion-----------------------------------------------------------
escludecols <- grep( "BEAD_STDERR$|Avg_NBEADS$", names(all_con), value = T)
# str(escludecols)
all_con <- all_con[ , names(all_con) %nin% escludecols]

## ----colseparat----------------------------------------------------------

# where does the colseparator come from?
colseparator

# what is the current colseparator?
colseparator_used

# change to "." if "." is not set in colseparator_used
if(colseparator_used != ".") {
  class(all_con)
  all_con <- data.table(all_con)
  numcol = grep( "AVG_Signal$|Detection.Pval$", names(all_con), value = T)
  for (j in numcol) set(all_con, j = j, value = as.numeric(gsub(pattern = colseparator_used, replacement = ".", all_con[[j]])))
  # table(showClassDF(all_con))
  all_con <- data.frame(all_con)
}

## ----checkalledacons-----------------------------------------------------
## unvollstaendige Files:
checkConvollstaendig3 <- dim(all_con[,grep("Detection.Pval", names(all_con), value = T)])
checkConvollstaendig3
if(checkConvollstaendig3[2] != dim(sample_overview_l4)[1]) {
  status4s01_allelinfosda = ("Problem found: Could not import control-expression information  (Detection.Pval) from all samples\n")
  stop(status4s01_allelinfosda)
} else {status4s01_allelinfosda = "Succesfully imported control-expression information from  (Detection.Pval) all samples, no problem found.\n"; message(status4s01_allelinfosda)}

## ----checkConallesigndacons-------------------------------------------------
# unvollstaendige Files:
checkConvollstaendig3b <- dim(all_con[,grep("AVG_Signal", names(all_con), value = T)])
checkConvollstaendig3b
if(checkConvollstaendig3b[2] != dim(sample_overview_l4)[1]) {
  status3bs01_allelinfosda = ("Problem found: Could not import expression information (AVG_Signal) from all samples for control probes\n")
  stop(status3bs01_allelinfosda)
} else  {status3bs01_allelinfosda = "Succesfully imported expression information (AVG_Signal) from all samples for control probes, no problem found.\n"; message(status3bs01_allelinfosda)}

## ----detailcheckConcons-----------------------------------------------------
all_detection <- grep("Detection", names(all_con), value = T)
all_avg <- grep("AVG", names(all_con), value = T)
all_detection <- stringr::str_replace(all_detection, ".Detection.Pval", "")
all_avg <- stringr::str_replace(all_avg, ".AVG_Signal", "")
qlist_import <- venn2(all_avg, all_detection, plotte = showVennplots)
# str(qlist_import)
stopifnot(length(qlist_import$q2) + length(qlist_import$q3) == 0)

## ----negatcheckCon----------------------------------------------------------
hh(all_con, 8)
dim(all_con)
min_nobk <- data.frame(colmin = sapply(all_con[, sapply(all_con, is.numeric)],
                                       function(x) min(x)),colnam = names(all_con[, sapply(all_con, is.numeric)]))
ht(min_nobk)
head(min_nobk[grep("AVG", min_nobk$colnam), "colmin"])
checkCon_na <- table( min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA = "ifany")      #false 2213
checkCon_na
table(is.na(names(checkCon_na)))
if (any(is.na(names(checkCon_na))) != 0){
  status4s02_NAexprwerte = ("Problem found: 'NA' as control-expression value not allowed, those control probe will be removed and checked again for NAs.\n")
  warning (status4s02_NAexprwerte)
} else {status4s02_NAexprwerte = "Did not find any 'NA' value within control-expression values, no problem identified.\n";
message(status4s02_NAexprwerte)}

## ----understand----------------------------------------------------------
if (grepl("^Problem found",status4s02_NAexprwerte)) {
  sm_checkConneg <- min_nobk[grepl("AVG", min_nobk$colnam) & (min_nobk$colmin <0 | is.na(min_nobk$colmin)), ]
  sm_checkConneg
  head(sm_checkConneg, 2)
  hh(all_con, 9)
  sm_checkConneg2 <- all_con[, c('ProbeID', 'TargetID', sm_checkConneg$colnam)] # CHANGED 20.8.18  TargetID instead of TargetID.x
  showNA(sm_checkConneg2)
  sm_checkConNA <- sm_checkConneg2[apply(sm_checkConneg2, 1, function(x) any(is.na(x))),,drop=F]
  sm_checkConNA  # NA verstanden, das sind die probes, wo eines fehlt

  # rauswerfen der missing probe
  # loaded3 = load("obj/s01_1_fordoku.RData")
  # loaded3
  # rm(list = c(setdiff(loaded3, "missing_ilmnCons")))
  # missing_ilmnCons
  # ilmn2remove = paste0("ILMN_",missing_ilmnCons$missing_ilmnCons) # TODO falls mal was gefunden wird, hier zeile checken
  hh(all_con)
  na_ilmns <- all_con[,grep("AVG", names(all_con))]
  hh(na_ilmns)
  na_ilmns <- as.matrix(na_ilmns)
  rownames(na_ilmns) <- all_con$ProbeID
  na_ilmns <- t(na_ilmns)
  hh(na_ilmns)
  na_ilmns_countNAs_con <- showNA(na_ilmns)
  ht(na_ilmns_countNAs_con)
  mytable(na_ilmns_countNAs_con$NAs)
  ilmn2remove_NAfound_con <- na_ilmns_countNAs_con[na_ilmns_countNAs_con$NAs>0,"var" ]
  ilmn2remove_NAfound_con
  message("removing NA- ILMNs from controls:\n", paste(ilmn2remove_NAfound_con, collapse = "\n"))

  s04_1_rausgeflogen_weilNA <- annotcon[ annotcon$Array_Address_Id %in%ilmn2remove_NAfound_con,]
  print(s04_1_rausgeflogen_weilNA)
  all_con_initial <- all_con
  all_con <- all_con[all_con$ProbeID %nin% ilmn2remove_NAfound_con, ]
  dim(all_con_initial)
  dim(all_con)
} else {
  na_ilmns_countNAs_con <- data.frame(NAs = character(0))
  s04_1_rausgeflogen_weilNA <- annotcon[0, ]
  sm_checkConNA <- character(0)
  ilmn2remove_NAfound_con <- character(0)
}

## ----negatcheckCon2---------------------------------------------------------
min_nobk <- data.frame(colmin = sapply(all_con[, sapply(all_con, is.numeric)],
                                       function(x) min(x)),colnam = names(all_con[, sapply(all_con, is.numeric)]))
checkConneg <- table( min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA = "ifany") #false 2213
checkConneg
if (any(grepl("TRUE", names(checkConneg)) & is.na(names(checkConneg)) == F)){
  status4s02_negativeexprwerte <- ("Problem found: negative control-expression values not allowed.")
  warning (status4s02_negativeexprwerte)
} else {status4s02_negativeexprwerte <- "Did not find any negative control-expression values, no problem identified.\n";
message(status4s02_negativeexprwerte)}

## ----understand2---------------------------------------------------------

if (grepl("^Problem found", status4s02_negativeexprwerte)) {
  sm_checkConneg <- min_nobk[grepl("AVG", min_nobk$colnam) & (min_nobk$colmin <0 | is.na(min_nobk$colmin)), ]
  sm_checkConneg
  sm_checkConneg2 <- all_con[, c('ProbeID', 'TargetID', sm_checkConneg$colnam)]  # CHANGED 20.8.18  TargetID instead of TargetID.x
  showNA(sm_checkConneg2)
  sm_checkConNeg3 <- sm_checkConneg2[apply(sm_checkConneg2[, sapply(sm_checkConneg2, is.numeric), drop = F],
                                     1, function(x) any(is.na(x) == F & x <0)), ]
  data.table(sm_checkConNeg3)
  neg_indCon <- unlist(sm_checkConNeg3[, grep("AVG", names(sm_checkConNeg3)), drop = F])
  neg_indCon <- neg_indCon[neg_indCon<0]
  neg_indCon <- unique(sapply(stringr::str_split(names(neg_indCon), "\\."), "[", 1))
  neg_indCon
  grep("AVG",names(sm_checkConNeg3), value = T)
  strange_probe_expr <- unlist(all_con[all_con$ProbeID %in% sm_checkConNeg3$ProbeID, grep("AVG", names(all_con))])
  strange_probe_expr[strange_probe_expr<0]
  # str(strange_probe_expr)
  # hist(strange_probe_expr, breaks = 111)

  # wie sieht das individuum insgesamt aus?
  # loaded3 <- load("obj/s01_1_fordoku.RData")
  # loaded3
  # rm(list = c(setdiff(loaded3, "pc123")))
  #
  # farbe <- as.numeric(rownames(pc123) %in% neg_indCon) + 1
  # try(rgl::plot3d(pc123[1:3], col = farbe, size = 7))
  # plot(pc123[1:2], col = farbe, cex = farbe, pch = "x")
} else neg_indCon = NULL
neg_indCon

## ----indraus-------------------------------------------------------------
all_minexpressval <- min(apply(all_con[, grep("Signal", names(all_con))], 2, function(x) min(x[x >= 0])))
all_minexpressval
sample_overview_l5 <- sample_overview_l4
if (length(neg_indCon) >0) negvalue_name2change_con <- grep(paste(neg_indCon, collapse = "|"), names(all_con), value = T) else negvalue_name2change_con = NULL
negvalue_name2change_con <- negvalue_name2change_con[grep("Signal", negvalue_name2change_con)]
if (length(neg_indCon) >0) {
  negprobes_con <- all_con$PROBE_ID[unique(as.vector(unlist(apply(all_con[, names(all_con) %in% negvalue_name2change_con, drop = F], 1,
                                                              function (x) which(x<0)))))] ## TODO unbedingt TESTEN
  all_con[, names(all_con) %in% negvalue_name2change_con][all_con[, names(all_con) %in% negvalue_name2change_con]<0] = all_minexpressval
} else negprobes_con = NULL
dim(all_con)
dim(sample_overview_l5)
dim(sample_overview_l4)

## ----whcheckConneg----------------------------------------------------------
min_nobk <- data.frame(colmin = sapply(all_con[, sapply(all_con, is.numeric)],
                                       function(x) min(x)), colnam = names(all_con[, sapply(all_con, is.numeric)]))
ht(min_nobk)
head(min_nobk[grep("AVG", min_nobk$colnam), "colmin"])
checkConneg <- table(min_nobk[grep("AVG", min_nobk$colnam), "colmin"] <0, useNA = "ifany") # false 2213
checkConneg
table(is.na(names(checkConneg)))
if(any(grepl("TRUE", names(checkConneg))) | any(is.na(names(checkConneg))) !=0){
  status4s02b_negativeexprwerte <- ("Problem found: negative control-expression values or control-NA values still exists but are not allowed.\n")
  warning (status4s02b_negativeexprwerte)
} else {status4s02b_negativeexprwerte <- "Negative control-expression values or NA control-expression values do not exist anymore, no problem identified anymore.\n";
if (grepl("^Problem found", status4s02_negativeexprwerte)) message(status4s02b_negativeexprwerte)
}
## CHANGED start 20.8.18 --> redundant columns are now removed earlier. no need for this chunk anymore
## ----redundati-----------------------------------------------------------
# if (dim(filestable)[1]>1) {
#   for( j in c("TargetID")){
#     redundant <- all_con[, grep(j, names(all_con))]
#
#     message("check whether similar fields in input data have indeed the same data (if '1', than test is ok, i.e. the same information was found accross those fields):")
#     uniquecheck <- table(length(unlist(apply(redundant, 1, function(x) unique(na.omit(x))))))
#     message(uniquecheck)   #883
#     if (uniquecheck != 1) {
# 	    message(names(redundant))
#     message("Those entries are (NA-filtered):")
#     message(table(unlist(lapply(redundant, function(x) length(na.omit(x))))))
#       stop (paste("similar fields in input data have NOT the same data, found in ", j))
#     }
#     # NA (nicht alle Werte von allen Kontrollproben gefunden)  betrifft folgende IDs
#     table((sapply(all_con, function(x) any(is.na(x) == T))))
#     dim(all_con)
#     all_con[1:11,1:5]
#     assign(j, unlist(apply(redundant, 1, function(x) if(sum(is.na(x) == T) == length(x)) x = NA else unique(na.omit(x)))))
#   }
#   all_con$TargetID <- TargetID
# }
## CHANGED end 20.8.18 --> redundant columns are now removed earlier. no need for this chung anymore

# file ohne redundante annotationbauen
genenames <- setdiff(names(all_con), grep("TargetID|ProbeID", names(all_con), value = T))
checkCon <- length(genenames)
if (length(genenames) != dim(sample_overview_l5)[1]*2)
  stop ("nicht alle gennamen wie erwartet")
all_con <- all_con[, c("ProbeID","TargetID", genenames)]
checkCon <- dim(all_con)
if (checkCon[2]-2 != dim(sample_overview_l5)[1]*2)
  stop ("reduktion von all_con nicht geklappt")

## ----reportMissings------------------------------------------------------
notCompliteProbesCon <- all_con[, "ProbeID"][(apply(as.matrix(all_con[, sapply(all_con, is.numeric)]),
                                                 1, function(x) any(is.na(x) ==T)))]
notCompliteProbesCon
table(all_con[all_con$ProbeID %in% notCompliteProbesCon, 2], useNA = "always")
table(annotcon[annotcon$Array_Address_Id %in% notCompliteProbesCon, 3])

##speichern

ht12object$rawdataOnlycons_joined = all_con

## ----saveercc------------------------------------------------------------
ercc <- annotcon[grep("ERCC", annotcon$Reporter_Group_Name), ]
ht(ercc, 2)
ercc$nuid <- lumi::IlluminaID2nuID(IlluminaID = ercc$Probe_Id)[, "nuID"]
ht(ercc, 2)
##speichern der "ERCC" proben
ht12object$ercc = ercc

## ----creationOneFile-----------------------------------------------------
all_nobkgd[1:11, 1:5]
all_con[1:11, 1:5]
dim(all_nobkgd)

##mergen
#fehlende spalten in all_con hinzufuegen und ueberfzaehlige loeschen
#kontrolliern, dass info von BEADS wirklich nicht mehr da ist
all_nobkgd2 <- all_nobkgd[, !names(all_nobkgd) %in% grep("Avg_NBEADS$|BEAD_STDERR", names(all_nobkgd), value = T)]
all_con$PROBE_ID <- annotcon[match(all_con$ProbeID, annotcon$Array_Address_Id), "Probe_Id"]
if(any(is.na(all_con$PROBE_ID) == T))
  stop("missings wo keine erwartet code 12444")
setdiff(names(all_nobkgd2), names(all_con))

#welche spalten haben NA eintraege?
na_entries <- all_nobkgd2[, 1][apply(as.matrix(all_nobkgd2[, sapply(all_nobkgd2, is.numeric)]),
                                     1, function(x) any(is.na(x) == T))]
na_entries
all_con$ProbeID <- NULL # das ist nicht PROBE_ID :)
all_con$TargetID <- NULL
qlistcheckCon <- venn2(names(all_nobkgd2), names(all_con), plotte = showVennplots)
# str(qlistcheckCon)
stopifnot(length(qlistcheckCon$q2) + length(qlistcheckCon$q3) == 0)

#sind alle eindeutig
qlist3 <- venn3(ercc$Probe_Id, all_con$PROBE_ID, all_nobkgd2$PROBE_ID, plotte = showVennplots)
# str(qlist3)
s33 <- annotcon[ annotcon$Probe_Id %in% qlist3$q3, ]
mytable(s33$Reporter_Group_Name)
s33 <- annotcon[ annotcon$Probe_Id %in% qlist3$q1, ]
mytable(stringr::str_sub(s33$Reporter_Group_Name, 1, 4))
columns_both <-qlistcheckCon$q1
qlist5 <- venn2(all_con$PROBE_ID, all_nobkgd2$PROBE_ID, plotte = showVennplots)
# str(qlist5)
ilmnboth <- qlist5$q1
all_con_overlap <- all_con[all_con$PROBE_ID %in% ilmnboth,columns_both]
all_con_overlap <- all_con_overlap[order(all_con_overlap$PROBE_ID),sort(names(all_con_overlap))]
all_nobkgd_overlap <- all_nobkgd2[all_nobkgd2$PROBE_ID %in% ilmnboth,columns_both]
all_nobkgd_overlap <- all_nobkgd_overlap[order(all_nobkgd_overlap$PROBE_ID),sort(names(all_nobkgd_overlap))]
all_con_overlap <- as.matrix(all_con_overlap)
dim(all_con_overlap)
all_nobkgd_overlap <- as.matrix(all_nobkgd_overlap)
dim(all_nobkgd_overlap)
hh(all_con_overlap)
hh(all_nobkgd_overlap)

as.vector(unlist(all_con_overlap) )[suppressWarnings(is.na(as.numeric(as.vector(unlist(all_con_overlap)) )))]
ercc_check1 <- data.frame(all_con = suppressWarnings(round(as.numeric(as.vector(unlist(all_con_overlap)) ),1)),
                          all_nobkgd =  suppressWarnings(round(as.numeric(as.vector(all_nobkgd_overlap)),1))
)
ht(ercc_check1)
# plot(ercc_check1$all_con, ercc_check1$all_nobkgd)
ercc_check2 <- sum(ercc_check1$all_con != ercc_check1$all_nobkgd)
sm <- ercc_check1[ercc_check1$all_con != ercc_check1$all_nobkgd,]
sm
ungleichkriterium <- 1.01*ercc_check1$all_con < ercc_check1$all_nobkgd & 0.99*ercc_check1$all_con > ercc_check1$all_nobkgd
sm2 = ercc_check1[ungleichkriterium, ]
sm2[is.na(sm2$all_con) == F, ]
ercc_check <- sum(ungleichkriterium, na.rm = T) == 0
if (ercc_check == F) {
  status4s05_ercc_check <- paste0("Problem found: ERCC expression levels differ  for identical probes in expression files and control files")
  stop (status4s05_ercc_check)
} else {status4s05_ercc_check = "ERCC expression levels are the same for identical probes in expression files and control files, no problem identified\n"; message(status4s05_ercc_check)
}
ilmncon <- all_con$PROBE_ID[!all_con$PROBE_ID %in% all_nobkgd2$PROBE_ID]
# str(ilmncon) #784 eintraege ok

# if(length(ilmncon) != 784) stop("haette hier eigentlich eine laenge von 784 erwartet...")
table(table(ilmncon))
all_con_new <- all_con[all_con$PROBE_ID %in% ilmncon, ]
dim(all_con_new)

##anfuegen von NA spalten fuer (eventuell noch) fehlnde control file, zB war das mal. "ControlProbeProfile_HSS207-HSS232.txt"
names_toadd <- setdiff(names(all_nobkgd2), names(all_con_new))
if (length(names_toadd) >0) {
  for (i in names_toadd) {
    # print(i)
    all_con_new$toadd <- NA
    names(all_con_new)[ names(all_con_new)=="toadd"] <- i
  }
} #TODO Reporten
names_toadd
qlist_prerbing <- venn2(names(all_nobkgd2), names(all_con_new), plotte = showVennplots)
# str(qlist_prerbing)
stopifnot(length(qlist_prerbing$q2) + length(qlist_prerbing$q3) == 0)
total_nobkgd <- rbind(all_nobkgd2, all_con_new)
dim(all_nobkgd2)
dim(total_nobkgd)
table(showNA2(total_nobkgd)$N)

#welche spalten haben NA eintraege?
ilmns_after_merge_with_NA_entries <- total_nobkgd[,1][apply(as.matrix(total_nobkgd[, sapply(total_nobkgd, is.numeric)]),
                                                            1, function(x) any(is.na(x) ==T))]
ilmns_after_merge_with_NA_entries
cols_after_merge_with_NA_entries <- names(total_nobkgd)[apply(as.matrix(total_nobkgd[, sapply(total_nobkgd, is.numeric)]),
                                                              2, function(x) any(is.na(x) ==T))]
cols_after_merge_with_NA_entries
table(annotcon[annotcon$Probe_Id %in% na_entries, 3])

#bauen eines Expression-Set Objektes

time5 <- Sys.time()
total_nobkgd_input_fn <- tempfile()
total_nobkgd_input_fn
message("Writing temporary file named",total_nobkgd_input_fn," in order to create expressionset of gene probes as well as control probes via R-package lumi... ")
data.table::fwrite(total_nobkgd, total_nobkgd_input_fn, sep = "\t",
                   row.names = F, col.names = T, quote = F)
vgl1 = venn2(total_nobkgd$PROBE_ID, all_nobkgd$PROBE_ID, plotte = showVennplots)
vgl2 = venn2(names(total_nobkgd), names(all_nobkgd), plotte = showVennplots)
# str(vgl2)
time5 <- Sys.time() - time5
time6 <- Sys.time()
message("Creating  expressionset of gene probes as well as control probes via lumi... ")
total_nobkgd_eset <- suppressWarnings(lumi::lumiR(total_nobkgd_input_fn, lib.mapping = "lumiHumanIDMapping", verbose = T,
                           columnNameGrepPattern = list(exprs='AVG_Signal|AVG_SIGNAL',
                                                        detection='Detection|DETECTION')))
class(total_nobkgd_eset) #"ExpressionSet"
total_nobkgd_eset

hh(detection(total_nobkgd_eset))
hh(exprs(total_nobkgd_eset))
head(names(total_nobkgd))

conprobes = vgl1$q2
## eset auch noch parallel einschraenken auf noncontrollen

all_nobkgd_eset = total_nobkgd_eset[featureNames(total_nobkgd_eset) %nin% conprobes,]
stopifnot(all(apply(exprs(all_nobkgd_eset), 2, function(x) any(is.na(x))) == F))
table(apply(detection(all_nobkgd_eset), 2, function(x) any(is.na(x))))
stopifnot(all(apply(detection(all_nobkgd_eset), 2, function(x) any(is.na(x))) == F))

all_nobkgd_eset_dim <- dim(all_nobkgd_eset)
all_nobkgd_eset_dim




## ----save2---------------------------------------------------------------
stopifnot(all(apply(exprs(total_nobkgd_eset), 2, function(x) any(is.na(x))) == F))
table(apply(detection(total_nobkgd_eset), 2, function(x) any(is.na(x))))
stopifnot(all(apply(detection(total_nobkgd_eset), 2, function(x) any(is.na(x))) == F))


# status fuer docku
stati = grep("^status", ls(), value = T)
stati
# for(i in stati) print(get(i))

# samples as used
sample_overview_l5doku <- sample_overview_l4_initial
mytable(sample_overview_l5doku$reason4exclusion)
mytable(sample_overview_l5doku$in_study)

sample_overview_l5 <- sample_overview_l5doku

total_nobkgd_eset_dim = dim(total_nobkgd_eset)

fordokuCon <- c(stati, "notCompliteProbesCon", "s04_1_rausgeflogen_weilNA",
                           "ilmns_after_merge_with_NA_entries",
                           "cols_after_merge_with_NA_entries",
                           "total_nobkgd_eset_dim", "neg_indCon")


## ----save, results='markup', echo=T, eval=T------------------------------
time6
time5


ht12object$chipsamples = sample_overview_l5
tosave = c(           "neg_ind",
           "sm_checkNA",
           "negvalue_name2change",
           "negvalue_name2change_con",
           "negprobes",
           "negprobes_con",
           "ilmn2remove_NAfound",
           'ilmn2remove_NAfound_con',
           "na_ilmns_countNAs",
           "na_ilmns_countNAs_con",
           "all_nobkgd_eset_dim",
           'sample_overview_l3', "checkConneg", fordokuCon)

ht12object$all_nobkgd_eset = all_nobkgd_eset
ht12object$total_nobkgd_eset = total_nobkgd_eset

stopifnot(sum(duplicated(tosave))==0)
ht12object$dokuobjects_createExpressionSet = lapply(tosave, function(x) get(x))


names(ht12object$dokuobjects_createExpressionSet) = tosave
# str(ht12object$dokuobjects_createExpressionSet)
ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object
  }
