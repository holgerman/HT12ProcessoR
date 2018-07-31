#' Filter samples with atypical number of expressed genes and annotate probes for their expression level
#'
#' @description This function filters from an HT12prepro object samples with atypical number of expressed genes and annotates probes for their respective expression level

#' @param ht12object A list object of class HT12prepro created with function createExpressionSet()
#' @param paramfile Path to the file specifying parameters.
#' @param filter1ind_expressedGenes Filter for extreme number of 'Detected.Genes..0.01.' at detection p-value 0.01. Valid is: 'Detected.Genes..0.01.' < Median - [value]x IQR AND > Median - [value]x IQR.  If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile
#' @param filter1probes_expressedProbes Attribute asigning an expression probe as 'expressed' if detected within at least [value]x n(individuals) at detection p-value 0.05.  If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @return A list object of class HT12prepro  where the slot with  sample-related attributes of the current processing-stage named `$chipsamples` is updated.  Excluded individual are characterized by column in_study ==F and reason4exclusion. Additionally, a slot with the detailed probe-related expression level information named `$genesdetail is created, and the slot with the history of the commands named `$history`` is updated.
#' @import data.table
#' @export

## debug
# paramfile = myparamfile
# filter1ind_expressedGenes= "from_paramfile"
# filter1probes_expressedProbes= "from_paramfile"
# ht12object =  prepro_sorbv2
# #
filterLowExpressed = function(ht12object,paramfile = NULL, filter1ind_expressedGenes= "from_paramfile", filter1probes_expressedProbes= "from_paramfile") {


### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

myparameters = match.call()
showVennplots = F


# status checken
historie =  ht12object$history$calls
if(any(grepl("createExpressionSet", historie))==F) stop("Function 'createExpressionSet()' has to be run before!")


#laden parameter
if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)


sample_overview_l5_initial <- ht12object$chipsamples
mytable(sample_overview_l5_initial$in_study)
sample_overview_l5 <- sample_overview_l5_initial[ sample_overview_l5_initial$in_study, ]
dim(sample_overview_l5_initial)
dim(sample_overview_l5)
table(table(sample_overview_l5$new_ID))
if(length(table(table(sample_overview_l5$new_ID))) != 1)
  stop("modify code, script expects ids in column new_ID!")

#laden expressionsets

all_nobkgd_eset = ht12object$all_nobkgd_eset
total_nobkgd_eset = ht12object$total_nobkgd_eset

#kontrollids laden
# annotcon_fn = "data/illumina_annot_control_beads_110407hk.txt"
# annotcon <- read.delim(annotcon_fn, as.is=T)
head(annotcon)

if('ilmn' %in% names(annotcon)) data.table::setnames(annotcon, 'ilmn', 'Probe_Id')
unique(annotcon$Reporter_Group_Name)

## ----calcoutlier---------------------------------------------------------
if(filter1ind_expressedGenes== "from_paramfile") numgen_cutoff <- as.numeric(getParam2("filter1ind_expressedGenes", myparam = param)) else numgen_cutoff = filter1ind_expressedGenes
numgen_cutoff
suppl_info <- vector(mode = 'list', length  = length(unique(sample_overview_l5$subgroup)))
names(suppl_info) <- unique(sample_overview_l5$subgroup)
# str(unique(sample_overview_l5$subgroup))
for( i in unique(sample_overview_l5$subgroup)) {
  # i = unique(sample_overview_l5$subgroup)[1]
  message("processing subgroup '", i, "' ...\n")
  bad_genenum_sub <- get_badExpressed_Inds(sample_overview_l5_subset =
                                          sample_overview_l5[sample_overview_l5$subgroup == i,, drop = F],
                                          numgen_cutoff = numgen_cutoff,
                                          subgroupname = i)
  suppl_info[[i]]['bad_genenum'] <- list(bad_genenum_sub$bad_genenum)
  suppl_info[[i]]['high_cutoff'] <- bad_genenum_sub$high_cutoff
  suppl_info[[i]]['low_cutoff'] <- bad_genenum_sub$low_cutoff
  suppl_info[[i]]['median'] <- bad_genenum_sub$median
  suppl_info[[i]]['genesfound'] <- list(bad_genenum_sub$genesfound)
  suppl_info[[i]]['IQR'] <- (bad_genenum_sub$high_cutoff - bad_genenum_sub$median)/numgen_cutoff
  # suppl_info[[i]]['genesfound_relIQR'] = list((suppl_info[[i]][['genesfound']]-bad_genenum_sub$median) / unique(suppl_info[[i]][['IQR']]))
}
bad_genenum <- unlist(lapply(suppl_info, function(x) x$bad_genenum))
bad_genenum
sample_overview_l5$failedFilter_ExprimGene <- sample_overview_l5$new_ID %in% bad_genenum
sample_overview_l5$in_study <- ifelse(sample_overview_l5$new_ID %in% bad_genenum, F, sample_overview_l5$in_study)
message("Remaining samples acrosss all subgroups:\n",sum(sample_overview_l5$in_study), ", i.e. ", ifelse(length(bad_genenum)==0, 0, length(bad_genenum)), " were filtered due to atypical number of expressed genes...")

## ----relgen--------------------------------------------------------------
sample_overview_l5 = data.table(data.frame(sample_overview_l5))
sample_overview_l5[,.N, .(in_study, reason4exclusion)]

sample_overview_l5[1]
sample_overview_l5[, Detected.Genes..0.01. := as.numeric(Detected.Genes..0.01.)]
sample_overview_l5[, mymedian := median(Detected.Genes..0.01.), by = subgroup]
sample_overview_l5[, iqr:= IQR(Detected.Genes..0.01.), by = subgroup]
sample_overview_l5[, Detected.Genes..0.01._relIQR := (mymedian - Detected.Genes..0.01.)/iqr, by = subgroup]
sample_overview_l5[, max(abs(Detected.Genes..0.01._relIQR)), by = .(in_study, subgroup)]



sample_overview_l5 = data.table(data.frame(sample_overview_l5))

sample_overview_l5[in_study == T, stopifnot(max(abs(Detected.Genes..0.01._relIQR)) <= numgen_cutoff)]
sample_overview_l5[, mymedian := NULL]
sample_overview_l5[, iqr:= NULL]
setDF(sample_overview_l5)

## ----validind------------------------------------------------------------
all_nobkgd_eset
total_nobkgd_eset # inkl kontrollinformation
# reduzieren auf die mit Expression

all_nobkgd_eset <- all_nobkgd_eset[, setdiff(pData(all_nobkgd_eset)$sampleID, bad_genenum)]
all_nobkgd_eset
total_nobkgd_eset <- total_nobkgd_eset[,setdiff(pData(total_nobkgd_eset)$sampleID, bad_genenum)]
total_nobkgd_eset

## ----probannot-----------------------------------------------------------
hh(exprs(all_nobkgd_eset))
genesdetail <- data.frame(nuid=row.names(exprs(total_nobkgd_eset)))
# str(genesdetail)

# kontrollstatus anfuegen an probenannotation
# annotcon nuid fuer ercc
ercc <- ht12object$ercc
ht(ercc,2)
par(mfrow <- c(1,1))
qlist456 <- venn3(genesdetail$nuid, ercc$Probe_Id, ercc$nuid, plotte = showVennplots)
annotcon[1,]
qlist456b <- venn3(annotcon$Probe_Id, ercc$Probe_Id, ercc$Array_Address_Id, plotte = showVennplots)
ht(annotcon,2)
annotcon$nuid <- suppressWarnings(lumi::IlluminaID2nuID(IlluminaID=annotcon$Probe_Id))[,"nuID"]
table(annotcon[is.na(annotcon$nuid), "Reporter_Group_Name"])
table(annotcon[is.na(annotcon$nuid) == F, "Reporter_Group_Name"])
annotcon[is.na(annotcon$nuid), "nuid"] <- annotcon[is.na(annotcon$nuid), "Probe_Id"]
annotcon <- rename(annotcon, c(Probe_Id = "ilmn"))
ht(annotcon, 2)
qlist456 <- venn3(genesdetail$nuid, annotcon$nuid, ercc$nuid, plotte = showVennplots)
par(mfrow = c(1,1))
housekeeping_samples <- annotcon[annotcon$Reporter_Group_Name == "housekeeping", "nuid"]
# str(housekeeping_samples)
qlist776 <- venn3(housekeeping_samples, annotcon$nuid, rownames(exprs(total_nobkgd_eset)), plotte = showVennplots)
# str(qlist776 )
genesdetail$is_purecontrol <- genesdetail$nuid %in% qlist776$q3
mytable(genesdetail$is_purecontrol)
genesdetail$is_ercc <- genesdetail$nuid %in% ercc$nuid
mytable(genesdetail$is_ercc)  ##herausnehmen der ercc-probe nicht mehr, nur filtervariable spezifizieren
xtabs_hk(~genesdetail$is_ercc + genesdetail$is_purecontrol)
genesdetail$is_housekeeping <- genesdetail$nuid %in% annotcon[annotcon$Reporter_Group_Name ==  "housekeeping", 'nuid']
mytable(genesdetail$is_housekeeping)
xtabs_hk(~genesdetail$is_housekeeping + genesdetail$is_purecontrol)
genesdetail$is_negative <- genesdetail$nuid %in% annotcon[annotcon$Reporter_Group_Name == "negative", 'nuid']
mytable(genesdetail$is_negative)
xtabs_hk(~genesdetail$is_negative + genesdetail$is_purecontrol)
ht(genesdetail)
table(table(rownames(exprs(total_nobkgd_eset))))
table(table(genesdetail$nuid))
qlist2111 <- venn2(rownames(exprs(all_nobkgd_eset)), genesdetail[genesdetail$is_purecontrol == F, "nuid"], plotte = showVennplots)
genesdetail$ilmn <- suppressWarnings(lumi::nuID2IlluminaID(genesdetail$nuid,  idType = 'Probe'))
ht(genesdetail)
genesdetail$ilmn[is.na(genesdetail$ilmn)] <- genesdetail$nuid[is.na(genesdetail$ilmn)]
if(any(is.na(genesdetail$ilmn)))
  stop("id-zuordnung lumi-ilmn failed somehow")
genesdetail <- moveColFront(genesdetail, 'ilmn')

## ----checkdirectionp-----------------------------------------------------
# str(lumi::detection(all_nobkgd_eset))
# lumi::detection(all_nobkgd_eset)[1:5,1:5]
housekeeping <- genesdetail[genesdetail$is_housekeeping, "nuid"]
housekeeping
lumi::detection(all_nobkgd_eset)[housekeeping,1:9]  #hoher lumi::detection p wert --> wird detetktiert
if(all(as.vector(lumi::detection(all_nobkgd_eset)[housekeeping,1:9]) > 0.5)==FALSE) {
  print(lumi::detection(all_nobkgd_eset)[housekeeping,1:9] )
  stop("The code assumes that high detection p-values means high expressed. Old versions of illumina do not have this. Here, housekeeping genes for the first 9 samples were checked and values smaller thatn 0.5 were found - not expected, stopping... ")} #cave, teilweise anders beschrieben, bei  paper archer_2010_briefings_in_bioinformatics_nov25_calling_algorithm_illumina.pdf

## ----calcexprtranscripts-------------------------------------------------

if(filter1probes_expressedProbes== "from_paramfile")  minexprimiert <- as.numeric(getParam2("filter1probes_expressedProbes", myparam = param)) else minexprimiert = filter1probes_expressedProbes
minexprimiert
getPresentProbes <- function (minexprimiert, all_nobkgd_eset, total_nobkgd_eset) {
  all_pvals <- as.vector(1 - lumi::detection(all_nobkgd_eset))
  #wieviel prozent der gene sind exprimiert
  check <- table(all_pvals < minexprimiert)
  check
  message("Across remaining ", dim(all_nobkgd_eset)[2], " samples from subgroup ", unique(sample_overview_l5[ sample_overview_l5$new_ID %in% Biobase::sampleNames(all_nobkgd_eset), "subgroup"]), ", ", proz(check/(sum(check)))[2], " probes are classified as expressed...")
  length(all_pvals)
  par(mfrow  = c(1,1))
  mytitle = paste0("Distribution 1-lumi::detection-values \n ",
                   dim(all_nobkgd_eset)[1],
                   " probes x ",
                   dim(all_nobkgd_eset)[2], " inds ")
  hist(1 - lumi::detection(all_nobkgd_eset), breaks = 100, main = mytitle)
  abline(v = 0.05, col = "red")

  # herausfinden  welche Transkripte in mehr als x% der Individuen transkribiert sind
  minind <-  minexprimiert*dim(total_nobkgd_eset)[2]
  minind <- ceiling(minind)
  minind

  # Aufschreiben der individuenzahl, in denen das jeweilige Gen transkribiert wird
  present <-  apply(lumi::detection(total_nobkgd_eset), 1, function(x) sum(x>=(1 - minexprimiert)))
  expressed <- present >= minind
  res <- c()
  res$present <- present
  res$expressed <- expressed
  res$minind <- minind
  res
}
unique(sample_overview_l5$subgroup)
for(i in unique(sample_overview_l5$subgroup)) {
  # i = unique(sample_overview_l5$subgroup)[1]
  # message("processing subgroup '", i, "' ...\n")
  subind <- sample_overview_l5[ sample_overview_l5$failedFilter_ExprimGene == F & sample_overview_l5$subgroup == i, 'new_ID']
  newdata <- getPresentProbes(minexprimiert, all_nobkgd_eset[,subind], total_nobkgd_eset[, subind])
  var_present <- paste0("present_", reformate_subgroup(i))
  genesdetail[var_present] <- newdata$present
  var_expressed <- paste0("expressed_", reformate_subgroup(i))
  genesdetail[var_expressed] <- newdata$expressed
  suppl_info[[i]]['minind'] <- newdata$minind
}
ht(genesdetail, 1)


## ----indfind-------------------------------------------------------------
# kleiner test
# dd = detectframe[1:5, 1:6]
# dd
# ee= exprmatrix[1:5, 1:6]
# ee
#
# ee2 = dd * ee
# controlzahl = apply(ee2, 1, function(x) sum(x!= 0))
# controlzahl
# ee2[ ee2 ==0] = NA
# ee2
# rowMeans(ee2, na.rm = T)
calcMeaneEtc <- function (minexprimiert, total_nobkgd_eset, genesdetail) {
  detectframe <- lumi::detection(total_nobkgd_eset)

  # detectframe[1:5,1:5]
  detectframe[detectframe>=(1-minexprimiert)] <- 1
  detectframe[detectframe<(1-minexprimiert)] <- 0

  # detectframe[1:5,1:5]
  exprmatrix <- exprs(total_nobkgd_eset)
  exprmatrix2 <- detectframe * exprmatrix
  controlzahl <- apply(exprmatrix2, 1, function(x) sum(x!= 0))
  exprmatrix2[ exprmatrix2 == 0] = NA
  meane <- apply(exprmatrix2, 1, function(x) ifelse(all(is.na(x)) == T, NA, mean(x, na.rm = T)))

  # calc relativ expression
  negativ_samples <- genesdetail[genesdetail$is_negative, "nuid"]
  # str(negativ_samples)
  negative_sometimesmissing <- setdiff(negativ_samples, row.names(exprs(total_nobkgd_eset)))
  negative_sometimesmissing
  negativ_present <- setdiff(negativ_samples, negative_sometimesmissing)
  length(negativ_present)
  negativeMean <- colMeans(exprs(total_nobkgd_eset[negativ_present,]))
  # str(negativeMean)
  exprmatrix2_rel <- t(t(exprmatrix2)/negativeMean)

  # exprmatrix2_rel[1:5,1:5]
  # exprs(total_nobkgd_eset)[1:5,1:5]
  # exprmatrix2[1:5,1:5]
  # negativeMean[1:5]         #haut hin
  MeanRatioExpressed <- apply(exprmatrix2_rel, 1, function(x) ifelse(all(is.na(x)) == T, NA,mean(x, na.rm = T)))
  # sm =  apply(exprmatrix2_rel, 1, function(x) ifelse(all(is.na(x)) == T, -1, length(na.omit(x))))
  # table(sm)
  MaxRatioExpressed <- suppressWarnings(apply(exprmatrix2_rel, 1, function(x) ifelse(all(is.na(x)) == T, NA, max(x, na.rm = T)))) #TODO warning aufloesen
  res <- c()
  res$meane <- meane
  res$controlzahl <- controlzahl
  res$negative_sometimesmissing <- negative_sometimesmissing
  res$MeanRatioExpressed <- MeanRatioExpressed
  res$MaxRatioExpressed <- MaxRatioExpressed
  res
}

for(i in unique(sample_overview_l5$subgroup)) {
  # i = unique(sample_overview_l5$subgroup)[1]
  # message("processing subgroup '", i, "' ...\n")
  subind <- sample_overview_l5[ sample_overview_l5$failedFilter_ExprimGene == F & sample_overview_l5$subgroup == i, 'new_ID']

  # mean ueber exprimierte
  newmeaneEtc <- calcMeaneEtc(minexprimiert, total_nobkgd_eset[,subind], genesdetail)
  stopifnot(all(genesdetail[, paste0('present_', reformate_subgroup(i))] == newmeaneEtc$controlzahl))
  var_meane <- paste0("meane_", reformate_subgroup(i))
  genesdetail[var_meane] <- newmeaneEtc$meane
  var_MeanRatioExpressed <- paste0("MeanRatioExpressed_", reformate_subgroup(i))
  genesdetail[var_MeanRatioExpressed] <- newmeaneEtc$MeanRatioExpressed
  var_MaxRatioExpressed <- paste0("MaxRatioExpressed_", reformate_subgroup(i))
  genesdetail[var_MaxRatioExpressed] <- newmeaneEtc$MaxRatioExpressed

  # mean ueber alle
  var_meanall <- paste0("meanall_", reformate_subgroup(i))
  genesdetail[var_meanall] <- rowMeans(exprs(total_nobkgd_eset[,subind]), na.rm=T)

  # dokuinfo
  if(exists('negative_sometimesmissing') )  {
    negative_sometimesmissing_old <- negative_sometimesmissing
    negative_sometimesmissing <- newmeaneEtc$negative_sometimesmissing
    stopifnot(negative_sometimesmissing_old  == negative_sometimesmissing)
  } else negative_sometimesmissing <- newmeaneEtc$negative_sometimesmissing
  suppl_info$negative_sometimesmissing <- negative_sometimesmissing
}
ht(genesdetail, 3)

## ----desind--------------------------------------------------------------
par(mfrow=c(1,1))
for(i in unique(sample_overview_l5$subgroup)) {
  # i = unique(sample_overview_l5$subgroup)[1]
  # message("processing subgroup '", i, "' ...\n")
  subind <- sample_overview_l5[sample_overview_l5$failedFilter_ExprimGene == F & sample_overview_l5$subgroup == i, 'new_ID']

  # mean ueber exprimierte
  var_present <- paste0("present_", reformate_subgroup(i))
  hist(genesdetail[genesdetail$is_purecontrol == F, var_present],
       breaks = 100,
       main = paste0("Number of probes, that are classified 'expressed' \nwithin the number of individuales shown at the x-axis (Total probes: ",
                   dim(all_nobkgd_eset)[1], ")"), cex.main=0.9)
  minind_sub <- suppl_info[[i]]$minind
  mtext(paste0("red: 'is_expressed'-cutoff: Probe had to be expressed in at lest ",
               minind_sub,
               " (",
               minexprimiert*100, "%) \nof all ",
               dim(all_nobkgd_eset)[2],
               " samples (Subgroup '", i, "')"),
        cex = 0.8, line = -0.7)
abline(v = minind_sub, col = "red")
}

## ----specprobes----------------------------------------------------------
# ercc vs. negatives vs. housekeeping
negativ_samples <- genesdetail[genesdetail$is_negative, "nuid"]
plotneg_ercc <- plyr::ddply(sample_overview_l5, 'subgroup', .progress = "text", function(df){

  # df= sample_overview_l5[ sample_overview_l5$subgroup ==sample_overview_l5$subgroup[1],]
  df <- df[df$failedFilter_ExprimGene == F, ]
  ercc_expr <- as.vector((exprs(total_nobkgd_eset[ercc$nuid,df$new_ID])))
  negative_expr <- as.vector((exprs(total_nobkgd_eset[negativ_samples, df$new_ID])))

  # str(negative_expr)
  housekeeping_expr <- as.vector((exprs(total_nobkgd_eset[ housekeeping_samples,df$new_ID])))
  plotneg_ercc_sub <- rbind(data.frame(gruppe = "ercc",
                                       raw_expression = ercc_expr),
                            data.frame(gruppe = "negative",
                                       raw_expression = negative_expr),
                            data.frame(gruppe = "housekeeping",
                                       raw_expression = housekeeping_expr))
  plotneg_ercc_sub
  }
)
ht(plotneg_ercc)
plotneg_ercc_bild <- ggplot(plotneg_ercc, aes(raw_expression, fill = gruppe)) +
  geom_density(alpha = 0.4) +
  scale_x_log10() +
  facet_grid(gruppe ~ ., scales='free_y') +
  facet_grid( gruppe ~ subgroup, scale = "free_y")
plotneg_ercc_bild

## ----save----------------------------------------------------------------

## speichern
ht(genesdetail, 3)

ht12object$genesdetail = genesdetail

# samples as used
sample_overview_l6doku <- sample_overview_l5_initial
qlist543 <- venn2(names(sample_overview_l5), names(sample_overview_l5_initial), plotte = showVennplots)
qlist543
for(i in qlist543$q2) {
  sample_overview_l6doku[i] <- sample_overview_l5[match_hk(sample_overview_l6doku$new_ID,
                                                           sample_overview_l5$new_ID), i]

}
sample_overview_l6doku[sample_overview_l6doku$new_ID %in% bad_genenum,
                       "in_study"] = F
sample_overview_l6doku[sample_overview_l6doku$new_ID %in% bad_genenum,
                       "reason4exclusion"] = "filterLowExpressed: extreme number expressed genes"
xtabs_hk(~sample_overview_l6doku$reason4exclusion + sample_overview_l6doku$subgroup)
mytable(sample_overview_l6doku$in_study)

sample_overview_l6 <- sample_overview_l6doku

ht12object$chipsamples = sample_overview_l6

# status fuer docku
dim_all_nobkgd_eset <- dim(all_nobkgd_eset)
fordoku <- c('suppl_info',
                       "plotneg_ercc_bild",
                       "genesdetail",
                       "sample_overview_l5",
                       "sample_overview_l6",
                       "numgen_cutoff",
                       "bad_genenum",
                       "minexprimiert",
                       "dim_all_nobkgd_eset",
             'filter1ind_expressedGenes', "filter1probes_expressedProbes")


stopifnot(sum(duplicated(fordoku))==0)
ht12object$dokuobjects_filterLowExpressed = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_filterLowExpressed) = fordoku
# str(ht12object$dokuobjects_filterLowExpressed)
ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object

}
