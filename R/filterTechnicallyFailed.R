#' @title Identify and filter samples with unusual combination of control probes in a HT12prepro object
#' @description log2 transformed and normalized  not heavily correlated control probes of the still valid samples are used to create a correlation-adjusted score (via the Mahalanobis distance of these parameters) reflecting technical accuracy of the expression measurement. Outlyers are identified defined as having a greater distance than `filter2ind_atypischIlmnKontroll` interquantile ranges. See vignette for an example.

#' @param ht12object A list object of class HT12prepro created with function transformNormalizeHT12object()
#' @param paramfile Path to the file specifying parameters
#' @param controlparameters2check Illuminas control features used for calculation of the Mahalanobis distance based outlyer criterium. Default ist 'hybrid_low, hybrid_med,  string_pm, string_mm, biotin, negative' (if spiked in erccc is used in ALL CHIPS use also 'ercc_1, ercc_2, ercc_3, ercc_4, ercc_5', if spked in artificial polyadenylated RNAs from Bacillus subtilis immediately preceding the reverse transcription step is used, also 'labeling'). In case of population based studies within a single tissue default is also 'housekeeping, Detected.Genes..0.01.'.  If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @param robustmethod_for_mahal Method used for calculation of robust correllation for Mahalanobis distance within function mdqc(). Previous default was "S-estimator" (which has a 25 percent breakdown point), now "MVE" is prefferred (i.e. the Minimum Volume Ellipsoid (MVE) which searches for the ellipsoid with the smallest volume that covers h data points). Alternative is "MCD" (i.e. the Minimum Covariance Determinant (MCD) estimator which looks for the subset of h data points whose covariance matrix has the smallest determinant). If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @param filter2ind_atypischIlmnKontroll Filter for extreme combination of Illuminas control features summarized as Mahalanobis distance to an artificial sample. This artificial sample having has average values for the selected quality control features. Valid is ln('Mahalanobis distance') < median(ln('Mahalanobis distance')) + [value] *  IQR(ln('Mahalanobis distance')). If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
#' @return A list object of class HT12prepro including updated slot `$chipsamples` with  sample-related attributes of the current processing-stage. Excluded individual are characterized by column in_study ==F and reason4exclusion. QC plots are shown.
#' @import data.table
#' @import mdqc
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# filter2ind_atypischIlmnKontroll = "from_paramfile"
# controlparameters2check = "from_paramfile"
# robustmethod_for_mahal = "from_paramfile"
# ht12object =  prepro_ht12


filterTechnicallyFailed = function(ht12object,paramfile = NULL,filter2ind_atypischIlmnKontroll = "from_paramfile", controlparameters2check = "from_paramfile", robustmethod_for_mahal = "from_paramfile") {

### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

myparameters = match.call()
showVennplots = F

# status checken
historie =  ht12object$history$calls
if(any(grepl("transformNormalizeHT12object", historie))==F) stop("Function 'transformNormalizeHT12object()' has to be run before!")

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

# laden annotation probes
genesdetail <- ht12object$genesdetail
subgroups <- unique(sample_overview_l6instudy$subgroup)
subgroups

# laden expressionsets
total_nobkgd_eset_ql =  ht12object$total_nobkgd_eset_ql

# annotcon
# kontrollids laden
head(annotcon)
if('ilmn' %in% names(annotcon)) data.table::setnames(annotcon, 'ilmn', 'Probe_Id')

unique(annotcon$Reporter_Group_Name)
annotcon$nuid <- suppressWarnings(lumi::IlluminaID2nuID(IlluminaID = annotcon$Probe_Id))[, "nuID"]
annotcon[is.na(annotcon$nuid), "nuid"] = annotcon[is.na(annotcon$nuid), "Probe_Id"]
annotcon <- rename(annotcon, c(Probe_Id = "ilmn"))
ht(annotcon, 2)

# details ercc laden
ht(ercc_det, 2)

## ----erccdetsilis--------------------------------------------------------
# plot(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.) ~ log10(ercc_det$concentration.in.Mix.2..attomoles.ul.))
qlist2 <- venn2(annotcon$Reporter_Group_Name, ercc_det$ERCC.ID, plotte = showVennplots)
matcher <- annotcon[annotcon$Reporter_Group_Name %in% ercc_det$ERCC.ID,]
ercc_det$ilmn <- matcher[match_hk(ercc_det$ERCC.ID, matcher$Reporter_Group_Name), "ilmn"]
gx_ercc <- exprs(total_nobkgd_eset_ql[genesdetail[genesdetail$is_ercc, "nuid"], ])
hh(gx_ercc)
rownames(gx_ercc) <- lumi::nuID2IlluminaID(rownames(gx_ercc))[, 'Probe_Id']
hh(gx_ercc)
ercc_dt <- data.frame(gx_ercc)
inds_in_ercc_dt <- names(ercc_dt)
ercc_dt$ilmn <- rownames(ercc_dt)
rownames(ercc_dt) <- NULL
ercc_dt <- moveColFront(ercc_dt, "ilmn")
hh(ercc_dt)
ercc_det$ercc_gruppe1 <- cut(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.),
                             breaks = c(-3, 0, 1, 2, 3, 4, 5), include.lowest = T)
mytable(ercc_det$ercc_gruppe1)
# boxplot(log10(ercc_det$concentration.in.Mix.1..attomoles.ul.) ~ercc_det$ercc_gruppe1)
ercc_det$ercc_gruppe2 <- cut(log10(ercc_det$concentration.in.Mix.2..attomoles.ul.),
                             breaks = c(-3, 0, 1, 2, 3, 4, 5), include.lowest = T)
mytable(ercc_det$ercc_gruppe2)
# boxplot(log10(ercc_det$concentration.in.Mix.2..attomoles.ul.) ~ercc_det$ercc_gruppe2)
for(i in setdiff(names(ercc_det), c("ilmn"))){
  ercc_dt[i] <- ercc_det[match_hk(ercc_dt$ilmn, ercc_det$ilmn), i]
}


ercc_det$ercc_gruppe1short <- paste0("ercc_", stringr::str_sub(sapply(stringr::str_split(ercc_det$ercc_gruppe1, ","), "[", 2), 1, 1))
xtabs(~ ercc_det$ercc_gruppe1short + ercc_det$ercc_gruppe1)

## ----read----------------------------------------------------------------
# noncorr = c( "hybrid_low", "string_pm", "string_mm", "biotin", "housekeeping",  "negative", "Detected.Genes..0.01.", "euklid", "labeling", "ercc","hybrid_med" )
if(controlparameters2check== "from_paramfile") noncorr <- getParam2("controlparameters2check", myparam = param) else noncorr = filter1ind_expressedGenes
noncorr
noncorr <- stringr::str_trim(unlist(stringr::str_split(noncorr, pattern = ",")))
noncorr
message("Using:\n", paste(noncorr, collapse = "\n"), " \nas technical attributes for the technical score...\n")
noncorr_full <- c( "hybrid_low", "hybrid_med", "hybrid_high",
                   "string_pm", "string_mm", "biotin",
                   "housekeeping",  "negative","labeling",
                   "ercc", "Detected.Genes..0.01.", "Detected.Genes..0.05.",
                   "samplemean",  paste0("ercc_",1:5))

# alle moeglichen, auch wenn quantile und mean bereits durch quantilenormalisierung abgedeckt sind,
# daher auch  samplemean bestraft die nicht-normalitaet der daten. Not that string_pm is actually
# average of (hybrid_med + hybrid_high).
if(any(noncorr %nin% noncorr_full))
  stop("Parameter ", setdiff(noncorr, noncorr_full),
       "not in valid strings for Parameter 'controlparameters2check'. These are \n",
       paste(noncorr_full, collapse = "\n"))

## ----process6------------------------------------------------------------
goodind  <- sample_overview_l6instudy$new_ID
# str(goodind)
total_nobkgd_eset_ql
total_nobkgd_eset_ql_goodind <- total_nobkgd_eset_ql[, goodind]
# total_nobkgd_eset_ql_goodind <- total_nobkgd_eset_ql_goodind[,goodind[1:50]]  # debug
total_nobkgd_eset_ql_goodind



## ----investigate---------------------------------------------------------
# str(annotcon)
unique(annotcon$Reporter_Group_Name)
unique(annotcon$Reporter_Group_id)
annotcon[annotcon$Reporter_Group_Name == "cy3_hyb", ]
table(table(annotcon$ilmn))
table(table(annotcon$Array_Address_Id))
groups <- (annotcon$Reporter_Group_Name)
table(groups[grep("ERCC", groups, invert = T)])
length(grep("^ERCC", unique(annotcon$Reporter_Group_id)))

# mittelwerte aus normalisierten file herauslesen auslesen und zugaenglich machen
# normalisierte Kontrollwerte holen
par(mfrow = c(1,1))
qlist55 <- venn2(annotcon$nuid, rownames(exprs(total_nobkgd_eset_ql)), plotte = showVennplots)
# str(qlist55)
expr_con <- t(exprs(total_nobkgd_eset_ql[qlist55$q1,]))
probenames <- colnames(expr_con)
newprobenames <- genesdetail[match_hk(probenames, genesdetail$nuid), "ilmn"]
ht(cbind(probenames, newprobenames))
if(any(is.na(newprobenames)))
  stop("Some probe assignment did not work, check code on tag 147_asdfun")
table(table(newprobenames))
colnames(expr_con) <- newprobenames
expr_con <- data.frame(new_ID = row.names(expr_con), expr_con)
hh(expr_con)
dim(expr_con)
if(any(is.na(expr_con$new_ID)))
  stop("Some probe assignment did not work, check code on tag 147_asdffasfun")

## ----newmatrices---------------------------------------------------------
# ercc details
expr_con$ercc_1 <- apply(as.matrix(expr_con[,ercc_det[ercc_det$ercc_gruppe1short == "ercc_1","ilmn"]]), 1, mean)
expr_con$ercc_2 <- apply(as.matrix(expr_con[,ercc_det[ercc_det$ercc_gruppe1short == "ercc_2","ilmn"]]), 1, mean)
expr_con$ercc_3 <- apply(as.matrix(expr_con[,ercc_det[ercc_det$ercc_gruppe1short == "ercc_3","ilmn"]]), 1, mean)
expr_con$ercc_4 <- apply(as.matrix(expr_con[,ercc_det[ercc_det$ercc_gruppe1short == "ercc_4","ilmn"]]), 1, mean)
expr_con$ercc_5 <- apply(as.matrix(expr_con[,ercc_det[ercc_det$ercc_gruppe1short == "ercc_5","ilmn"]]), 1, mean)

# hybridisation  = cy3_hyb
expr_con$hybrid_low <- apply(as.matrix(expr_con[,c("ILMN_1343050", "ILMN_1343052")]), 1, mean)
expr_con$hybrid_med <- apply(as.matrix(expr_con[,c("ILMN_2038768", "ILMN_2038771")]), 1, mean)
expr_con$hybrid_high <- apply(as.matrix(expr_con[,c("ILMN_2038769", "ILMN_2038770")]), 1, mean)

# plotting
myplot <- reshape(expr_con[, c("new_ID","hybrid_low", "hybrid_med", "hybrid_high")],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = c("hybrid_low",
                            "hybrid_med",
                            "hybrid_high"),
                  timevar = "bybridisationcategori",
                  varying = list(c("hybrid_low", "hybrid_med", "hybrid_high")))
head(myplot)
tail(myplot)
myplot$bybridisationcategori <- stringr::str_replace(myplot$bybridisationcategori, c("hybrid_low"), c("hybrid_1_low"))
myplot$bybridisationcategori <- stringr::str_replace(myplot$bybridisationcategori, c("hybrid_med"), c("hybrid_2_med"))
myplot$bybridisationcategori <- stringr::str_replace(myplot$bybridisationcategori, c("hybrid_high"), c("hybrid_3_high"))
dim(myplot)
myplot <- myplot[myplot$new_ID %in% goodind,]
dim(myplot)
hybridisationsplot <- lattice::bwplot(hybrid_low~factor(bybridisationcategori),
                                      data = myplot,
                                      xlab = "hybybridisation levels",
                                      ylab = "average expression",
                                      main = "illumina qc - overall cy3 hyb")
# plot(hybridisationsplot)
expr_con$hybrid_medrel <- (expr_con$hybrid_med - expr_con$hybrid_low)
expr_con$hybrid_highrel <- (expr_con$hybrid_high - expr_con$hybrid_low)

# plotting
myplot <- reshape(expr_con[, c("new_ID", "hybrid_medrel", "hybrid_highrel")],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = c("hybrid_medrel", "hybrid_highrel"),
                  timevar="rel_hybridisation",
                  varying=list(c("hybrid_medrel", "hybrid_highrel")))
head(myplot)
tail(myplot)
hybridplot_rel <- lattice::bwplot(hybrid_medrel~factor(rel_hybridisation),
                                  data = myplot,
                                  xlab = "rel hybybridisation levels",
                                  ylab = "average expression",
                                  main = "illumina qc - overall cy3 hyb Diff log2 Werte")
# plot(hybridplot_rel)

# mismatch vs perfectmatch = low_stringency_hyb
expr_con$string_pm <- apply(as.matrix(expr_con[,c("ILMN_2038769", "ILMN_2038770", "ILMN_2038768", "ILMN_2038771")]), 1, mean)
expr_con$string_mm <- apply(as.matrix(expr_con[,c("ILMN_1343061", "ILMN_1343063", "ILMN_1343062", "ILMN_1343064")]), 1, mean)

# plotting
myplot <- reshape(expr_con[, c("new_ID", "string_pm", "string_mm")],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = c("string_pm", "string_mm"),
                  timevar = "mm_vs_pm",
                  varying = list(c("string_pm", "string_mm")))
head(myplot)
tail(myplot)
perfectmismatchplot <- lattice::bwplot(string_pm~factor(mm_vs_pm),
                                       data = myplot,
                                       xlab = "mm_vs_pm",
                                       ylab = "average expression",
                                       main = "illumina qc - overall mm_vs_pm")
# plot(perfectmismatchplot)
expr_con$string_pmrel <- (expr_con$string_pm - expr_con$string_mm)
par(mfrow = c(1,1))
# boxplot(expr_con$string_pmrel,
# ylab = "(expr_con$string_pm - expr_con$string_mm)",
# main = "stringency_controls_pm_mm_relative")

## biotin
get_mean_expressionlevel <- function (congruppe) {
  congruppe_name = deparse(substitute(congruppe))
  specific_samples <- annotcon[grep(congruppe, annotcon$Reporter_Group_Name, ignore.case = T), "ilmn"]
  message("Info: - control-category ", congruppe_name, " ", length(specific_samples), " includes following probes: ",paste(specific_samples, collapse = ", "), "\n")
  apply(as.matrix(expr_con[, specific_samples]), 1, mean)
}
expr_con$biotin <- get_mean_expressionlevel("biotin")

# housekeeping
expr_con$housekeeping <- get_mean_expressionlevel("housekeeping")

# labeling
expr_con$labeling <- get_mean_expressionlevel("labeling")

# ercc
expr_con$ercc <-  get_mean_expressionlevel("ercc")

# negative
specific_samples <- annotcon[annotcon$Reporter_Group_Name == "negative", "ilmn"]
specific_samples
negative_sometimesmissing <- setdiff(specific_samples, names(expr_con))
negative_sometimesmissing
specific_samples <- setdiff(specific_samples, negative_sometimesmissing)
length(specific_samples)
na_check <- apply(as.matrix(expr_con[, specific_samples]), 2, function(x) any(is.na(x) == T))
table(na_check)
expr_con$negative <- apply(as.matrix(expr_con[, specific_samples]), 1, mean)

# quantile und mittelwert
#checken, dass gleiche reihenfolge existiert
if(identical(as.character(colnames(exprs(total_nobkgd_eset_ql))), as.character(row.names(expr_con))) != T)
  stop("sample_overview_l6 leute aus eset muessen gematcht werden, reihenfolge sonst falsch")

#expression set bauen, was keine kontroll-DNA entfaellt und als exprimiert gilt
total_nobkgd_eset_ql

# kontrollen herauswerfen
expressed_cols <- grep("^expressed_", names(genesdetail), value = T)
expressed_cols
genesdetail$expressed_in_any_subgroup <- apply(genesdetail[, expressed_cols, drop = F], 1,
                                               function(x) any(x == T))
goodprobesNocons <- genesdetail [genesdetail$is_purecontrol == F & genesdetail$expressed_in_any_subgroup, "nuid"]
qlist55 <- venn2(goodprobesNocons, rownames(exprs(total_nobkgd_eset_ql)), plotte = showVennplots)

all_nobkgd_expressed_eset_ql <- total_nobkgd_eset_ql[goodprobesNocons, ]
all_nobkgd_expressed_eset_ql

# mean berechnen
dim(exprs(all_nobkgd_expressed_eset_ql))
expr_con$samplemean <- esApply(all_nobkgd_expressed_eset_ql, 2, mean)

# quantile berechnen
expr_con$p05 <- esApply(all_nobkgd_expressed_eset_ql, 2, function(x) quantile(x, probs = 0.05, names = F))
expr_con$p25 <- esApply(all_nobkgd_expressed_eset_ql, 2, function(x) quantile(x, probs = 0.25, names = F))
expr_con$p50 <- esApply(all_nobkgd_expressed_eset_ql, 2, function(x) quantile(x, probs = 0.50, names = F))
expr_con$p75 <- esApply(all_nobkgd_expressed_eset_ql, 2, function(x) quantile(x, probs = 0.75, names = F))
expr_con$p95 <- esApply(all_nobkgd_expressed_eset_ql, 2, function(x) quantile(x, probs = 0.95, names = F))
summary(expr_con$p05)
summary(expr_con$p25)
summary(expr_con$p50)
summary(expr_con$p75)
summary(expr_con$p95)

# daten von existierender qc hinzufuegen - detected genes
dim(sample_overview_l6instudy)
names(sample_overview_l6instudy)
dim(expr_con)
hh(expr_con)
expr_con <- merge(expr_con,
                  sample_overview_l6instudy[,c("new_ID", "Detected.Genes..0.01.","Detected.Genes..0.05.")],
                  by.x = "new_ID",
                  by.y = "new_ID",
                  all.x = T,
                  sort = F) # change statt 0.05
dim(expr_con)
summary(expr_con$Detected.Genes..0.01.)

# noise negative background labelling housekeeping all genes
expr_con$housekeeping2neg <- (expr_con$housekeeping - expr_con$negative)
if(any(is.na(expr_con$housekeeping2neg)))
  stop("NA bei berechnung")
expr_con$samplemean2neg <- (expr_con$samplemean - expr_con$negative)
if(any(is.na(expr_con$samplemean2neg)))
  stop("NA bei berechnung")
expr_con$labeling2neg <- (expr_con$labeling - expr_con$negative)
if(any(is.na(expr_con$labeling2neg)))
  stop("NA bei berechnung")
expr_con$ercc2neg <- (expr_con$ercc - expr_con$negative)
if(any(is.na(expr_con$ercc2neg)))
  stop("NA bei berechnung")
expr_con$p952neg <- (expr_con$p95 - expr_con$negative)
if(any(is.na(expr_con$p952neg)))
  stop("NA bei berechnung")
expr_con$p052neg <- (expr_con$p05 - expr_con$negative)
if(any(is.na(expr_con$p052neg)))
  stop("NA bei berechnung")
expr_con$p952p05 <- (expr_con$p95 / expr_con$p05)
if(any(is.na(expr_con$p952p05)))
  stop("NA bei berechnung")
expr_con$biotin2neg <- (expr_con$biotin - expr_con$negative)
if(any(is.na(expr_con$biotin2neg)))
  stop("NA bei berechnung")

## ----plotting------------------------------------------------------------
grep("2neg", names(expr_con), value = T)
plotcategories <- c("biotin2neg",
                    "housekeeping2neg",
                    "samplemean2neg",
                    "ercc2neg",
                    "labeling2neg",
                    "p952neg",
                    "p952p05",
                    "p052neg")
myplot <- reshape(expr_con[, c("new_ID",plotcategories)],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = plotcategories,
                  timevar = "exprelevels",
                  varying = list(plotcategories))
head(myplot)
tail(myplot)
expressionintensplot_rel <- lattice::bwplot(biotin2neg ~ factor(exprelevels),
                                            data = myplot,
                                            xlab = "genecat",
                                            ylab = " expression ratios",
                                            main = "illumina qc - expression in subclasses",
                                            scales = list(x = list(rot = 20)))
# plot(expressionintensplot_rel)

# plottin ohne bezug zu negativ
plotcategories <- c("biotin",
                    "housekeeping",
                    "samplemean",
                    "ercc",
                    "labeling",
                    "p95",
                    "p95",
                    "p05",
                    "negative")
myplot <- reshape(expr_con[, c("new_ID",plotcategories)],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = plotcategories,
                  timevar = "exprelevels",
                  varying = list(plotcategories))
head(myplot)
tail(myplot)
expressionintensplot_abs <- lattice::bwplot(biotin ~ factor(exprelevels),
                                            data = myplot,
                                            xlab = "genecat",
                                            ylab = "expression",
                                            ain = "illumina qc - expression in subclasses",
                                            scales = list(x = list(rot = 20)))
# plot(expressionintensplot_abs)

# allequantile quantile
plotcategories <- c("p95",
                    "p75",
                    "p50",
                    "samplemean",
                    "p25",
                    "p05",
                    "negative")
myplot <- reshape(expr_con[, c("new_ID",plotcategories)],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs ,
                  times = plotcategories,
                  timevar = "exprelevels",
                  varying = list(plotcategories))
head(myplot)
tail(myplot)
quantilplot <- lattice::bwplot(p95 ~ factor(exprelevels),
                               data = myplot,
                               xlab = "genecat",
                               ylab = "expression",
                               main = "illumina qc - expression in subclasses",
                               scales = list(x = list(rot = 20)))
# plot(quantilplot)

## ----calcMahal5, fig.width=12, fig.height=12-----------------------------
# covarianz
# liste aller maasse erzeugen ohne die relativen massen
names_save <- names(expr_con)[!names(expr_con) %in% grep("^ILMN", names(expr_con), value = T)]
names_save <- names_save[grep("2neg", names_save, invert = T)]
names_save <- names_save[grep("2p05", names_save, invert = T)]
names_save <- names_save[grep("rel$", names_save, invert = T)]
names_save

# a = unique(annotcon$Reporter_Group_id)
# names_save = unique(annotcon[annotcon$Reporter_Group_id %in% c("phage_lambda_genome:mm2",
# "housekeeping", "phage_lambda_genome", "phage_lambda_genome:high", "phage_lambda_genome:pm",
# "phage_lambda_genome:low", "phage_lambda_genome:med"),"ilmn"])
expr_con_save <- expr_con[, names_save]
# str(expr_con_save)
hh(expr_con_save)
row.names(expr_con_save) <- expr_con_save$new_ID

# cov berechnen
expr_con_save_n <- expr_con_save
expr_con_save_n$new_ID <- NULL
expr_con_save_n <- scale(expr_con_save_n)
expr_con_save_n_cov <- abs(cov(expr_con_save_n))
dim(expr_con_save_n_cov)
dim(expr_con_save_n)

## ----compadd-------------------------------------------------------------
newQC <- setdiff(names_save, names(sample_overview_l6))
newQC
dim(sample_overview_l6)
sample_overview_l6save <- sample_overview_l6
sample_overview_l6 <- merge(sample_overview_l6,
                            expr_con_save[,c("new_ID", newQC)],
                            by.x="new_ID",
                            by.y="new_ID",
                            all.x = T,
                            sort = F)
setnames(sample_overview_l6, newQC, paste0(newQC, "_s06"))
class(sample_overview_l6)
ht(sample_overview_l6, 1)
dim(sample_overview_l6)
names(sample_overview_l6)
# plot(sample_overview_l6[grep("biot", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("p05", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("p25", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("p50", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("p75", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("p95", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("negat", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("label", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])
# plot(sample_overview_l6[grep("house", names(sample_overview_l6), ignore.case = T)][sample_overview_l6$in_study, ])

## ----covid, fig.width = 12, fig.height = 12------------------------------
hh(expr_con_save_n_cov)
mexpr_con_save_n_cov <- melt(expr_con_save_n_cov)
ht(mexpr_con_save_n_cov,4)
mexpr_con_save_n_cov <- plyr::ddply(mexpr_con_save_n_cov, plyr::.(Var2), transform, rescale = scales::rescale(value))
allecor_plot <- ggplot2::ggplot(mexpr_con_save_n_cov, ggplot2::aes(Var1, Var2, label = round(value,2))) +
  ggplot2::geom_tile(ggplot2::aes(fill = rescale), colour = "white") +
  ggplot2::scale_fill_gradient(low = "steelblue", high = "yellow") +
  ggplot2::geom_text(size=I(5)) +
  ggplot2::theme(axis.text.x = element_text(angle= -90, hjust = 0, size = 12), axis.text.y = element_text(size = 12)) +
  ggplot2::ggtitle("Kovaraianz Illumina QC Paramter")
# allecor_plot

# verwendete QC Parameter:
noncorr
phagenprobe <- annotcon[grep("phage", annotcon$Reporter_Group_id, ignore.case = T), c("ilmn", "Reporter_Group_Name", "Reporter_Group_id")]
phagenprobe$dupliilmn <- duplicated(phagenprobe$ilmn)
phagenprobe <- phagenprobe[order(phagenprobe$Reporter_Group_id, phagenprobe$ilmn), ]
phagenprobe

# Not that string_pm is actually average of (hybrid_med + hybrid_high).
# I dropped hybrid high as this is after quantilenormalisation way above
# expression levels of relevant samples.
noncorr
subset <- mexpr_con_save_n_cov[mexpr_con_save_n_cov$Var1 %in% noncorr & mexpr_con_save_n_cov$Var2 %in% noncorr, ]
plotcor2 <- ggplot2::ggplot(subset, ggplot2::aes(Var1, Var2, label = round(value,3))) +
  ggplot2::geom_tile(ggplot2::aes(fill = rescale), colour = "white") +
  ggplot2::scale_fill_gradient(low = "steelblue", high = "yellow") +
  ggplot2::geom_text(size = I(5)) +
  ggplot2::theme(axis.text.x = element_text(angle = -90, hjust = 0, size = 12),axis.text.y = element_text(size = 12)) +
  ggplot2::ggtitle("Covariance of included Illumina control-QC Paramters") + ggplot2::guides(fill = F)
plot(plotcor2 )

# allequantile quantile
setdiff(noncorr_full, names(expr_con))
plotcategories <- c(noncorr_full, "p95", "p75", "p50", "p25", "p05")
myplot <- reshape(expr_con[, c("new_ID", plotcategories)],
                  direction = "long",
                  idvar = "new_ID",
                  ids = expr_con$new_IDs,
                  times = plotcategories,
                  timevar = "exprelevels",
                  varying = list(plotcategories))
names(myplot) = c('new_ID', 'variable', 'value')
ht(myplot)

allcategplot2 <- ggplot2::ggplot(myplot[myplot$variable %in% noncorr,],
                                 ggplot2::aes(variable, value, fill = variable)) +
  ggplot2::geom_boxplot() +
  ggplot2::theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.1, size = 12),
                 axis.text.y = element_text(size = 12)) + ggplot2::guides(fill = F) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(10)) + ggplot2::ggtitle("Distribution of attributes used for filtering technically failed chips")
  # allcategplot2

## ----calcmahal-----------------------------------------------------------
expr_con_mahal <- expr_con[, noncorr]
expr_con_mahal <- scale(expr_con_mahal)
row.names(expr_con_mahal) <- expr_con$new_ID
hh(expr_con_mahal)
head(expr_con_mahal)
table(showNA(expr_con_mahal))
dim(expr_con_mahal)
set.seed(1902) # manchmal scheint mdqc inconistent zu sein, kann mich aber taeuschen. starte aber mit defnierten zufallsyahlen


if(robustmethod_for_mahal== "from_paramfile") robustmethod <- getParam2("robustmethod_for_mahal", myparam = param) else robustmethod = robustmethod_for_mahal
robustmethod


myqc <- mdqc::mdqc(expr_con_mahal[, ], method = "nogroups", robust = robustmethod)
# str(myqc)

# mdqc_vals <- myqc$mdqcValues[[1]]
# head(mdqc_vals)
# summary(myqc)
mahalnum <- (myqc$mdqcValues[[1]])
mahalnum <- log(mahalnum)
mahal <- data.frame(mahal = mahalnum,
                    new_ID = names(mahalnum))
ht(mahal)

if(filter2ind_atypischIlmnKontroll== "from_paramfile") mahalparam <- as.numeric(getParam2("filter2ind_atypischIlmnKontroll", myparam = param)) else mahalparam = filter2ind_atypischIlmnKontroll
mahalparam

mahalparam
badmahalcutoff <- mahalparam*IQR(mahalnum) + median(mahalnum)
badmahalcutoff

## ----visi----------------------------------------------------------------
par(mfrow = c(1, 1))
plotbereich <- c(1.2 * min(mahalnum, badmahalcutoff), 1.2 * max(c(mahalnum, badmahalcutoff), na.rm = T))
hist(mahalnum, breaks = 100, main = "Distribution Mahalanobis-Distance QC-parameter", xlim = plotbereich)
abline(v = median(mahalnum), col = "green")
abline(v = badmahalcutoff, col = "orange")
legend("topright", legend = c("green = median", paste0("orange = cutoff (", mahalparam, "x IQR+median)")))
# boxplot(mahalnum, main = "Verteilung Mahalanobis-Distance QC parameter", ylim = plotbereich)
# abline(h = median(mahalnum), col = "green")
# abline(h = badmahalcutoff, col="orange")
# legend("topright", legend = c("gruen = median", paste0("orange = ", mahalparam, "x IQR+median")))

## ----badind--------------------------------------------------------------
bad_mahal_inds <- mahal[mahal$mahal > badmahalcutoff, "new_ID"]
# venn2(bad_mahal_inds2, bad_mahal_inds)
length(bad_mahal_inds)
length(bad_mahal_inds) / length(pData(total_nobkgd_eset_ql)$sampleID)
bad_mahal_inds

# anfuegen relevanter attribute
recaclulatedvars <- c("ercc_1", "ercc_2", "ercc_3",
                      "ercc_4", "ercc_5", "hybrid_low",
                      "hybrid_med", "hybrid_high", "string_pm",
                      "string_mm", "biotin", "housekeeping",
                      "labeling", "ercc", "negative",
                      "samplemean", "p05", "p25",
                      "p50", "p75", "p95")
noncorr2add <- setdiff(c(paste0(noncorr[noncorr %in% recaclulatedvars], "_s06"),
                         noncorr[noncorr %nin% recaclulatedvars]),
                       names(sample_overview_l6))
noncorr2add
stopifnot(length(noncorr2add) == 0)

# for(i in noncorr2add) {
#   sample_overview_l6[i] = expr_con[match_hk(sample_overview_l6$new_ID,expr_con$new_ID), i]
# }
# setnames(sample_overview_l6, noncorr2add, paste0(noncorr2add, "_s06"))
sample_overview_l6$mahal <- mahal[match_hk(sample_overview_l6$new_ID, mahal$new_ID), "mahal"]
ht(sample_overview_l6, 1)
bad_mahal_indsdet <- sample_overview_l6[sample_overview_l6$new_ID %in% bad_mahal_inds, c("new_ID", "old_ID")]
bad_mahal_indsdet
ht(myplot, 1)
myplot$techn_problems <- myplot$new_ID %in% bad_mahal_inds
plotdata <- myplot[myplot$variable %nin% c('Detected.Genes..0.01.', 'Detected.Genes..0.05.'), ]
plotdata <- plotdata[order(plotdata$techn_problems), ]
allcategplot3 <- ggplot2::ggplot(plotdata, ggplot2::aes(variable, value, fill = variable) ) +
  ggplot2::geom_boxplot() +
  ggplot2::geom_point(ggplot2::aes(col = techn_problems, shape = techn_problems, size = techn_problems)) +
  ggplot2::theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.2, size = 12), axis.text.y = element_text(size = 12)) +
  ggplot2::guides(fill = F) +
  ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(10)) +
  ggplot2::scale_size_manual(values=c(2,4)) + ggplot2::ggtitle("Distribution of attributes used for filtering technically failed chips")
plot(allcategplot3)

# wie verteilen sich die schlechten ueber die experimente
table(sample_overview_l6[sample_overview_l6$new_ID %in% bad_mahal_inds, "fileset_id"])
table(as.character(sample_overview_l6[sample_overview_l6$new_ID %in% bad_mahal_inds, "Sentrix.Barcode"])) # changed 21.8.18 to fix case where no bad_mahal_inds exists


## ----save6---------------------------------------------------------------


sample_overview_l6[sample_overview_l6$new_ID %in% bad_mahal_inds,"in_study"] = F
xtabs_hk(~sample_overview_l6$in_study + sample_overview_l6$reason4exclusion)

# relative mahal berechnen start

sample_overview_l6 = data.table::data.table(data.frame(sample_overview_l6))

sample_overview_l6[, mymedian := median(mahal, na.rm = T)]
sample_overview_l6[, iqr:= IQR(mahal, na.rm = T)]
sample_overview_l6[, mahal_relIQR := (mahal - mymedian)/iqr]
sample_overview_l6[, suppressWarnings(range(mahal_relIQR, na.rm = T)), by = in_study] # changed 21.8.18 to fix case where no bad_mahal_inds exists
mahalparam
sample_overview_l6[in_study == T & (mahal_relIQR > mahalparam)]
sample_overview_l6[in_study == T, stopifnot(max(mahal_relIQR) <= mahalparam)]
sample_overview_l6[,mymedian := NULL]
sample_overview_l6[,iqr:= NULL]
setDF(sample_overview_l6)

# relative mahal berechnen ende

message("In total, ", ifelse(length(bad_mahal_inds)==0, 0, bad_mahal_inds), " samples were removed due to atypical combination of control probe levels, ", sum(sample_overview_l6$in_study), " samples are still valid...")

sample_overview_l6[sample_overview_l6$new_ID %in% bad_mahal_inds,
                   "reason4exclusion"] = "s06 extreme combination of expression control features"
mytable(sample_overview_l6$reason4exclusion)
mytable(sample_overview_l6$in_study)



## ----speichern---------------------------------------------------------





ht(sample_overview_l6, 1)
dim(sample_overview_l6)
names(sample_overview_l6)

ht12object$chipsamples = sample_overview_l6

# status fuer docku
plottts <- setdiff(ls()[grep("plot", ls())],
                   c( "barplothk2", 'mulitplot','myplot','norm_plot', "plotRanges"))
plottts


dim_all_nobkgd_expressed_eset_ql <- dim(all_nobkgd_expressed_eset_ql)
dim_all_nobkgd_expressed_eset_ql
fordoku = c(plottts,
                      "dim_all_nobkgd_expressed_eset_ql",
                      "bad_mahal_inds",
                      "sample_overview_l6",
                      "sample_overview_l6save",
                      "mahalnum",
                      "badmahalcutoff",
                      "mahalparam",
                      "phagenprobe",
                      "noncorr", "robustmethod")


stopifnot(sum(duplicated(fordoku))==0)


ht12object$dokuobjects_filterTechnicallyFailed = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_filterTechnicallyFailed) = fordoku

# recalculated control parameter
dim(expr_con)
ht12object$recalculated_illu_qc_parameter = expr_con

ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))

ht12object

}
