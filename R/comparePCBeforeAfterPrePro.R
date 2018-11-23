#' @title Visualize PCA before and after preprocessing expression levels
#'
#' @description PCA on expression data is done using good samples surviving preprocessing. Following PCAs are calculated and shown as figures: PCA of all transcripts before preprocessing, including additional scaling, PCA of gene expression after preprocessing using only expression probes that are expressed, not overinflated regarding batch affects and beeing specific to the human genome when considering a remapping approach according to Barbosa-Moralis et al. (2010)
#'

#' @param ht12object A list object of class HT12prepro created with function visualizePreprocessing
#' @param paramfile Path to the file specifying parameters
#' @return A list object of class HT12prepro including updated slots, e.g   the slot with the history of the commands named `$history`` is updated. Plots are also stored in slot `$`
#' @import data.table
#' @export

## debug
# paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
# ht12object =  prepro_ht12
# showPlots = T

comparePCBeforeAfterPrePro = function(ht12object,paramfile = NULL,showPlots=T) {

### Do you want to automatically convert strings to factor variables in a data.frame? WARNING!!! This makes your code less portable/reproducible.
options(stringsAsFactors=FALSE)

  myparameters = match.call()
  showVennplots = F

  # status checken
  historie =  ht12object$history$calls
  if(any(grepl("visualizePreprocessing", historie))==F) stop("Function 'visualizePreprocessing()' has to be run before!")

  #laden parameter

  if(is.null(paramfile)==F) param <- data.frame(data.table::fread(paramfile))




head(ilmnAnnot014allgInfos)

## ----functions-----------------------------------------------------------

plotte3D <- function (tocolor, sample_overview_l10, matrix3spalt, mylabels, mysize= 1.5) {

  farbe = factor(sample_overview_l10[match_hk(rownames(matrix3spalt), sample_overview_l10$new_ID), tocolor])
  farbentopf = rainbow(length(unique(farbe)))
  farbe = as.character(factor(farbe, labels = farbentopf) )
  try(threejs::scatterplot3js(matrix3spalt, color=farbe, labels=mylabels, size=mysize ))
}


sample_overview_l10 <- ht12object$chipsamples
mytable(sample_overview_l10$in_study)
sample_overview_l10instudy <- sample_overview_l10[ sample_overview_l10$in_study,]
dim(sample_overview_l10)
dim(sample_overview_l10instudy)
table(table(sample_overview_l10instudy$new_ID))
if(length(table(table(sample_overview_l10instudy$new_ID))) != 1)
  stop("IDs (column new_ID) must be unique....stopping...")

# laden annotation probes
genesdetail <- ht12object$genesdetail
mytable(genesdetail$expressed)
ht(genesdetail, 2)
mytable(stringr::str_sub(genesdetail$ilmn, 1, 4))


#laden expressionsets nach combat
total_nobkgd_eset_ql_combat = ht12object$total_nobkgd_eset_ql_combat
total_nobkgd_eset_ql_combat

#laden expressionsets VOR combat
total_nobkgd_eset= ht12object$total_nobkgd_eset
total_nobkgd_eset

## ----good----------------------------------------------------------------
goodind <- sample_overview_l10instudy$new_ID
total_nobkgd_eset
total_nobkgd_eset <- total_nobkgd_eset[, goodind]
total_nobkgd_eset
total_nobkgd_eset_ql_combat

## ----pcinitial-----------------------------------------------------------
data.pc__allinitial <- Biobase::exprs(total_nobkgd_eset[, goodind])
data.pc__allinitial <- t(data.pc__allinitial)
hh(data.pc__allinitial)
dim(data.pc__allinitial)
pc.data_allinitial <- prcomp(data.pc__allinitial) # this is not scaled!
#screeplot(pc.data_allinitial, main = "Eigenvalues")
# str(pc.data_allinitial)
form <- as.numeric(factor(sample_overview_l10[match_hk(rownames(data.pc__allinitial),
                                                       sample_overview_l10$new_ID), "strangebatch"]))
farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc__allinitial),
                                             sample_overview_l10$new_ID), "subgroup"])
n_good__allinitial <- dim(data.pc__allinitial)[2]
(cumsum((pc.data_allinitial$sdev)^2) / sum(pc.data_allinitial$sdev^2))[1:10]
# if(showPlots) pairs(pc.data_allinitial$x[, 1:7], col = farbe,
      # pch = form,
      # main = paste0("all initial values log2 only ",
                    # n_good__allinitial, " transcripts"))
texte = paste(c("Colored according subgroup ",
                paste(levels(farbe), collapse = ", "),
                " shown in ",
                paste(palette()[1:length(levels(farbe))],
                      collapse = ", "),
                ", respectively. form = specialbatch"),
              collapse = "")
# if(showPlots) mtext(texte, line = 1)
beschriftung13a = data.frame("ID1:",
	new_ID = rownames(pc.data_allinitial$x[,1:3]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,1:3]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,1:3]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,1:3]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,1:3]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data_allinitial$x[,1:3], mylabels = beschriftung13a)
beschriftung46 = data.frame("ID1:",
	new_ID = rownames(pc.data_allinitial$x[,4:6]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,4:6]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,4:6]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,4:6]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data_allinitial$x[,4:6]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data_allinitial$x[,4:6], mylabels = beschriftung46)

## ----pcinitial_scaled----------------------------------------------------
data.pc__allinitial_scaled <- Biobase::exprs(total_nobkgd_eset[, goodind])
data.pc__allinitial_scaled <- t(data.pc__allinitial_scaled)
hh(data.pc__allinitial_scaled)
dim(data.pc__allinitial_scaled)
pc.data_allinitial_scaled <- prcomp(data.pc__allinitial_scaled,   center = T, scale. = T)
#screeplot(pc.data_allinitial_scaled, main = "Eigenvalues")
# str(pc.data_allinitial_scaled)
form <- as.numeric(factor(sample_overview_l10[match_hk(rownames(data.pc__allinitial_scaled), sample_overview_l10$new_ID), "strangebatch"]))
farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc__allinitial_scaled), sample_overview_l10$new_ID), "subgroup"])
n_good__allinitial_scaled <- dim(data.pc__allinitial_scaled)[2]
(cumsum((pc.data_allinitial_scaled$sdev)^2) / sum(pc.data_allinitial_scaled$sdev^2))[1:10]
if(showPlots) pairs(pc.data_allinitial_scaled$x[, 1:7], col = farbe,
      pch = form,
      main = paste0("all initial_scaled values log2 only ",
      n_good__allinitial_scaled,
      " transcripts"))
texte <- paste(c("Colored according subgroup ",
                 paste(levels(farbe),
                       collapse = ", "),
                 " shown in ",
                 paste(palette()[1:length(levels(farbe))],
                       collapse = ", "),
                 ", respectively."),
               collapse = "")
if(showPlots) mtext(texte,  outer = T, line = 1.5)


farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc__allinitial_scaled), sample_overview_l10$new_ID), "processingbatch"])

if(showPlots) pairs(pc.data_allinitial_scaled$x[, 1:7], col = farbe,
                    pch = form,
                    main = paste0("all initial_scaled values log2 only ",
                                  n_good__allinitial_scaled,
                                  " transcripts"))
texte <- paste("Colored according processingbatch ")
if(showPlots) mtext(texte,  outer = T, line = 1.5)






beschriftung13 = data.frame("ID1:",
	new_ID = rownames(pc.data_allinitial_scaled$x[,1:3]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,1:3]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,1:3]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,1:3]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,1:3]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data_allinitial_scaled$x[,1:3], mylabels = beschriftung13)
beschriftung46 = data.frame("ID1:",
	new_ID = rownames(pc.data_allinitial_scaled$x[,4:6]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,4:6]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,4:6]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,4:6]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data_allinitial_scaled$x[,4:6]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data_allinitial_scaled$x[,4:6], mylabels = beschriftung46)


## ----pca-----------------------------------------------------------------
all_goodexpressedprobe <- grep("goodexpressed", names(genesdetail), value = T)
all_goodexpressedprobe
genesdetail$goodexpressedprobe <- apply(genesdetail[, all_goodexpressedprobe, drop = F], 1, function(x) all(x))
mytable(genesdetail$goodexpressedprobe)
goodprobes <- genesdetail[genesdetail$goodexpressedprobe, "nuid"]
# str(goodprobes)
data.pc <- Biobase::exprs(total_nobkgd_eset_ql_combat[goodprobes, goodind])
data.pc <- t(data.pc)
hh(data.pc)

pc.data <- prcomp(data.pc) # this is not scaled! I want it as high expresssed appears more relevant and quantile normalisztion is already done
# plot(pc.data, main = "Eigenvalues")
# str(pc.data)
form <- as.numeric(factor(sample_overview_l10[match_hk(rownames(data.pc), sample_overview_l10$new_ID), "strangebatch"]))
farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc), sample_overview_l10$new_ID), "subgroup"])
(cumsum((pc.data$sdev)^2) / sum(pc.data$sdev^2))[1:10]
#pairs(pc.data$x[, 1:6], col = farbe,  pch = form, main  = paste0("good expressed ", dim(data.pc)[2], " transcripts"))
texte = paste(c("Colored according subgroup ",
                paste(levels(farbe), collapse = ", "),
                " shown in ",
                paste(palette()[1:length(levels(farbe))], collapse = ", "),
                ", respectively."),
              collapse = "")
#mtext(texte,   outer = T, line = 1.5)
#pairs(pc.data$x[, 7:12], col = farbe,  pch = form, main  = paste0("good expressed ", dim(data.pc)[2], " transcripts"))
texte = paste(c("Colored according subgroup ",
                paste(levels(farbe), collapse = ", "),
                " shown in ",
                paste(palette()[1:length(levels(farbe))], collapse = ", "),
                ", respectively."),
              collapse = "")
#mtext(texte,  outer = T, line = 1.5)
beschriftung13 = data.frame("ID1:",
	new_ID = rownames(pc.data$x[,1:3]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data$x[,1:3]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data$x[,1:3]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data$x[,1:3]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data$x[,1:3]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data$x[,1:3], mylabels = beschriftung13)
beschriftung46 = data.frame("ID1:",
	new_ID = rownames(pc.data$x[,4:6]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.data$x[,4:6]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.data$x[,4:6]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.data$x[,4:6]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.data$x[,4:6]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.data$x[,4:6], mylabels = beschriftung46)


## ----pcalltranscripts----------------------------------------------------
data.pcall <- Biobase::exprs(total_nobkgd_eset_ql_combat[, goodind])
data.pcall <- t(data.pcall)
hh(data.pcall)
pc.dataall <- prcomp(data.pcall) # this is not scaled! I want it as high expresssed appears more relevant and quantile normalisztion is already done
#screeplot(pc.dataall, main = "Eigenvalues")
# str(pc.dataall)
form <- as.numeric(factor(sample_overview_l10[match_hk(rownames(data.pc), sample_overview_l10$new_ID), "strangebatch"]))
farbe <- factor(sample_overview_l10[ match_hk(rownames(data.pc), sample_overview_l10$new_ID), "subgroup"])
(cumsum((pc.dataall$sdev)^2) / sum(pc.dataall$sdev^2))[1:10]
#pairs(pc.dataall$x[,1:7], col = farbe,  pch = form, main  = paste0("all ", dim(data.pcall)[2], " transcripts"))
texte <- paste(c("Colored according subgroup ",
                 paste(levels(farbe), collapse = ", "),
                 " shown in ",
                 paste(palette()[1:length(levels(farbe))], collapse = ", "),
                 ", respectively."),
               collapse = "")
#mtext(texte, line = 1)
beschriftung13 = data.frame("ID1:",
	new_ID = rownames(pc.dataall$x[,1:3]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.dataall$x[,1:3]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.dataall$x[,1:3]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.dataall$x[,1:3]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.dataall$x[,1:3]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.dataall$x[,1:3], mylabels = beschriftung13)
beschriftung46 = data.frame("ID1:",
	new_ID = rownames(pc.dataall$x[,4:6]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.dataall$x[,4:6]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.dataall$x[,4:6]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.dataall$x[,4:6]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.dataall$x[,4:6]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.dataall$x[,4:6], mylabels = beschriftung46)

## ----pcagood-------------------------------------------------------------
genesdetail$badprobe_dunning <- ilmnAnnot014allgInfos[match_hk(genesdetail$nuid, ilmnAnnot014allgInfos$nuid), "badprobe_dunning"]
mytable(genesdetail$badprobe_dunning)
goodprobes_inkldunning <- genesdetail[genesdetail$goodexpressedprobe & (genesdetail$badprobe_dunning == F), "nuid"]
# str(goodprobes_inkldunning)
table(is.na(goodprobes_inkldunning))
data.pc_dunning <- Biobase::exprs(total_nobkgd_eset_ql_combat[ goodprobes_inkldunning, goodind])
data.pc_dunning <- t(data.pc_dunning)
hh(data.pc_dunning)
dim(data.pc_dunning)
pc.datadunning <- prcomp(data.pc_dunning ) # this is not scaled! I want it as high expresssed appears more relevant and quantile normalisztion is already done
#screeplot(pc.datadunning, main = "Eigenvalues")
# str(pc.datadunning)
form <- as.numeric(factor(sample_overview_l10[match_hk(rownames(data.pc_dunning), sample_overview_l10$new_ID), "strangebatch"]))
farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc_dunning), sample_overview_l10$new_ID), "subgroup"])
n_good_dunning <- dim(data.pc_dunning)[2]
(cumsum((pc.datadunning$sdev)^2) / sum(pc.datadunning$sdev^2))[1:10]
if(showPlots) pairs(pc.datadunning$x[, 1:7],
      col = farbe,
      pch = form,
      main = paste0("good preproc & good mapping ", n_good_dunning, " transcripts"))
texte <- paste(c("Colored according subgroup ",
                 paste(levels(farbe), collapse = ", "),
                 " shown in ",
                 paste(palette()[1:length(levels(farbe))], collapse = ", "),
                 ", respectively."), collapse = "")
if(showPlots) mtext(texte,  outer = T, line = 1.5)


farbe <- factor(sample_overview_l10[match_hk(rownames(data.pc_dunning), sample_overview_l10$new_ID), "processingbatch"])
if(showPlots) pairs(pc.datadunning$x[, 1:7],
                    col = farbe,
                    pch = form,
                    main = paste0("good preproc & good mapping ", n_good_dunning, " transcripts"))
texte <- paste("Colored according processingbatch")
if(showPlots) mtext(texte,  outer = T, line = 1.5)










beschriftung13d = data.frame("ID1:",
	new_ID = rownames(pc.datadunning$x[,1:3]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,1:3]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,1:3]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,1:3]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,1:3]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
if(showPlots) plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.datadunning$x[,1:3], mylabels = beschriftung13d)



beschriftung46 = data.frame("ID1:",
	new_ID = rownames(pc.datadunning$x[,4:6]), "ID2:",
	old_ID = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,4:6]), sample_overview_l10$new_ID), "old_ID"],  "procbatch:",
	fileset_id = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,4:6]), sample_overview_l10$new_ID), "fileset_id"],"Sentrix:",
	Sentrix = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,4:6]), sample_overview_l10$new_ID), "Sentrix.Barcode"], "Numgens:",
	genes01 = sample_overview_l10[match_hk(rownames(pc.datadunning$x[,4:6]), sample_overview_l10$new_ID), "Detected.Genes..0.01."])
# plotte3D(tocolor = "Sentrix.Barcode", sample_overview_l10, matrix3spalt = pc.datadunning$x[,4:6], mylabels = beschriftung46)

## ----speichern-----------------------------------------------------------


# for doku
n_maxsave <- min(c(100, dim(pc.dataall$x)[2]))
pc__all_rotated_1to100 <- pc.dataall$x[, 1:n_maxsave]
pc__all_sd <- pc.dataall$sdev
pc__good_rotated_1to100 <- pc.data$x[, 1:n_maxsave]
pc__good_sd <- pc.data$sdev
pc__gooddunning_rotated_1to100 <- pc.datadunning$x[, 1:n_maxsave]
pc__gooddunning_sd <- pc.datadunning$sdev
pc__all_initial_rotated_1to100 <- pc.data_allinitial$x[, 1:n_maxsave]
pc__all_initial_sd <- pc.data_allinitial$sdev
pc__all_initial_scaled_rotated_1to100 <- pc.data_allinitial_scaled$x[, 1:n_maxsave]
pc__all_initial_scaled_sd <- pc.data_allinitial_scaled$sdev
pcs <- grep('^pc\\.', ls(), value = T)
pcs_doku <- grep('^pc__', ls(), value = T)
pcs_doku

fordoku =  c(pcs, pcs_doku, 'goodprobes', 'n_good_dunning')

stopifnot(sum(duplicated(fordoku))==0)


ht12object$dokuobjects_comparePCBeforeAfterPrePro = lapply(fordoku, function(x) get(x))


names(ht12object$dokuobjects_comparePCBeforeAfterPrePro) = fordoku


ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
ht12object$history
ht12object


}

