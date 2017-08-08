#' @title compareEsetExpressions
#' @description FUNCTION_DESCRIPTION
#' @param all_nobkgd_eset_ql_combat PARAM_DESCRIPTION
#' @param adjusted_eset PARAM_DESCRIPTION
#' @param firstmessage PARAM_DESCRIPTION, Default: NULL
#' @param showPlots PARAM_DESCRIPTION, Default: T
#' @return OUTPUT_DESCRIPTION

#' @importFrom Biobase exprs
#' @importFrom Biobase exprs
#' @import data.table
#' @export

compareEsetExpressions <- function(all_nobkgd_eset_ql_combat, adjusted_eset, firstmessage = NULL, showPlots=T, compmethod = "spearman") {
    
  # 22/8 # 30/8 generalisisert 5/7/14 ausgabe erweitert
  eset1_names = deparse(substitute(all_nobkgd_eset_ql_combat))
  eset2_names = deparse(substitute(adjusted_eset))
  if (is.null(firstmessage)) message("comparing two expression sets ", eset1_names, " and ", eset2_names) else message(firstmessage)
    eset1name = deparse(substitute(all_nobkgd_eset_ql_combat))
    eset2name = deparse(substitute(adjusted_eset))
    
    probes1 = featureNames(all_nobkgd_eset_ql_combat)
    probes2 = featureNames(adjusted_eset)
    probesBoth = intersect(probes1, probes2)
    
    
    message("Found ", sum(probes1 %in% probes2), " of ", length(probes1), " probes of 1st eset in 2nd eset...\n  ")
    message("Found ", sum(probes2 %in% probes1), " of ", length(probes2), " probes of 2nd eset in 1st eset...\n  ")
    message("Working with intersect of ", length(probesBoth), " probes\n  ")
    
    inds1 = row.names(pData(all_nobkgd_eset_ql_combat))
    inds2 = row.names(pData(adjusted_eset))
    indsBoth = intersect(inds1, inds2)
    
    
    message("Found ", sum(inds1 %in% inds2), " of ", length(inds1), " inds of 1st eset in 2nd eset...\n  ")
    message("Found ", sum(inds2 %in% inds1), " of ", length(inds2), " inds of 2nd eset in 1st eset...\n  ")
    message("Working with intersect of ", length(indsBoth), " samples...\n  ")
    
    
    
    all_nobkgd_eset_ql_combat = all_nobkgd_eset_ql_combat[probesBoth, indsBoth]
    adjusted_eset = adjusted_eset[probesBoth, indsBoth]
    allcor_t <- c()
    time1 <- Sys.time()
    counter <- 1
    
    for (i in row.names(Biobase::exprs(all_nobkgd_eset_ql_combat))) {
        allcor_t[counter] <- cor(Biobase::exprs(all_nobkgd_eset_ql_combat)[i, ], Biobase::exprs(adjusted_eset)[i, ], method = compmethod)
        counter <- counter + 1
    }
    names(allcor_t) = rownames(Biobase::exprs(all_nobkgd_eset_ql_combat))
    head(allcor_t)
    
    
    
    ### plotten
    histgroupfak = round(length(allcor_t)/20, -1)
    histgroupfak = ifelse(histgroupfak == 0, "Sturges", histgroupfak)
    message("caclulated hist_groupfak probes=", histgroupfak, "\n\n")
    
    if(showPlots==T) {
    par(mfrow = c(2, 2))
    mytitle = paste0(compmethod, " correlation expression probes")
    
    hist(allcor_t, breaks = histgroupfak, main = mytitle, col = "lightblue")
    mtext(paste(eset1name, "\n", eset2name), cex = 0.7)
    
    boxplot(allcor_t, main = mytitle, col = "lightblue")
    mtext(paste(eset1name, "\n", eset2name), cex = 0.7)
    
    
    
    ## plotten vor adjustierung hinter adjustierung der correlation bzgl. der Individuen
    mytitle = paste0(compmethod, " correlation samples")
    
    allcor <- c()
    counter <- 1
    
    for (i in colnames(Biobase::exprs(all_nobkgd_eset_ql_combat))) {
        allcor[counter] <- cor(Biobase::exprs(all_nobkgd_eset_ql_combat)[, i], Biobase::exprs(adjusted_eset)[, i], method =compmethod)
        counter <- counter + 1
    }
    
    names(allcor) = colnames(Biobase::exprs(all_nobkgd_eset_ql_combat))
    head(allcor)
    
    
    ## plotten
    histgroupfak = max(round(length(allcor)/20, -1), 10)
    histgroupfak = ifelse(histgroupfak == 0, "Sturges", histgroupfak)
    
    hist(allcor, breaks = histgroupfak, main = mytitle, col = "orange")
    mtext(paste(eset1name, "\n", eset2name), cex = 0.7)
    
    boxplot(allcor, main = mytitle, col = "orange")
    mtext(paste(eset1name, "\n", eset2name), cex = 0.7)
    time1 <- Sys.time() - time1
    # print(time1)
    }
    res = list(cor_ind = allcor, cor_transkript = allcor_t)
    return(res)
    
}
