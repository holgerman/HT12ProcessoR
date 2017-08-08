makeCategorPseudoSNPDaten <- function(eset_indinfo, esetspalte = "hybridisierungchipserialnumber") {
    eset_indinfo2 = copy(eset_indinfo)
	eset_indinfo2 = data.table::data.table(data.frame(eset_indinfo2))
    spaltennamen = eset_indinfo2$sampleID
    sentrix_indinfo = as.matrix(data.table(sentrix = eset_indinfo2[, as.numeric(factor(get(esetspalte)))]))
    sentrix_indinfo = data.table(t(sentrix_indinfo))
    hh(sentrix_indinfo)
    setnames(sentrix_indinfo, spaltennamen)
    sentrix_indinfo = cbind(data.table(snpid = esetspalte), sentrix_indinfo)
    numcateg = eset_indinfo2[, length(unique(get(esetspalte)))]
    

    sentrix_indinfo = sentrix_indinfo[, `:=`(snpid, NULL)]
    sentrix_indinfo = as.matrix(sentrix_indinfo)
    row.names(sentrix_indinfo) = esetspalte
    hh(sentrix_indinfo)
    sliced_snps = MatrixEQTL::SlicedData$new()
    sliced_snps = sliced_snps$CreateFromMatrix(sentrix_indinfo)
    
    
    res = c()
    res$numcateg = numcateg
    res$sliced_snps = sliced_snps
    res$varname = esetspalte
    res
}
