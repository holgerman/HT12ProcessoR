makeCategorCovarDaten <- function(eset_indinfo, esetspalte = "subgroup") {
 
    eset_indinfo2 = data.table::copy(eset_indinfo)
    eset_indinfo2 = data.table::data.table(data.frame(eset_indinfo2))
    spaltennamen = eset_indinfo2$sampleID
    data.table::setnames(eset_indinfo2, esetspalte, "covar")
    
    covar = model.matrix(~covar, eset_indinfo2)
    covar = covar[, grep("covar", colnames(covar), value = T), drop=F]
    covar_indinfo = covar
    covar_indinfo = data.table::data.table(t(covar_indinfo))
    hh(covar_indinfo)
    data.table::setnames(covar_indinfo, spaltennamen)
    covar_indinfo = cbind(data.table::data.table(snpid = colnames(covar)), covar_indinfo)
    
    
    covar_indinfo = covar_indinfo[, `:=`(snpid, NULL)]
    covar_indinfo = as.matrix(covar_indinfo)
    row.names(covar_indinfo) = colnames(covar)
    hh(covar_indinfo)
    sliced_snps = MatrixEQTL::SlicedData$new()
    sliced_snps =sliced_snps$CreateFromMatrix(covar_indinfo)
    
    
    res = c()
    res$numcateg = NULL
    res$sliced_snps = sliced_snps
    res$varname = esetspalte
    res
}
