makeMetrischPseudoSNPDaten <- function(pdata_eset, esetspalte = "Alter") {
    # pdata_eset = pData(eset_preproc)
    
    pdata_eset2 = copy(pdata_eset)
    #setDT(pdata_eset2)
	pdata_eset2 = data.table::data.table(data.frame(pdata_eset2))
        
    spaltennamen = pdata_eset2$sampleID
    
    metrische_var_indinfo = as.matrix(data.table(metrische_var = pdata_eset2[, esetspalte, with = F]))
    metrische_var_indinfo = data.table(t(metrische_var_indinfo))
    hh(metrische_var_indinfo)
    setnames(metrische_var_indinfo, spaltennamen)
    
   
    metrische_var_indinfo = as.matrix(metrische_var_indinfo)
    row.names(metrische_var_indinfo) = esetspalte
    hh(metrische_var_indinfo)
    sliced_snps = MatrixEQTL::SlicedData$new()
    sliced_snps = sliced_snps$CreateFromMatrix(metrische_var_indinfo)
    
    
    res = c()
    res$sliced_snps = sliced_snps
    res$varname = esetspalte
    res
}
