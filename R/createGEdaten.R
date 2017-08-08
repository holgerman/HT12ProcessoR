createGEdaten <- function(ge_vor_combat) {
    ge_vor_combat = data.table(geneid = rownames(ge_vor_combat), ge_vor_combat)
    hh(ge_vor_combat, 9)
    zeillennamen = unlist(ge_vor_combat[, "geneid", with = F])
    ge_vor_combat = ge_vor_combat[, `:=`(geneid, NULL)]
    hh(ge_vor_combat)
    ge_vor_combat = as.matrix(ge_vor_combat)
    row.names(ge_vor_combat) = zeillennamen
    
    sliced_genes = MatrixEQTL::SlicedData$new()
    sliced_genes = sliced_genes$CreateFromMatrix(ge_vor_combat)
    sliced_genes
}
