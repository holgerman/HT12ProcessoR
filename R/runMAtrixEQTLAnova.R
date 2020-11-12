runMAtrixEQTLAnova <- function(eset, output_file_name, myesetspalte, mycovarspalte = NULL, round4ANOVAcheck) {
    # eset = total_nobkgd_eset_ql; output_file_name = 'obj/s08_2_anova_sentrix_beforeCOMBAT'; myesetspalte = 'Sentrix.Barcode'; mycovarspalte = 'subgroup' eset =
    # total_nobkgd_eset_ql; output_file_name = 'obj/s08_2_anova_bigbatch_beforeCOMBAT'; myesetspalte = 'bigbatch' ; mycovarspalte = 'subgroup'
  # eset = total_nobkgd_eset_ql; output_file_name = NULL; myesetspalte = "Sentrix.Barcode";  mycovarspalte = "subgroup"

    ge_vor_combat = exprs(eset)  #[1:200, 1:100]
    hh(ge_vor_combat, 9)


    eset_indinfo2 = pData(eset)  #[1:100,]
    eset_indinfo = copy(eset_indinfo2)
    #data.table::setDT(eset_indinfo)
	eset_indinfo = data.table::data.table(data.frame(eset_indinfo))


    ht(eset_indinfo)

    stopifnot(identical(eset_indinfo$sampleID, colnames(ge_vor_combat)))

    sentrixdaten = makeCategorPseudoSNPDaten(eset_indinfo, esetspalte = myesetspalte)

    sentrixdaten$varname
    sentrixdaten$numcateg
    sentrixdaten$sliced_snps

    sliced_genes = createGEdaten(ge_vor_combat)
    sliced_genes$ResliceCombined(sliceSize = 5050)


    options(MatrixEQTL.ANOVA.categories = sentrixdaten$numcateg)


    if (length(mycovarspalte) > 0) {

        covarclass = class(eset_indinfo[, get(mycovarspalte)])
        message("covar is of class ", covarclass)
        if (covarclass %in% c("integer", "numeric"))
            covardaten = makeMetrischPseudoSNPDaten(pdata_eset = eset_indinfo, esetspalte = mycovarspalte)
        if (covarclass %in% c("character", "factor"))
            covardaten = makeCategorCovarDaten(eset_indinfo, esetspalte = mycovarspalte)
        if (covarclass %nin% c("character", "factor", "integer", "numeric"))
            stop("covar class must be integer, numeric, character or factor but is", covarclass)


        covardaten$varname
        covardaten$numcateg
        covardaten$sliced_snps

        sentrix_vor_combat = MatrixEQTL::Matrix_eQTL_engine(snps = sentrixdaten$sliced_snps, gene = sliced_genes, cvrt = covardaten$sliced_snps, output_file_name = output_file_name,
            pvOutputThreshold = 1, useModel = MatrixEQTL::modelANOVA, noFDRsaveMemory = F)

        # plot(sentrix_vor_combat)
        sentrix_vor_combat_pvals = sentrix_vor_combat$all$eqtls
        ht(sentrix_vor_combat_pvals, 5)

        names(sentrix_vor_combat_pvals$pvalue) = sentrix_vor_combat_pvals$gene

        df_vgl1 = data.frame(gx = ge_vor_combat[sentrix_vor_combat_pvals$gene[4], ], eset_indinfo[, mycovarspalte, with = F], confounder = eset_indinfo[, factor(get(myesetspalte))])
        vgl1_zeile = paste0("anova(lm( gx ~  ", paste(mycovarspalte, collapse = " + "), " + confounder, data=df_vgl1))")
        vgl1 = eval(parse(text = vgl1_zeile))
        vgl1
        check1 = identical(round(vgl1["confounder", "F value"], 5), round(sentrix_vor_combat_pvals$statistic[4], 5))
        # stopifnot(check1 & (is.na(vgl1$`F value`[1])==F))

        df_vgl2 = data.frame(gx = ge_vor_combat[sentrix_vor_combat_pvals$gene[nrow(sentrix_vor_combat_pvals)], ], eset_indinfo[, mycovarspalte, with = F], confounder = eset_indinfo[,
            factor(get(myesetspalte))])
        vgl2_zeile = paste0("anova(lm( gx ~  ", paste(mycovarspalte, collapse = " + "), " + confounder, data=df_vgl2))")
        vgl2 = eval(parse(text = vgl2_zeile))
        vgl2
        check2 = identical(round(vgl2["confounder", "F value"], 5), round(sentrix_vor_combat_pvals$statistic[nrow(sentrix_vor_combat_pvals)], 5))

        # stopifnot(check2 & (is.na(vgl2$`F value`[2])==F))

        if (check1 & (is.na(vgl1$`F value`[1]) == F) & check2 & (is.na(vgl2$`F value`[2]) == F)) {
            message("For two examples, R ANOVA and MatrixEQTL ANOVA incl. covariate are identical. Great.")
        } else {

            message("For two examples, R ANOVA and MatrixEQTL ANOVA incl. covariate were not identical. Rerunning with classical anova")

            allgenes = sentrix_vor_combat_pvals$gene
            stopifnot(identical(colnames(ge_vor_combat), eset_indinfo$sampleID))
            oldrownames = rownames(sentrix_vor_combat_pvals)
            rownames(sentrix_vor_combat_pvals) = allgenes

            pb <- txtProgressBar(min = 0, max = length(allgenes), style = 3)

            for (mygenenum in seq(along = allgenes)) {
                mygene = allgenes[mygenenum]
                # df_vgl1 = data.frame(gx = ge_vor_combat[mygene,], eset_indinfo[,mycovarspalte, with = F], confounder = eset_indinfo[,factor(get(myesetspalte))])
                resi = anova(lm(ge_vor_combat[mygene, ] ~ eset_indinfo[, get(mycovarspalte)] + eset_indinfo[, factor(get(myesetspalte))]))
                # str(resi)
                sentrix_vor_combat_pvals[mygene, c("statistic", "pvalue", "FDR")] = c(resi$`F value`[2], resi$`Pr(>F)`[2], NA)
                setTxtProgressBar(pb, mygenenum)
            }

            message("Done rerun classcial ANOVA")

        }


    } else {


        sentrix_vor_combat = MatrixEQTL::Matrix_eQTL_engine(snps = sentrixdaten$sliced_snps, gene = sliced_genes, output_file_name = output_file_name, pvOutputThreshold = 1, useModel = MatrixEQTL::modelANOVA,
            noFDRsaveMemory = F)

        # plot(sentrix_vor_combat)
        sentrix_vor_combat_pvals = sentrix_vor_combat$all$eqtls
        ht(sentrix_vor_combat_pvals, 5)

        names(sentrix_vor_combat_pvals$pvalue) = sentrix_vor_combat_pvals$gene

        vgl1 = anova(lm(ge_vor_combat[sentrix_vor_combat_pvals$gene[4], ] ~ factor(unlist(eset_indinfo[, myesetspalte, with = F]))))
        vgl1
        anova_f1 = vgl1$`F value`[1]
        matrixeqtl_f1 = sentrix_vor_combat_pvals$statistic[4]
        check1 = identical(round(anova_f1, round4ANOVAcheck), round(matrixeqtl_f1, round4ANOVAcheck))
        message('Checking Example 1 F value Anova :', anova_f1)
        message('Checking Example 1 F value MatrixEQTLS :', matrixeqtl_f1)
        stopifnot(check1 & (is.na(vgl1$`F value`[1]) == F))

        vgl2 = anova(lm(ge_vor_combat[sentrix_vor_combat_pvals$gene[nrow(sentrix_vor_combat_pvals)], ] ~ factor(unlist(eset_indinfo[, myesetspalte, with = F]))))
        vgl2

        anova_f2 = vgl2$`F value`[1]
        matrixeqtl_f2 = sentrix_vor_combat_pvals$statistic[nrow(sentrix_vor_combat_pvals)]
        check2 = identical(round(anova_f2, round4ANOVAcheck), round(matrixeqtl_f2, round4ANOVAcheck))
        message('Checking Example 2 F value Anova :', anova_f2)
        message('Checking Example 2 F value MatrixEQTLS :', matrixeqtl_f2)

        check2 = identical(round(anova_f2, round4ANOVAcheck), round(matrixeqtl_f2, round4ANOVAcheck))
        stopifnot(check2 & (is.na(vgl2$`F value`[1]) == F))
        if (check1 & check2)
            message("For two examples, R ANOVA and MatrixEQTL in  ANOVA without Covariate are sufficiently identical (rounding both F values by ",round4ANOVAcheck,"). Great.")
    }
    sentrix_vor_combat_pval = sentrix_vor_combat_pvals$pvalue
    names(sentrix_vor_combat_pval) = sentrix_vor_combat_pvals$gene
    sentrix_vor_combat_pval
}
