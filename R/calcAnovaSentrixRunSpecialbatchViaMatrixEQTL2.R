calcAnovaSentrixRunSpecialbatchViaMatrixEQTL2 <- function(total_nobkgd_eset_ql, total_nobkgd_eset_ql_combat, genesdetail, subgroups, anova_with_special_batch, strangebatch,
                                                          adjust4subgroup = F) {

  ##debug
  # adjust4subgroup = T


  orinames = names(genesdetail)
  par(mfrow = c(2, 2))

  n_sentrix = length(unique(pData(total_nobkgd_eset_ql)$Sentrix.Barcode))
  message("Found ", n_sentrix, " different sentrix IDs...")
  n_bigbatch = length(unique(pData(total_nobkgd_eset_ql)$bigbatch))
  if(n_bigbatch==1) message("No different batches present for the processingbatch") else message("Found ", n_bigbatch, " different processingbatches ...")

  n_strangebatch = length(unique(pData(total_nobkgd_eset_ql)$strangebatch))
  if(n_strangebatch==1) message("No different batches present assigned as specialbatches") else message("Found ", n_strangebatch, " different specialbatches (i.e. custom annotated in s02...Rmd) ...")


  if (adjust4subgroup == T) {

    message("... acounting for covariate `subgroup` in ANOVA ...")

    if (n_sentrix <= 1) {
      all_pval_before_sentrix = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql)))
      all_pval_after_sentrix = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql_combat)))

      names(all_pval_before_sentrix) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_sentrix) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_sentrix == 2) {


      ## berechnen ANOVA PREcombatvariante Sentrix.Barcode

      pData(total_nobkgd_eset_ql)$Sentrix.Barcode = factor(pData(total_nobkgd_eset_ql)$Sentrix.Barcode)
      all_pval_before_sentrix = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "Sentrix.Barcode", adjust4subgroup = T, mycovarspalte = "subgroup")
      # str(all_pval_before_sentrix)


      ## berechnen ANOVA POSTcombatvariante Sentrix.Barcode

      pData(total_nobkgd_eset_ql_combat)$Sentrix.Barcode = factor(pData(total_nobkgd_eset_ql_combat)$Sentrix.Barcode)

      all_pval_after_sentrix = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "Sentrix.Barcode", adjust4subgroup = T, mycovarspalte = "subgroup")



      # str(all_pval_after_sentrix)




    }
    if (n_sentrix > 2) {
      message("Running ANOVA in the form of high-efficient eQTL analysis via package MatrixEQTL ....\n")
      message("Calculating ANOIVA for Sentrix.Barcode before Combat...")
      all_pval_before_sentrix = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name = NULL, myesetspalte = "Sentrix.Barcode",
                                                   mycovarspalte = "subgroup")
      # str(all_pval_before_sentrix)
      message("Calculating ANOIVA for Sentrix.Barcode after Combat...")
      all_pval_after_sentrix = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = NULL, myesetspalte = "Sentrix.Barcode",
                                                  mycovarspalte = "subgroup")
      # str(all_pval_after_sentrix)
    }



    if (n_bigbatch == 1) {

      all_pval_before_bigbatch = rep(x = 1, times = length(all_pval_before_sentrix))
      all_pval_after_bigbatch = rep(x = 1, times = length(all_pval_after_sentrix))

      names(all_pval_before_bigbatch) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_bigbatch) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_bigbatch == 2) {
      ## berechnen ANOVA PREcombatvariante bigbatch i.e. run

      pData(total_nobkgd_eset_ql)$bigbatch = factor(pData(total_nobkgd_eset_ql)$bigbatch)
      all_pval_before_bigbatch = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "bigbatch", adjust4subgroup = T, mycovarspalte = "subgroup")
      # str(all_pval_before_bigbatch)


      ## berechnen ANOVA POSTcombatvariante bigbatch i.e. run

      pData(total_nobkgd_eset_ql_combat)$bigbatch = factor(pData(total_nobkgd_eset_ql_combat)$bigbatch)

      all_pval_after_bigbatch = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "bigbatch", adjust4subgroup = T, mycovarspalte = "subgroup")



      # str(all_pval_after_bigbatch)

    }
    if (n_bigbatch > 2) {
      message("Calculating ANOIVA for Processing Batch (i.e. fileset_id) before Combat...")
      all_pval_before_bigbatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name = NULL, myesetspalte = "bigbatch",
                                                    mycovarspalte = "subgroup")
      # str(all_pval_before_bigbatch)

      message("Calculating ANOIVA for Processing Batch (i.e. fileset_id) after Combat...")
      all_pval_after_bigbatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = NULL, myesetspalte = "bigbatch",
                                                   mycovarspalte = "subgroup")
      # str(all_pval_after_bigbatch)
    }

    if (n_strangebatch == 1) {
      all_pval_before_strangebatch = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql)))
      all_pval_after_strangebatch = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql_combat)))

      names(all_pval_before_strangebatch) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_strangebatch) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_strangebatch == 2) {
      ## berechnen ANOVA PREcombatvariante strangebatch

      pData(total_nobkgd_eset_ql)$strangebatch = factor(pData(total_nobkgd_eset_ql)$strangebatch)
      all_pval_before_strangebatch = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "strangebatch", adjust4subgroup = T, mycovarspalte = "subgroup")
      # str(all_pval_before_strangebatch)


      ## berechnen ANOVA POSTcombatvariante strangebatch

      pData(total_nobkgd_eset_ql_combat)$strangebatch = factor(pData(total_nobkgd_eset_ql_combat)$strangebatch)

      all_pval_after_strangebatch = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "strangebatch", adjust4subgroup = T, mycovarspalte = "subgroup")



      # str(all_pval_after_strangebatch)
    }
    if (n_strangebatch > 2) {
      message("Calculating ANOIVA for userdefined strangebatch before Combat...")
      all_pval_before_strangebatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name = "obj/s08_2_anova_strangebatch_beforeCOMBAT", myesetspalte = "strangebatch",
                                                        mycovarspalte = "subgroup")
      # str(all_pval_before_strangebatch)
      message("Calculating ANOIVA for userdefined strangebatch after Combat...")
      all_pval_after_strangebatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = "obj/s08_2_anova_strangebatch_afterCOMBAT", myesetspalte = "strangebatch",
                                                       mycovarspalte = "subgroup")
      # str(all_pval_after_strangebatch)
    }

  } else {
    message("... not adjusting ANOVA for subgroup ...")

    if (n_sentrix <= 1) {
      all_pval_before_sentrix = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql)))
      all_pval_after_sentrix = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql_combat)))

      names(all_pval_before_sentrix) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_sentrix) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_sentrix == 2) {
      pData(total_nobkgd_eset_ql)$Sentrix.Barcode = factor(pData(total_nobkgd_eset_ql)$Sentrix.Barcode)
      all_pval_before_sentrix = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "Sentrix.Barcode")
      # str(all_pval_before_sentrix)

      pData(total_nobkgd_eset_ql_combat)$Sentrix.Barcode = factor(pData(total_nobkgd_eset_ql_combat)$Sentrix.Barcode)
      all_pval_after_sentrix = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "Sentrix.Barcode")
      # str(all_pval_after_sentrix)

    }

    if (n_sentrix > 2) {
      message("Calculating ANOIVA for Sentrix.Barcode (n_sentrix > 2) before Combat...")
      all_pval_before_sentrix = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name = NULL, myesetspalte = "Sentrix.Barcode")
      # str(all_pval_before_sentrix)
      message("Calculating ANOIVA for Sentrix.Barcode (n_sentrix > 2) after Combat...")
      all_pval_after_sentrix = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = NULL, myesetspalte = "Sentrix.Barcode")
      # str(all_pval_after_sentrix)
    }


    if (n_bigbatch == 1) {
      all_pval_before_bigbatch = rep(x = 1, times = length(all_pval_before_sentrix))
      all_pval_after_bigbatch = rep(x = 1, times = length(all_pval_after_sentrix))

      names(all_pval_before_bigbatch) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_bigbatch) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_bigbatch == 2) {
      pData(total_nobkgd_eset_ql)$bigbatch = factor(pData(total_nobkgd_eset_ql)$bigbatch)
      all_pval_before_bigbatch = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "bigbatch")
      # str(all_pval_before_bigbatch)

      pData(total_nobkgd_eset_ql_combat)$bigbatch = factor(pData(total_nobkgd_eset_ql_combat)$bigbatch)
      all_pval_after_bigbatch = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "bigbatch")
      # str(all_pval_after_bigbatch)

    }
    if (n_bigbatch > 2) {
      message("Calculating ANOIVA for Processing Batch (i.e. fileset_id, n_fileset_id >2) before Combat...")
      all_pval_before_bigbatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name =NULL, myesetspalte = "bigbatch")
      # str(all_pval_before_bigbatch)

      message("Calculating ANOIVA for Processing Batch (i.e. fileset_id, n_fileset_id >2) after Combat...")
      all_pval_after_bigbatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = NULL, myesetspalte = "bigbatch")
      # str(all_pval_after_bigbatch)
    }

    if (n_strangebatch == 1) {
      all_pval_before_strangebatch = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql)))
      all_pval_after_strangebatch = rep(x = 1, times = nrow(exprs(total_nobkgd_eset_ql_combat)))

      names(all_pval_before_strangebatch) = rownames(exprs(total_nobkgd_eset_ql))
      names(all_pval_after_strangebatch) = rownames(exprs(total_nobkgd_eset_ql_combat))
    }
    if (n_strangebatch == 2) {
      pData(total_nobkgd_eset_ql)$strangebatch = factor(pData(total_nobkgd_eset_ql)$strangebatch)
      all_pval_before_strangebatch = calcAnova(eset = total_nobkgd_eset_ql, toadjust.hyb = "strangebatch")
      # str(all_pval_before_strangebatch)

      pData(total_nobkgd_eset_ql_combat)$strangebatch = factor(pData(total_nobkgd_eset_ql_combat)$strangebatch)
      all_pval_after_strangebatch = calcAnova(eset = total_nobkgd_eset_ql_combat, toadjust.hyb = "strangebatch")
      # str(all_pval_after_strangebatch)
    }
    if (n_strangebatch > 2) {

      message("Calculating ANOIVA for userdefined strangebatch  (n_strangebatch > 2) before Combat...")
      all_pval_before_strangebatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql, output_file_name = NULL, myesetspalte = "strangebatch")
      # str(all_pval_before_strangebatch)

      message("Calculating ANOIVA for userdefined strangebatch  (n_strangebatch > 2) after Combat...")
      all_pval_after_strangebatch = runMAtrixEQTLAnova(eset = total_nobkgd_eset_ql_combat, output_file_name = NULL, myesetspalte = "strangebatch")
      # str(all_pval_after_strangebatch)
    }
  }







  ## zusammenfassen ergebnisframe
  eigenschaften = data.frame(probeID = rownames(total_nobkgd_eset_ql))

  eigenschaften$all_pval_after_sentrix = all_pval_after_sentrix[match_hk(eigenschaften$probeID, names(all_pval_after_sentrix))]
  eigenschaften$all_pval_after_bigbatch = all_pval_after_bigbatch[match_hk(eigenschaften$probeID, names(all_pval_after_bigbatch))]
  eigenschaften$all_pval_after_strangebatch = all_pval_after_strangebatch[match_hk(eigenschaften$probeID, names(all_pval_after_strangebatch))]
  eigenschaften$all_pval_before_sentrix = all_pval_before_sentrix[match_hk(eigenschaften$probeID, names(all_pval_before_sentrix))]
  eigenschaften$all_pval_before_bigbatch = all_pval_before_bigbatch[match_hk(eigenschaften$probeID, names(all_pval_before_bigbatch))]
  eigenschaften$all_pval_before_strangebatch = all_pval_before_strangebatch[match_hk(eigenschaften$probeID, names(all_pval_before_strangebatch))]


  ############################################################################### identifizieren der Ueberinflationierten Eigenschaften der Adjustierung merken
  bonf_pwert <- 0.05/length(all_pval_after_sentrix)
  bonf_pwert

  all_pval = eigenschaften[, 2:7]
  head(all_pval)

  all_pval_bonf = all_pval <= bonf_pwert
  colnames(all_pval_bonf) = paste0(colnames(all_pval_bonf), "_bf")
  ht(all_pval_bonf, 3)
  colSums(all_pval_bonf)
  colSums(all_pval_bonf)/(dim(all_pval_bonf)[1])



  eigenschaften = data.frame(cbind(eigenschaften, all_pval_bonf))

  ############################################################################### anfuegen an probe annotation

  qlist55 = venn2(eigenschaften$probeID, genesdetail$nuid, plotte = F)
  for (i in setdiff(colnames(eigenschaften), "probeID")) {
    genesdetail[i] = eigenschaften[match_hk(genesdetail$nuid, eigenschaften$probeID), i]
  }

  ht(genesdetail, 2)
  xtabs(~genesdetail$all_pval_after_sentrix_bf + genesdetail$all_pval_after_bigbatch_bf + genesdetail$is_purecontrol)
  # xtabs(~genesdetail$all_pval_after_sentrix_bf + genesdetail$all_pval_after_bigbatch_bf +genesdetail$expressed_Leber ) xtabs(~genesdetail$all_pval_before_sentrix_bf
  # + genesdetail$all_pval_before_bigbatch_bf +genesdetail$expressed_Leber )


  ############################################################## definition guter probe inkl special batch


  for (i in subgroups) {
    # i = subgroups[1]
    message("Define probes with not over-inflated ANOVA for subgroup ", i, "...\n")
    expressed_var = paste0("expressed_", reformate_subgroup(i))
    goodexpressedprobe_var = paste0("goodexpressedprobe_", reformate_subgroup(i))



    genesdetail[goodexpressedprobe_var] = genesdetail$all_pval_after_sentrix_bf == F & genesdetail$all_pval_after_bigbatch_bf == F & genesdetail$all_pval_after_strangebatch_bf ==
      F & genesdetail$is_purecontrol == F & genesdetail[, expressed_var] == T

    mytable(genesdetail[goodexpressedprobe_var])
  }



  ### Ausgabe
  message("Created followin new sample attributes:\n", paste(setdiff(names(genesdetail), orinames), collapse = "\n"))
  erg = c()
  erg$bonf_pwert = bonf_pwert
  erg$genesdetail = genesdetail
  erg


}
