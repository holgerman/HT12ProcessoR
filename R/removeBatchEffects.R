  #' @title Remove Batch effects from expressoin data of a HT12prepro object
  #'
  #' @description This functions removes via ComBat from packge sva batch effects origniating from the Chip (Sentrix.Barcode) of the sample.  If second_combat_withvar is provided, a second round of removeBatchEffects() is done.
  #'
  #'

  #' @param ht12object A list object of class HT12prepro created with function filter4MinBatchsize
  #' @param paramfile Path to the file specifying parameters
  #' @param save_subgroupcontrast Protect the subgroup contrast (column subgroup in the table stored in the slot $chipsamples of the HT12prepro object)from beeing affected in batch ajdustment, CAVE risk of false positives when subgroups are unbalanced with batches, hence FALSE is in this case recommended. If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.

  #' @param excludeERCC exclude ERCC samples from following steps, this is suggested if some samples have spiked-in ERCC samples but others have not. If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
  #' @param second_combat_withvar A second round of removeBatchEffects() can be intended. In this case, batches from the first round checking Sentrix.Barcode batches as well as batches from this second round specified as a column in ht12object$chipsamples named according to this variable are removed sequentially. If "", no 2nd combat will be done. If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
  #' @param additional_normalisation_after_combat Indicator (TRUE or FALSE), whether to normalize data after batch-correction with Combat using the normalization method specified in parameter 'normalisation_method'. If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
  #' @param normalisation_method Method used for normalisation if parameter additional_normalisation_after_combat = TRUE. Either 'quantile' (quantile normalisation) or 'rsn' (robust spline normalisation). If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.
  #' @return A list object of class HT12prepro where a slot with an expression set with normalized and transformed and batch-corrected data including control probe information is created. This slot is named `$total_nobkgd_eset_ql_combat`. Additionally, a graph is shown showing probe-wise as well as sample-wise correlation of expression data before and after batch adjustment.

  #' @import data.table
  #' @export

  ## debug
  # paramfile = "/mnt/ifs1_projekte/genstat/02_projekte/1704_boettcher_ge_ht12/01_prepro/input_parameter_007.txt"
  # save_subgroupcontrast = "from_paramfile"
  # excludeERCC = "from_paramfile"
  # ht12object =  prepro_ht12
  # second_combat_withvar= "from_paramfile"
  # normalisation_method = "from_paramfile"
  # additional_normalisation_after_combat= "from_paramfile"
  # showPlot_additional_normalisation_after_combat=T

  removeBatchEffects = function(ht12object,
                                paramfile = NULL,
                                excludeERCC = "from_paramfile",
                                save_subgroupcontrast = "from_paramfile",
                                second_combat_withvar= "from_paramfile",
                                additional_normalisation_after_combat= "from_paramfile",
                                normalisation_method = "from_paramfile",
                                showPlot_additional_normalisation_after_combat=T) {

    # ht12object=prepro_ht12;paramfile = myparamfile = paste0(prepro_folder, "/input_parameter_010.txt");excludeERCC = "from_paramfile"; save_subgroupcontrast = "from_paramfile"; second_combat_withvar= "from_paramfile"; additional_normalisation_after_combat= "from_paramfile"; normalisation_method = "from_paramfile"; showPlot_additional_normalisation_after_combat=T


    ### strings are imported as strings and not as factors
    options(stringsAsFactors=FALSE)

    myparameters = match.call()
    showVennplots = F

    # status checken
    historie =  ht12object$history$calls
    if(any(grepl("filter4MinBatchsize", historie))==F) stop("Function 'filter4MinBatchsize()' has to be run before!")

    #laden parameter
    if(is.null(paramfile)==F) param <- data.frame(data.table::fread(paramfile))




    sample_overview_l7 <-ht12object$chipsamples
    mytable(sample_overview_l7$in_study)
    sample_overview_l7instudy <- sample_overview_l7[ sample_overview_l7$in_study, ]
    dim(sample_overview_l7)
    dim(sample_overview_l7instudy)
    table(table(sample_overview_l7instudy$new_ID))
    if(length(table(table(sample_overview_l7instudy$new_ID))) != 1)
      stop("IDs (column new_ID) must be unique....stopping...")

    # laden annotation probes

    genesdetail <- ht12object$genesdetail
    ht(genesdetail, 2)

    #laden expressionsets
    total_nobkgd_eset_ql = ht12object$total_nobkgd_eset_ql

    ## ----goodind-------------------------------------------------------------
    # auszahlen aktualisiere barcode
    table(table(sample_overview_l7$Sentrix.Barcode))
    # table(sample_overview_l7[grep("mis", sample_overview_l7$Sentrix.Barcode), 5])
    Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber <- sample_overview_l7[match_hk(Biobase::pData(total_nobkgd_eset_ql)$sampleID,
                                                                                                       sample_overview_l7$new_ID), "Sentrix.Barcode" ]
    check <- showNA(Biobase::pData(total_nobkgd_eset_ql))$NAs
    if(sum(check) != 0)
      stop("not all smples appear to ahave a Sentrix ID!")
    ht(Biobase::pData(total_nobkgd_eset_ql), 1)

    # gute individuen
    goodind  <- sample_overview_l7instudy$new_ID
    # str(goodind)
    total_nobkgd_eset_ql
    total_nobkgd_eset_ql <- total_nobkgd_eset_ql[,goodind]
    total_nobkgd_eset_ql

    ## ----criteria2-----------------------------------------------------------
    tabled <- table(Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber)
    table(tabled)
    singlbarcoders <- names(tabled[tabled == 1])
    singlbarcoders
    singlbarcoders_ind <- sample_overview_l7instudy[sample_overview_l7instudy$Sentrix.Barcode %in% singlbarcoders, "new_ID"]
    singlbarcoders_ind
    # str(goodind)
    goodind <- setdiff(goodind, singlbarcoders_ind)
    # str(goodind)
    total_nobkgd_eset_ql <- total_nobkgd_eset_ql[,goodind]
    total_nobkgd_eset_ql



    ## ----throwercc-----------------------------------------------------------

    if(excludeERCC== "from_paramfile") excludeERCC_used <- getParam2("excludeERCC", myparam = param) else excludeERCC_used = excludeERCC
    excludeERCC_used <- as.logical(excludeERCC_used)
    excludeERCC_used

    if(excludeERCC_used == T){
      ht(genesdetail,1)
      ercc <- na.omit(genesdetail[genesdetail$is_ercc == T, "nuid"])
      ercc
      total_nobkgd_eset_ql
      total_nobkgd_eset_ql <- total_nobkgd_eset_ql[rownames(Biobase::exprs(total_nobkgd_eset_ql)) %nin% ercc,]
      total_nobkgd_eset_ql
    }

    ## ----combat8b------------------------------------------------------------
    head(Biobase::pData(total_nobkgd_eset_ql))
    Biobase::pData(total_nobkgd_eset_ql)$subgroup <- sample_overview_l7[ match_hk(Biobase::pData(total_nobkgd_eset_ql)$sampleID,
                                                                                  sample_overview_l7$new_ID), "subgroup"]
    showNA(Biobase::pData(total_nobkgd_eset_ql))
    singularcheck <- data.table(Biobase::pData(total_nobkgd_eset_ql))
    singularcheck2 <- singularcheck[,.N, by = list(hybridisierungchipserialnumber, subgroup) ]
    singularcheck3 <- singularcheck2[allDuplicatedEntries(hybridisierungchipserialnumber)]
    singularcheck3

    if(save_subgroupcontrast== "from_paramfile") save_subgroupcontrast_used <- getParam2("save_subgroupcontrast", myparam = param) else save_subgroupcontrast_used = save_subgroupcontrast

    save_subgroupcontrast_used = as.logical(save_subgroupcontrast_used)


    # check whether batches are present and skip combat if only one batch is foud
    if(length( unique(Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber))==1) {
      message("Skipping adjustment for Chip batcheffect as only a single chip is included...\n")
      total_nobkgd_eset_ql_combat = total_nobkgd_eset_ql
      status_combatMitCovar = F
    } else {


      # ComBat without subgroup =================================================

      if(nrow(singularcheck3) == 0 | save_subgroupcontrast_used == F) {
        message("NOT additionally protecting contrast `subgroup` in Combat...")

        batch <- Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber

        # ComBat with covariates --------------------------------------------------

        # change modcombat based on covars provided in params file
        if(as.logical(getParam2("ComBat_adjust_covariates", myparam = param))) {

          # add columns provided to modcombat
          combat.covars <- getParam2("ComBat_covariates_to_adjust", myparam = param)

          # format
          combat.covars <- gsub(pattern = " ", replacement = "", x = combat.covars)
          combat.covars <- unlist(strsplit(x = combat.covars, split = ",", fixed = T))

          # get and match data
          # Add data
          pData(total_nobkgd_eset_ql)[, combat.covars] <- NA
          Biobase::pData(
            total_nobkgd_eset_ql)[
              , combat.covars
              ] <-  sample_overview_l7[
                match_hk(
                  Biobase::pData(
                    total_nobkgd_eset_ql
                    )$sampleID, sample_overview_l7$new_ID
                  ), combat.covars
                ]

          # filter the individuals for NAs in the phenotypic data and set a filter for removal of those inds both in the phenotypes and in the gx data
          good.ids <- na.omit(Biobase::pData(total_nobkgd_eset_ql)[, c("sampleID", combat.covars)])$sampleID
          batch <- Biobase::pData(total_nobkgd_eset_ql)[pData(total_nobkgd_eset_ql)$sampleID %in% good.ids, "hybridisierungchipserialnumber"]
          message("Removing the individuals ", paste(
            setdiff(pData(total_nobkgd_eset_ql)$sampleID, good.ids), collapse = ", "), " due to missings in the coviarate data, which is not permitted in ComBat.")

          # sample_overview_l7[
          #       match_hk(
          #         Biobase::pData(
          #           total_nobkgd_eset_ql
          #           )$sampleID, sample_overview_l7$new_ID
          #         ), c("in_study", "subgroup")]

          # create formula
          combat.formula <- formula(
            paste0("~ 1 + ", paste(combat.covars, collapse = " + "))
            )

          # add the data in the model matrix
          modcombat <- model.matrix(combat.formula, data = Biobase::pData(total_nobkgd_eset_ql))
          status_combatMitCovar = T

          # ComBat without Covariates -----------------------------------------------

        } else {

          # empty filter in case of no covariates, i.e. nothing to filter for
          good.ids <- Biobase::pData(total_nobkgd_eset_ql)$sampleID

          # use the model matrix without covariates
          modcombat <- model.matrix(~1, data = Biobase::pData(total_nobkgd_eset_ql))
          status_combatMitCovar = F
        }

        # run combat with or without covariates depending on previous if-else
        combat_edata <- sva::ComBat(dat=Biobase::exprs(total_nobkgd_eset_ql)[, good.ids], # this filters out samples with missings in covariate data
                                    batch = batch,
                                    mod = modcombat,
                                    par.prior = TRUE,
                                    prior.plots = T)
        total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql
        stopifnot(identical(rownames(combat_edata), rownames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        stopifnot(identical(colnames(combat_edata), colnames(Biobase::exprs(total_nobkgd_eset_ql_combat)[,good.ids])))
        total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql_combat[, good.ids]
        Biobase::exprs(total_nobkgd_eset_ql_combat)[, good.ids] <- combat_edata[, good.ids]

      } else {


        # ComBat with Subgroup ====================================================

        # ComBat with Covariates --------------------------------------------------

        # batch adjust including subgroup with and without covariates
        # change modcombat based on covars provided in params file
        if(as.logical(getParam2("ComBat_adjust_covariates", myparam = param))) {

          # add columns provided to modcombat
          combat.covars <- getParam2("ComBat_covariates_to_adjust", myparam = param)

          # format
          combat.covars <- gsub(pattern = " ", replacement = "", x = combat.covars)
          combat.covars <- unlist(strsplit(x = combat.covars, split = ",", fixed = T))

          # get and match data
          # Add data
          pData(total_nobkgd_eset_ql)[, combat.covars] <- NA
          Biobase::pData(
            total_nobkgd_eset_ql)[
              , combat.covars
              ] <-  sample_overview_l7[
                match_hk(
                  Biobase::pData(
                    total_nobkgd_eset_ql
                  )$sampleID, sample_overview_l7$new_ID
                ), combat.covars
                ]

          # filter the individuals for NAs in the phenotypic data and set a filter for removal of those inds both in the phenotypes and in the gx data
          good.ids <- na.omit(Biobase::pData(total_nobkgd_eset_ql)[, c("sampleID", combat.covars)])$sampleID
          batch <- Biobase::pData(total_nobkgd_eset_ql)[pData(total_nobkgd_eset_ql)$sampleID %in% good.ids, "hybridisierungchipserialnumber"]
          message("Removing the individuals ", paste(
            setdiff(
              pData(
                total_nobkgd_eset_ql
                )$sampleID, good.ids
              ), collapse = ", "
            ), " due to missings in the coviarate data, which is not permitted in ComBat.")


          # create formula
          combat.formula <- formula(
            paste0("~ 1 + subgroup + ", paste(combat.covars, collapse = " + "))
          )

          # add the data in the model matrix
          modcombat <- model.matrix(combat.formula, data = Biobase::pData(total_nobkgd_eset_ql))
          status_combatMitCovar <- "subgroup + covars"

        } else {

          # ComBat without Covariates -----------------------------------------------

          # empty filter in case of no covariates, i.e. nothing to filter for
          good.ids <- Biobase::pData(total_nobkgd_eset_ql)$sampleID
          batch <- Biobase::pData(total_nobkgd_eset_ql)$hybridisierungchipserialnumber

          # formula with additional subgroup
          combat.formula <- formula("~ 1 + subgroup")
          status_combatMitCovar <- "subgroup"

          # use the model matrix without covariates
          modcombat <- model.matrix(combat.formula, data = Biobase::pData(total_nobkgd_eset_ql))
        }

        # ComBat Batch adjustment -------------------------------------------------

        message("Additional protecting contrast `subgroup` in Combat...")
        # modcombat <- model.matrix(~subgroup, data=Biobase::pData(total_nobkgd_eset_ql))
        subgroup <- Biobase::pData(total_nobkgd_eset_ql)[pData(total_nobkgd_eset_ql)$sampleID %in% good.ids, "subgroup"]
        balanciertcheck <- xtabs_hk(~ batch + subgroup)
        # balanciertcheck <- xtabs_hk(~ batch + Biobase::pData(total_nobkgd_eset_ql)$subgroup)
        message(balanciertcheck)
        message(chisq.test(balanciertcheck))
        if(chisq.test(balanciertcheck)$p.value <= 0.05)
          warning("subgroups not balanced, but contrast saved in combat - might result in false positives...", immediate. = T)
        combat_edata <- sva::ComBat(dat = Biobase::exprs(total_nobkgd_eset_ql)[, good.ids],
                                    batch = batch,
                                    mod = modcombat,
                                    par.prior = TRUE,
                                    prior.plots = T)
        total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql
        stopifnot(identical(rownames(combat_edata), rownames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        stopifnot(identical(colnames(combat_edata), colnames(Biobase::exprs(total_nobkgd_eset_ql_combat)[, good.ids])))

        # replace original data with adjusted data
        total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql_combat[, good.ids]
        Biobase::exprs(total_nobkgd_eset_ql_combat)[, good.ids] <- combat_edata[, good.ids]

      }
    }

    if(second_combat_withvar== "from_paramfile") second_combat_withvar_used <- getParam2("second_combat_withvar", myparam = param) else second_combat_withvar_used = second_combat_withvar
    second_combat_withvar_used

    if(second_combat_withvar_used != "") {
      message("Adjusting genexpression data in a second combat round for `",
              second_combat_withvar_used, "`...")

      # Add data
      Biobase::pData(total_nobkgd_eset_ql_combat)[,second_combat_withvar_used] <-  sample_overview_l7[match_hk(Biobase::pData(total_nobkgd_eset_ql_combat)$sampleID,
                                                                                                               sample_overview_l7$new_ID),
                                                                                                      second_combat_withvar_used]
      check <- showNA(Biobase::pData(total_nobkgd_eset_ql_combat))
      if(sum(check) != 0)
        stop("Not all individuals have a Sentrix ID!")
      ht(Biobase::pData(total_nobkgd_eset_ql_combat), 1)

      # Exclude single barcoders
      tabled <- table(Biobase::pData(total_nobkgd_eset_ql_combat)[, second_combat_withvar_used])
      table(tabled)
      newsinglbarcoders <- names(tabled[tabled == 1])
      newsinglbarcoders
      newsinglbarcoders_ind <- sample_overview_l7instudy[sample_overview_l7instudy$Sentrix.Barcode %in% newsinglbarcoders, "new_ID"]
      newsinglbarcoders_ind
      # str(goodind)
      goodind <- setdiff(goodind, newsinglbarcoders_ind)
      # str(goodind)
      total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql_combat[, goodind]
      total_nobkgd_eset_ql_combat
      singlbarcoders_ind <- unique(c(singlbarcoders_ind, newsinglbarcoders_ind))
      singlbarcoders <- unique(c(singlbarcoders, newsinglbarcoders))

      ## singularcheck 2nd round
      singularcheck2ndround <- data.table(Biobase::pData(total_nobkgd_eset_ql_combat))
      singularcheck2ndround2 <- singularcheck2ndround[,.N, by = list(get(second_combat_withvar_used),
                                                                     subgroup)]
      singularcheck2ndround3 <- singularcheck2ndround2[allDuplicatedEntries(get)]
      singularcheck2ndround3

      # combat 2nd round
      if(nrow(singularcheck2ndround3) == 0 | save_subgroupcontrast_used == F) {
        message("NOT additional protecting contrast `subgroup` in Combat in 2nd round of adjusting with",
                second_combat_withvar_used,"...")
        batch <- Biobase::pData(total_nobkgd_eset_ql_combat)[, second_combat_withvar_used]
        modcombat <- model.matrix(~1, data = Biobase::pData(total_nobkgd_eset_ql_combat))
        combat_edata <- sva::ComBat(dat = Biobase::exprs(total_nobkgd_eset_ql_combat),
                                    batch = batch,
                                    mod = modcombat,
                                    par.prior = TRUE,
                                    prior.plots = T)
        stopifnot(identical(rownames(combat_edata), rownames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        stopifnot(identical(colnames(combat_edata), colnames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        Biobase::exprs(total_nobkgd_eset_ql_combat) <- combat_edata
        status_2ndcombatMitCovar <- F
        # also falls bei adjustierung contrast subgroup protected werden soll:
      } else {
        message("Additional protecting contrast `subgroup` in Combat...")
        batch <- Biobase::pData(total_nobkgd_eset_ql_combat)[, second_combat_withvar_used]
        modcombat <- model.matrix(~subgroup, data = Biobase::pData(total_nobkgd_eset_ql_combat))
        balanciertcheck <- xtabs_hk(~batch + Biobase::pData(total_nobkgd_eset_ql_combat)$subgroup)
        message(balanciertcheck)
        message(chisq.test(balanciertcheck))
        stopifnot(chisq.test(balanciertcheck)$p.value > 0.05)
        combat_edata <- sva::ComBat(dat = Biobase::exprs(total_nobkgd_eset_ql_combat),
                                    batch = batch,
                                    mod = modcombat,
                                    par.prior = TRUE, prior.plots = T)
        stopifnot(identical(rownames(combat_edata), rownames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        stopifnot(identical(colnames(combat_edata), colnames(Biobase::exprs(total_nobkgd_eset_ql_combat))))
        Biobase::exprs(total_nobkgd_eset_ql_combat) <- combat_edata
        status_2ndcombatMitCovar = "subgroup"
      }
    } else { # also falls keine 2. runde von combat stattfinden soll:
      message("Not doing a second adjustment round in combat")
      singularcheck2ndround3 = NA
      status_2ndcombatMitCovar = NA
    }
    total_nobkgd_eset_ql_combat


    ## ----addnorm-------------------------------------------------------------


    if(additional_normalisation_after_combat== "from_paramfile") do_additional_normalisation_after_combat <- getParam2("additional_normalisation_after_combat", myparam = param) else do_additional_normalisation_after_combat = additional_normalisation_after_combat
    do_additional_normalisation_after_combat




    all_cor_values_2nd_norm <- NA
    if(do_additional_normalisation_after_combat){



      if(normalisation_method == "from_paramfile")  methodtransform <- getParam2("normalisation_method", myparam = param) else  methodtransform = normalisation_method
      methodtransform

      total_nobkgd_eset_ql_combat2 <- transformNormalize(total_nobkgd_eset_ql_combat,
                                                         methodtransform = methodtransform,
                                                         dolog2 = F)

      hh(exprs(total_nobkgd_eset_ql_combat2))



      if(showPlot_additional_normalisation_after_combat==T) plotNormTransform(total_nobkgd_eset_ql_combat2, total_nobkgd_eset_ql_combat, png_fn = NULL)

      all_cor_values_2nd_norm <- compareEsetExpressions(total_nobkgd_eset_ql_combat, total_nobkgd_eset_ql_combat2, firstmessage = "Comparing expression set before and after second normalisation after removing batch effects via ComBat...\n", showPlots=showPlot_additional_normalisation_after_combat)
      total_nobkgd_eset_ql_combat <- total_nobkgd_eset_ql_combat2
    }

    ## ----coorellationCombat8b------------------------------------------------
    if(identical(Biobase::pData(total_nobkgd_eset_ql_combat)$sampleID,
                 Biobase::pData(total_nobkgd_eset_ql)[
                   Biobase::pData(
                     total_nobkgd_eset_ql
                     )$sampleID %in% good.ids,
                   ]$sampleID) == F) {
      stop("Something went wrong - code 33443361")
      }



    # hinzufuegen zu transkriptsample_overview_l7uten
    all_cor_values <- compareEsetExpressions(total_nobkgd_eset_ql_combat, total_nobkgd_eset_ql,firstmessage = "Comparing expression set before and after  removing batch effects via ComBat...\n")
    ht(genesdetail, 1)
    # str(all_cor_values$cor_transkript)
    genesdetail$spearmancor_combat <- all_cor_values$cor_transkript[match_hk(genesdetail$nuid,
                                                                             names(all_cor_values$cor_transkript))]
    ht(genesdetail, 5)
    sample_overview_l7$spearmancor_combat <- all_cor_values$cor_ind[match_hk(sample_overview_l7$new_ID,
                                                                             names(all_cor_values$cor_ind))]
    ht(sample_overview_l7, 1)
    head(sample_overview_l7[order(sample_overview_l7$spearmancor_combat), ])

    ## ----save8b--------------------------------------------------------------

    # probeattribs
    dim(genesdetail)
    dim(ht12object$genesdetail)
    ht12object$genesdetail = genesdetail

    # individuenattrib
    sample_overview_l7[sample_overview_l7$in_study & !(sample_overview_l7$new_ID %in% good.ids), ]$reason4exclusion <- "tmp"
    sample_overview_l7[!(sample_overview_l7$new_ID %in% good.ids),"in_study"] = F

    # don't overwrite the other reason4exclusions
    sample_overview_l7[sample_overview_l7$reason4exclusion %in% "tmp" ,"reason4exclusion"] <- "Missings in covariate data. This is not allowed with batch adjustment via ComBat()"
    # sample_overview_l7[!(sample_overview_l7$new_ID %in% good.ids),"reason4exclusion"] <- "Missings in covariate data. This is not allowed with batch adjustment via ComBat()"
    sample_overview_l7[sample_overview_l7$new_ID %in% singlbarcoders_ind,"in_study"] = F
    sample_overview_l7[sample_overview_l7$new_ID %in% singlbarcoders_ind,"reason4exclusion"] <- "Batch/effect correction via ComBat not possible - only 1 Ind. with this Sentrix ID"

    mytable(sample_overview_l7$reason4exclusion)
    mytable(sample_overview_l7$in_study)
    ht12object$chipsamples =  sample_overview_l7

    ht12object$total_nobkgd_eset_ql_combat = total_nobkgd_eset_ql_combat

    # status fuer docku
    dim_total_nobkgd_eset_ql_combat <- dim(total_nobkgd_eset_ql_combat)
    dim_total_nobkgd_eset_ql_combat
    fordoku =c("singularcheck3",
               "singularcheck2ndround3",
               "status_combatMitCovar",
               "status_2ndcombatMitCovar",
               "dim_total_nobkgd_eset_ql_combat",
               "singlbarcoders_ind",
               "singlbarcoders",
               "sample_overview_l7",
               "genesdetail",
               "all_cor_values",
               "all_cor_values_2nd_norm",
               'second_combat_withvar_used',
               'save_subgroupcontrast_used')




    stopifnot(sum(duplicated(fordoku))==0)


    ht12object$dokuobjects_removeBatchEffects = lapply(fordoku, function(x) get(x))


    names(ht12object$dokuobjects_removeBatchEffects) = fordoku


    ht12object$history = rbind(ht12object$history, data.frame(calls = paste(Sys.time(), deparse(myparameters))))
    ht12object$history
    ht12object


  }

