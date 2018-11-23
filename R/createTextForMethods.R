#' @title Create Text for Material and Methods to report preprocessing
#' @description Creates Text for Material and Methods to report preprocessing
#' @param ht12object A list object of class HT12prepro
#' @param paramfile Path to the file specifying parameters
#' @return an object with a markdown formated short and long text describing preprocessing

#' @import data.table
#' @export

## debug
# ht12object = prepro_ht12; paramfile=myparamfile

createTextForMethods <- function(ht12object,paramfile) {

 param <- data.frame(data.table::fread(paramfile))

  chipsamples = ht12object$chipsamples
  data.table::setDT(chipsamples)

  sample_overview_l3 = ht12object$dokuobjects_createExpressionSet$sample_overview_l3
  data.table::setDT(sample_overview_l3)
  givenind_n = sample_overview_l3[is.na(reason4exclusion),.N ]

  noncorr = ht12object$dokuobjects_filterTechnicallyFailed$noncorr
  head(categmahal)
  stopifnot(all(noncorr %in% categmahal$short))
  categmahalused = categmahal[ short %in% noncorr,sort(unique(long))]

  bad_mahal_inds = ht12object$dokuobjects_filterTechnicallyFailed$bad_mahal_inds
  badmahal_n = length(bad_mahal_inds)

  singlbarcoders_ind = ht12object$dokuobjects_filter4MinBatchsize$singlbarcoders_ind
  singlbarcoders_ind_n = length(singlbarcoders_ind)
  num_bad_euklid = ht12object$dokuobjects_filterAtypicalExpressed$num_bad_euklid
  badeuklid_n = sum(suppressWarnings(na.omit(as.numeric(unlist(num_bad_euklid)))))


  probeannot =ht12object$genesdetail
  data.table::setDT(probeannot)

  all_probes_n = dim(probeannot)[1]



  all_geneprobes_n = sum(grepl("transcr", probeannot$Probe_Class))
  all_controlprobes_n = sum(grepl("control", probeannot$Probe_Class))
  all_geneprobes_extracted_n = probeannot[is.na(is_purecontrol)==F, sum(grepl("transcr", Probe_Class))]
  all_controlprobes_extracted_n = probeannot[is.na(is_purecontrol)==F, sum(grepl("control", Probe_Class))]


  mytable(probeannot$Probe_Class)

  xtabs_hk(~probeannot$Probe_Class +probeannot$is_purecontrol )

  expressedprobecols <- grep(paste(c("^expressed"), collapse = "|"), names(probeannot), value = T)
  expressedprobecols
  showme <- data.frame(eval(expr = parse(text = paste0("with(probeannot, xtabs_hk(~is_purecontrol + is_ercc +",
                                                       paste(expressedprobecols, collapse = "+"), "))"))))
  showme[order(showme$is_purecontrol, showme$is_ercc,showme$Freq), ]
  probeannot$expressedprobe_allsubgroups = apply(probeannot[, expressedprobecols,
                                                            with = F], 1,
                                                 function(x) all(x == T))


  transcript_expressedprobe_allsubgroups_n = sum(probeannot$expressedprobe_allsubgroups & (probeannot$is_purecontrol ==F), na.rm = T)


  transcript_excluded_ilmn = probeannot[(probeannot$all_pval_after_bigbatch_bf | probeannot$all_pval_after_sentrix_bf) & probeannot$is_purecontrol==F, ilmn]
  transcript_excluded_ilmn_n = length(transcript_excluded_ilmn)

  subgroups = chipsamples[,unique(subgroup)]
  zahlsugbrupp = table(sample_overview_l3[ is.na(sample_overview_l3$reason4exclusion),subgroup])

  euklidtextfak = paste0('This was done separately for subgroups ',paste(names(zahlsugbrupp), collapse = " / "), ", initially including ",paste(huebsch(as.numeric(zahlsugbrupp)), collapse = " / "), ' individuals, respectively.')

  expresstext = paste0('This was done separately for subgroups ',  paste(subgroups, collapse = ' / '), ". ")

  combattextfak = paste0('Thereby, subgroups ',  paste(subgroups, collapse = '/ '), ' were included as covariates in order to protect differences between subgroups in the expression data. This was possible as no strong imbalancies between batches and subgroups were observed [PMID:26272994]')
  anovatextfak = paste0('Again, this analysis used subgroups ',  paste(subgroups, collapse = ' / '), ' as covariates.')

  zahlsugbrupp_fin = table(chipsamples[ chipsamples$in_study ==T,subgroup])

  finasubgroupind_fak = paste0('Finally included samples comprised of subgroups ',paste(names(zahlsugbrupp_fin), collapse = " / "), " including ",paste(huebsch(as.numeric(zahlsugbrupp_fin)), collapse = " / "), ' individuals, respectively.')


  anova_with_special_batch = ht12object$dokuobjects_checkBatchEffect$anova_with_special_batch

  more_than_1subgroup = length(subgroups) > 1

  dim_total_nobkgd_eset_ql_combat = dim(ht12object$total_nobkgd_eset_ql_combat)

  bad_genenum = ht12object$dokuobjects_filterLowExpressed$bad_genenum
  numgen_cutoff = ht12object$dokuobjects_filterLowExpressed$numgen_cutoff
  minexprimiert =  ht12object$dokuobjects_filterLowExpressed$minexprimiert
  normmethod = ht12object$dokuobjects_transformNormalizeHT12object$normmethod
  second_combat_withvar = ht12object$dokuobjects_removeBatchEffects$second_combat_withvar_used
  robustmethod_for_mahal = getParam2("robustmethod_for_mahal", myparam = param)
  filter2ind_atypischIlmnKontroll =  ht12object$dokuobjects_filterTechnicallyFailed$mahalparam
  save_subgroupcontrast = ht12object$dokuobjects_removeBatchEffects$save_subgroupcontrast_used
  outlierCriteriumEuklid  =ht12object$dokuobjects_filterAtypicalExpressed$outlierCriteriumEuklid
  good_entrez = ht12object$dokuobjects_writeFilesTosend$good_entrez


longtext_1 = paste0(" Raw data of all ",huebsch(all_geneprobes_n)," gene-expression probes and ", huebsch(all_controlprobes_n), " control probes was extracted by Illumina GenomeStudio without additional background correction. ",huebsch(dim_total_nobkgd_eset_ql_combat[1])," probes (",huebsch(all_geneprobes_extracted_n)," gene-expression probes and ",huebsch(all_controlprobes_extracted_n)," control probes) could be successfully imputed in a total of ",huebsch(givenind_n)," included samples.\n Data was further processed within R/ Bioconductor R `[PMID:15461798]`.\n\n")

longtext_2 = paste0("Initially, ", huebsch(length(bad_genenum))," (",proz(length(bad_genenum)/givenind_n),") samples having an extreme number of expressed genes (defined as median ± ",numgen_cutoff," x interquartile ranges (IQR) of the cohort's values) were excluded. ",if(more_than_1subgroup) capture.output(pander::pander(euklidtextfak)), "\nTranscripts not expressed at p = 0.05 (as defined by Illumina and implemented in the R / Bioconductor package `lumi` `[PMID:18467348]`) in at least ",proz(as.numeric(minexprimiert))," of all samples were excluded from further analysis. ", if(more_than_1subgroup) capture.output(pander::pander(expresstext)), huebsch(transcript_expressedprobe_allsubgroups_n), " (", proz(transcript_expressedprobe_allsubgroups_n/all_geneprobes_n), ") probes remained in analysis ",if(more_than_1subgroup) capture.output(pander::pander('in all subgroups.')), ".\n\n")

longtext_3 = paste0("Expression values were ",normmethod, "-normalised and log2-transformed `[PMID:20525181]`.\n\nFurthermore, we defined for each sample a combined quantitative measure combining quality control features available for head-12 v4 (i.e. ",capture.output(pander::pander(categmahalused)), "). We calculated Mahalanobis-distance using R / Bioconductor package `mdqc` `[PMID:17933854]` between all samples and an artificial sample having average values for these quality control features applying robust method `",robustmethod_for_mahal, "` . ", huebsch(badmahal_n), " (", proz(as.numeric(badmahal_n)/givenind_n,2), ") samples with a distance larger than median + ", filter2ind_atypischIlmnKontroll, " x IQR were excluded.\n\n")


longtext_4 = paste0("Transcript levels were adjusted for the known batch Sentrix barcode (i.e. expression chip-ID) ", if(identical(second_combat_withvar, "")==F) capture.output(pander::pander(paste(" and",second_combat_withvar ))), " using an empirical Bayes method as described `[PMID:16632515]`. The empirical Bayes method required that at least two samples for each batch are provided. This excluded ", huebsch(singlbarcoders_ind_n), " (", proz(as.numeric(singlbarcoders_ind_n)/givenind_n,2), ") samples. ", if(more_than_1subgroup & save_subgroupcontrast==T) capture.output(pander::pander(combattextfak)), " Success of adjustment was checked using ANOVA for the Sentrix barcode, ", if(!anova_with_special_batch %in% c(FALSE, "")) capture.output(pander::pander(paste0(", ",anova_with_special_batch, ","))), " as well as the processing batch (in a processing batch, several expression chips were jointly processed, in consequence, within a processing-batch, several Sentrix barcodes are completely nested). ", if(more_than_1subgroup) capture.output(pander::pander(anovatextfak))," A total of ", huebsch(transcript_excluded_ilmn_n), " (", proz(transcript_excluded_ilmn_n/all_geneprobes_n), ") gene-expression probes still over-inflated following Bonferroni-correction were excluded from following preprocessing steps and ", huebsch(sum(probeannot$goodexpressedprobe_allsubgroups, na.rm=T)), " (", proz(sum(probeannot$goodexpressedprobe_allsubgroups, na.rm=T)/all_geneprobes_n,2), ") gene-expression probes remained,", if(more_than_1subgroup) capture.output(pander::pander(' in all subgroups')), ".\n\n")

longtext_4b =paste0("For further outlier detection, we calculated the Euclidian distance between all samples ", if(more_than_1subgroup) capture.output(pander::pander('of the same subgroup')), " and an artificial sample. This sample was defined as the average of samples after removing 10% samples farthest away from the average of all samples  ", if(more_than_1subgroup) capture.output(pander::pander('done separately for each subgroup')), " (implemented in the R / Bioconductor package `lumi` `[PMID:18467348]`). ", huebsch(badeuklid_n), " (", proz(as.numeric(badeuklid_n)/givenind_n), ") samples with a distance larger than median + ", outlierCriteriumEuklid, " x IQR were excluded resulting in ", huebsch(sum( chipsamples$in_study, na.rm=T))," remaining samples. ",if(more_than_1subgroup) capture.output(pander::pander(finasubgroupind_fak)), "\n\n")


longtext_5 = paste0("Mapping of genes corresponding to expression probes and assignment of gene names was done using information of a remapping approach `[PMID:19923232]` applying gene-information of the Entrez gene database of the National Center for Biotechnology Information (NCBI), available at http://www.ncbi.nlm.nih.gov/gene/.  This information was retrieved using the R add-on package from Bioconductor `illuminaHumanv4.db_1.14.0` that relates to NCBI data dated on 2012-March-7 (hg19).\n\n")

longtext_5b = paste0("This remapping approach resulted in a total of ",  huebsch(sum(probeannot$perfectprobe_allsubgroups, na.rm=T)), " valid gene-expression probes corresponding to ",huebsch(length(good_entrez)), " valid unique genes", if(more_than_1subgroup) capture.output(pander::pander(' available in all subgroups')), ". Thereby, the probe annotation quality score `[PMID:19923232]` had to be at least _good_.\n\nEntrez-gene IDs were used to retrieve information for the abbreviated gene names (HGNC identifier `[PMID:25361968]`) and transcription start- and end- site of corresponding genes via Bioconductor package `org.Hs.eg.db_2.7.1`. This package is based on hg19 coordinates retrieved from Golden Path data provided by UCSC Genome Bioinformatics at ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19 with a date stamp of 2010-Mar22.\n\nChromosomal mapping of genes was done with R / Bioconductor package `org.Hs.eg.db_2.7.1`. This package is based on hg19 coordinates retrieved from Golden Path data provided by UCSC Genome Bioinformatics at ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19 with a date stamp of 2010-Mar22. Mapping information `[PMID:19923232]` was only accepted if distance between chromosomal coordinates of expression probes and TSS / TSE was smaller than 1Mb.  We found ", huebsch(sum(probeannot$perfectprobe_allsubgroups & probeannot$mappings == 1 & is.na(probeannot$mappings)==F, na.rm=T)), " probes", if(more_than_1subgroup) capture.output(pander::pander(' available in all subgroups')), ", corresponding to ",huebsch(length(unique(na.omit(probeannot[probeannot$perfectprobe_allsubgroups & probeannot$mappings == 1 & is.na(probeannot$mappings)==F,EntrezReannotated])))), " genes, to map to a single position in the human genome (hg19).")


longtext = paste0(longtext_1, longtext_2, longtext_3,longtext_4,longtext_4b, longtext_5, longtext_5b)



shorttext_1 = paste0("Raw data of all ", huebsch(all_geneprobes_n), " gene-expression probes was extracted by Illumina GenomeStudio without additional background correction. Data was further processed within R / Bioconductor. Expression values were log2-transformed and ", normmethod,"-normalised `[PMID:18467348] [PMID:20525181]`. Batch effects of expression BeadChips", if(identical(second_combat_withvar, "")==F) paste(" and",second_combat_withvar ), " were corrected using an empirical Bayes method `[PMID:16632515]`.\n\n")

shorttext_2 = paste0("Within pre-processing, gene-expression probes detected by Illumina GenomeStudio as expressed in less than ", proz(as.numeric(minexprimiert)), " of the ",  if(more_than_1subgroup) pander("subgroup-specific"), "samples were excluded as well as probes still found to be significantly associated with batch effects after Bonferroni-correction. Furthermore, gene-expression probes with poor mapping on the human trancriptome `[PMID:19923232]`  were also excluded.\nIn summary, these filters resulted in ", huebsch(sum(probeannot$perfectprobe_allsubgroups, na.rm=T)), " valid gene-expression probes corresponding to ",  huebsch(length(good_entrez)), " unique genes", if(more_than_1subgroup) capture.output(pander::pander(' found in all subgroups')), ". Among those probes, ", huebsch(sum(probeannot$perfectprobe_allsubgroups & probeannot$mappings == 1 & is.na(probeannot$mappings)==F, na.rm=T)), " probes corresponding to ", huebsch(length(unique(na.omit(probeannot[probeannot$perfectprobe_allsubgroups & probeannot$mappings == 1 & is.na(probeannot$mappings)==F,EntrezReannotated])))), " genes mapped to a single position in the human genome (hg19).\n\n")

shorttext_3 = paste0("Three criteria were used to remove samples of low quality: First, the number of detected gene-expression probes of a sample was required to be within ± ", numgen_cutoff, " interquartile ranges (IQR) from the median. Second, the Mahalanobis distance of several quality characteristics of each sample (", capture.output(pander::pander(categmahalused)), ")`[PMID:17933854]` had to be within median + ", filter2ind_atypischIlmnKontroll, " x IQR. Third, Euclidean distances of expression values as described `[PMID:18467348]` had to be within ", outlierCriteriumEuklid, " x IQR from the median.\nOverall, of the assayed ", huebsch(givenind_n), " samples, ", givenind_n -sum( chipsamples$in_study, na.rm=T), " samples were excluded for quality reasons.  ", if(more_than_1subgroup) capture.output(pander::pander(finasubgroupind_fak)))

shorttext = paste0(shorttext_1, shorttext_2, shorttext_3)

res = c()
res$extended_version = longtext
res$short_version = shorttext
res
}


