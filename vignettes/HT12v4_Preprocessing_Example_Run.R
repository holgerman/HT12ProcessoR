## ----knitr, cache = F, results = "hide", echo = F, warning = T----------------
# save options
options.backup <- options()

# change line wrap to something more readable
options(width = 80)

# Set the global knitr options
knitr::opts_chunk$set(cache = F,
                      results = "hide",
                      echo = T,
                      include = T,
                      message = T,
                      warning = T,
                      fig.width = 9,
                      fig.height = 6)

## ----initiate-----------------------------------------------------------------
# start time measurement
time.start <- Sys.time()

# load required packages
for (i in c(
  "HT12ProcessoR",
  "toolboxH",
  "lumi",
  "knitr",
  "here",
  "data.table",
  "evaluate",
  "plotly",
  "ggplot2",
  "pander"
)) {
  suppressPackageStartupMessages(library(i, character.only = TRUE))
}

# set some basic options
panderOptions('table.split.table', Inf)
options(stringsAsFactors=FALSE)

## ----prepare.examples, echo = T-----------------------------------------------
# check the root of your project
# here::dr_here(show_reason = T) #  the here package makes the navigation in directories easier but works best from within R-Projects

# get the project root directory
# prepro.folder <- here() 
# setwd(/path/to/project)
prepro.folder <- getwd() # please choose a working directory

# create the object with the names of the input files ##
# IDs can be chosen freely, but using terms identifying the data files might be advisable 
my.fileset_id <- c('590-651', 'HSS106-HSS155')

# create the results-folder if necessary
dir.create(paste0(prepro.folder, "/tosend/"))

# also create an obj/ folder for intermediary results
dir.create(paste0(prepro.folder, "/obj/"))

# control probe data file
my.con_f <- c(system.file("extdata", "Run12_ControlProbeProfile_590-651.txt", package = "HT12ProcessoR"), # control probe measurements + ercc (spike-in -- class of control probes)
              system.file("extdata", "ControlProbeProfile_HSS106-HSS155.txt", package = "HT12ProcessoR")) #  second example (processing batch 1 & 2)

# gene probe file
my.nobkgd_f <- c(system.file("extdata", "Run12_nonorm_nobkgd_590-651.txt", package = "HT12ProcessoR"), # raw probe data of human gene expression 
                 system.file("extdata", "Einzelanalyse_nonorm_nobkgd_HSS106-HSS155.txt", package = "HT12ProcessoR"))

# sample file
my.sample_f <- c(system.file("extdata", "Run12_SamplesTable_590-651.txt", package = "HT12ProcessoR"), # sample annotation (each measured indiviuum)
                 system.file("extdata", "SamplesTable_HSS106-HSS155.txt", package = "HT12ProcessoR"))

# consolidate in input-file
input <- data.table(fileset_id = my.fileset_id, 
                    con_f = my.con_f,
                    nobkgd_f = my.nobkgd_f,
                    sample_f = my.sample_f)

# how should the input file be called?
input_filename <- paste0(prepro.folder, "/obj/01_file_overview.txt")

# save the input file for later use
write.table(x = input,
            file = input_filename,
            quote = F,
            col.names = T,
            row.names = F,
            sep = "\t")

## ----parameter.file-----------------------------------------------------------

myparams <- read.delim(system.file("extdata", "input_parameter_010.txt", package = "HT12ProcessoR"),
                       as.is = T)

## ----customize.parameter.file-------------------------------------------------
# where is the input file we created located?
myparams[myparams$variable == "file_names_of_files_and_folders", "value"] <- input_filename

# where should the results be saved?
myparams[myparams$variable == "datafolder_results", "value"] <- paste0(prepro.folder, "/tosend/")

# where should the file be saved (overwrite the previously read-in file)
myparamfile <- paste0(prepro.folder, "/obj/input_parameter_010.txt")

# save the file again with the above changes
write.table(myparams,
            myparamfile,
            quote = F,
            col.names = T,
            row.names = F,
            sep="\t")

## ----pre.pro.step.1-----------------------------------------------------------
# execute function for the checks described above
prepro_ht12 <- checkExtractChipsamples(paramfile = myparamfile) # check the output of the function for details

# print object in HTML-report
pander(prepro_ht12$history)

# check the class of the object
class(prepro_ht12)

# another file preview
pander(ht(prepro_ht12$chipsamples, 1))

## ----pre.pro.step.2-----------------------------------------------------------
# preview first and last lines
pander(ht(prepro_ht12$chipsamples, 2))

# define samples to be excluded, e.g. when different projects are processed on the same chip
exclude.this <- c("614", "638", "HSS 108", "HSS 112", "HSS 115", "HSS 125", "HSS 127", "HSS 133", "HSS 137", "HSS 138", "HSS 146", "HSS 151","HSS 155")

# Exclude Samples by modifying the $chipsamples element
prepro_ht12$chipsamples <- prepro_ht12$chipsamples[!prepro_ht12$chipsamples$old_ID  %in% exclude.this, ]

# Define Subgroups -  in this example, sampels are from Peripheral Blood Mononuclear Cells (PBMC)
prepro_ht12$chipsamples$subgroup <- "PBMC"
table(prepro_ht12$chipsamples$subgroup)

# show which IDs were processed in which batch
xtabs(~prepro_ht12$chipsamples$Sentrix.Barcode + prepro_ht12$chipsamples$processingbatch)

# do we have strange batches on top of regular batch structure
xtabs(~prepro_ht12$chipsamples$strangebatch + prepro_ht12$chipsamples$processingbatch)

## Quick Check for major differences based on Attributes from the provided ILLUMINA-sample-files
initial_pca <- calcInitialPCA(prepro_ht12)

# plot the results
pca.plot <- ggplot(initial_pca$scores, aes(PC1,
                                           PC2,
                                           pch = subgroup,
                                           col = processingbatch,
                                           size = Detected.Genes..0.01.,
                                           alpha = in_study,
                                           label = paste(old_ID, new_ID, sep = "\n"))) +
  geom_point() +
  scale_alpha_manual(values = c(0.7, 0.3))

# plot a static version
pca.plot

# interactive visualisation - check homogeinity etc.
pca.plot %>% ggplotly()

## ----pre.pro.step.3-----------------------------------------------------------
# create the ExpressionSet 
prepro_ht12 <- createExpressionSet(prepro_ht12, paramfile = myparamfile) # watch the function output

# explore file output
# show the command history
kable(prepro_ht12$history)

# preview expressionset including samples and controls
prepro_ht12$total_nobkgd_eset

# detection levels per probe
hh(lumi::detection(prepro_ht12$total_nobkgd_eset))

# expression lvl per probe
hh(lumi::exprs(prepro_ht12$total_nobkgd_eset))

## ----pre.pro.step.4-----------------------------------------------------------
# execute step 4 by running this function
prepro_ht12 <- filterLowExpressed(prepro_ht12, paramfile = myparamfile)

# show command history
kable(prepro_ht12$history)

# preview function output
head(prepro_ht12$genesdetail)

# show excluded samples
setDT(prepro_ht12$chipsamples)
prepro_ht12$chipsamples[,.N,.(in_study, reason4exclusion)]

# and reset class of this object
setDF(prepro_ht12$chipsamples)

## ----pre.pro.step.5-----------------------------------------------------------
# Transform and normalize the data with this function
prepro_ht12 <- transformNormalizeHT12object(prepro_ht12, paramfile = myparamfile)

# execute the filter step described above with this function
prepro_ht12 <- filterTechnicallyFailed(prepro_ht12, paramfile = myparamfile)

# show command execution history
kable(prepro_ht12$history)

## ----pre.pro.step.6-----------------------------------------------------------
# filter for batch-size
prepro_ht12 <- filter4MinBatchsize(prepro_ht12, paramfile = myparamfile)

# repeat transformation and normalization
# This is only necessary if samples were removed in previous step
prepro_ht12 <- transformNormalizeHT12object(prepro_ht12, paramfile = myparamfile)

# execute batch- adjustment via sva::ComBat()
prepro_ht12 <- removeBatchEffects(prepro_ht12, paramfile = myparamfile)

# review command history
kable(prepro_ht12$history)

## ----pre.pro.step.7-----------------------------------------------------------------------------------------------------------------------------------------------------
# check for remaining batch-effects
prepro_ht12 <- checkBatchEffects(prepro_ht12, paramfile = myparamfile)

# preview the newly created object
ht(prepro_ht12$genesdetail, 2)

# review command history
kable(prepro_ht12$history)

## ----pre.pro.step.8-----------------------------------------------------------------------------------------------------------------------------------------------------
# filter for atypical expression levels
prepro_ht12 <- filterAtypicalExpressed(prepro_ht12, paramfile = myparamfile)

# preview newly created output
ht(prepro_ht12$chipsamples, 2)

# review command history
kable(prepro_ht12$history)

## ----pre.pro.step.9-----------------------------------------------------------------------------------------------------------------------------------------------------
# create QC-plots
prepro_ht12 <- visualizePreprocessing(prepro_ht12, paramfile = myparamfile)

# PCA before/after preprocessing
prepro_ht12 <- comparePCBeforeAfterPrePro(prepro_ht12, paramfile = myparamfile)

# review command history
kable(prepro_ht12$history)

## ----pre.pro.step.10----------------------------------------------------------------------------------------------------------------------------------------------------
## Set parameter to use provided renaming file
renaming_fn <- system.file("extdata", "renamingsamples.txt", package = "HT12ProcessoR")

# edit parameter file with location of renaming-file
myparams[ myparams$variable == "file_renaming_samples_tosent", "value"] <- renaming_fn

# write the parameter file
write.table(myparams,
            myparamfile,
            quote = F,
            col.names = T,
            row.names = F,
            sep="\t")

# write result files into the tosend/ folder
prepro_ht12 = writeFilesTosent(prepro_ht12, paramfile = myparamfile)

# preview first and last line of $chipsamples object 
ht(prepro_ht12$chipsamples, 1)

# review command history
kable(prepro_ht12$history)

## ----save.rdata---------------------------------------------------------------------------------------------------------------------------------------------------------
# save an .RData object containing the final prepro_ht12 object
save(prepro_ht12, file = paste0(prepro.folder, "/obj/prepro_example1_HT12object.RData"))

## ----methods.text, results = "asis"-------------------------------------------------------------------------------------------------------------------------------------
# create text-object
mm_text <- createTextForMethods(prepro_ht12, paramfile = myparamfile)

# preview created text as a long version
cat('### Extended Version')
cat(mm_text$extended_version)

# preview created text as a short version
cat('### Short Version')
cat(mm_text$short_version)

## ----session.info, results = 'markup'-----------------------------------------------------------------------------------------------------------------------------------
# output script run time
total.time <- round(as.numeric(Sys.time() - time.start), 2)
message("This script ran for a total of ", total.time, " minutes.")

# output session info
sessionInfo()

