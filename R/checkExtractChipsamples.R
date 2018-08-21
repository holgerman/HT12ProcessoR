#' @title Check Illumina HT12v4 expression data
#' @description This function checks the data and extracts the sample-IDs including information from Illuminas sample-files. See vignette for an example.
#' @param paramfile Path to the file specifying preprocessing parameters, can be NULL
#' @param file_names_of_files_and_folders Path to a tab-delimeted file with information of files created in Illuminas Genome Studio.  Following columns are required:   `fileset_id`: a short ID for a fileset,  `con_f`: full path and filename of the control probe data file,  `nobkgd_f`: full path and filename of the gene probe data file,  `sample_f`: full path and filename of the control probe data file. If NULL, information is used from parameterfile specified in 'paramfile', Default: 'from_paramfile'
#' @param colseparator Separator splitting numbers into units of 1000. Typically "." or  ",". Default: 'from_paramfile'
#' @param renameDublettes If any sample Ids is found more than once in the sample file provided by Illumina, suffixes .1 .2 etc are added to the newly created new_ID used to designate samples within the preprocessing. If FALSE, the function will stop when duplicated Sample Ids are found. Default: 'from_paramfile'
#' @param prefix_new_samplename Prefix of new_ID, the new R-compatibel sample names used to designate samples within the preprocessing. Default: 'from_paramfile'
#' @return a list of class `HT12prepro` including a slot with sample-related attributes of the current processing-stage named `$chipsamples`

#' @importFrom plyr rename
#' @importFrom stringr str_locate str_length str_sub str_split str_trim str_replace_all str_replace
#' @import data.table
#' @export


# debug
# paramfile = myparamfile
# file_names_of_files_and_folders = "from_paramfile"; colseparator= "from_paramfile"; renameDublettes= "from_paramfile"; prefix_new_samplename= "from_paramfile"

checkExtractChipsamples = function(paramfile = NULL, file_names_of_files_and_folders = "from_paramfile",  colseparator= "from_paramfile", renameDublettes= "from_paramfile",prefix_new_samplename= "from_paramfile") {


### strings are imported as strings and not as factors
  options(stringsAsFactors=FALSE)

  myparameters = match.call()

## ----addfunct------------------------------------------------------------
showVennplots = F

getIDTag <- function (file_con, endstr, id_tag_at_theend = T, general_prefix) {
  #3/8/2012
  # to get ID tag from Knut filenames
  # endstr = string_controlfiles
  nummeri = ifelse(id_tag_at_theend, 2, 1)
  ende = stringr::str_locate(string=file_con, pattern=endstr)[,nummeri]
  laenge = stringr::str_length(file_con)
  if(id_tag_at_theend) file_con <- stringr::str_sub(string=file_con, ende+2, laenge)
  if(id_tag_at_theend ==F) file_con <- stringr::str_sub(string=file_con, stringr::str_locate(string=file_con, pattern=general_prefix)[,1], ende-2)
  file_con
}

getIDsFromKnut <- function(x) {
  #20/7/2012
  #funktion, die die IDs der PRoben aus den Dateien HT12v4 sausliest
  #Diese funktioniert nur f?r nobkgd und bkgd und control files
  #
  # 24.6.14 umgestellt, das suffix statt punkt der ID erkenner ist, damit punkt in der ID auch erfasst werden kann
  # 7.10.14 angepasst
  #    x = "20140626-PROGRESS I Wdhl 115_Sample Probe Profile_FinalReport.txt"

  #   x = j
  test <- readLines(x, n=1)
  test = ifelse(test == "[Header]", readLines(x, n=9)[9], test)
  test <- unlist(stringr::str_split(test, "\t"))
  idstring = "\\.AVG|\\.Detection|\\.BEAD|\\.Avg|\\.NARRAYS|\\.ARRAY_STDEV"
  test <- (grep(idstring, test, value=T))
  test <- sapply(stringr::str_split(test, idstring), "[",1)
  table(test)
  #checken, dass alle eintr?ge f?r proben gleich h?ufig da sind, e.g. immer nbead, standardabweichung etc
  if(length(table(table(grep("[0-9]", test, value=T))))!=1) stop("Fuer manche als Zahlen kodierte Eintraege nicht alle Spaltenwerte gegeben!")
  if(length(table(table(test)))!=1) stop("Fuer manche als automatisch mit punkt kodierte Eintraege nicht alle Spaltenwerte gegeben!")
  samples_incl <- unique(test)
  as.character(samples_incl)
}

getIDsFromKnutSamplefiles <- function(x,file_overview=file_overview) {
  #20/7/2012
  #funktion, die die IDs der PRoben aus den Dateien HT12v4 ausliest
  #Diese funktioniert nur fuer files  derklasse sampleinfo
  # x = na.exclude(unlist(file_overview[,c("sample_f")]))[1]

  # x = i
  fn = grep(x, file_overview$sample_f, value  = T, fixed = T)
  vortest <- readLines( fn, n=1)
  lines2skip = ifelse(vortest == "[Header]", 8,  0)
  test <- read.delim( fn, skip = lines2skip, as.is = T)
  samples_incl <- test$Sample.ID
  as.character(samples_incl)
}

getAllInfosFromKnutSamplefiles3 <- function(x, docheck=T, mydec = colseparator_used) {   #  x="r:/genstat/01_daten/1603_eqtl_fettgewebe_kinder/ge_daten/Illumina_HT12_FG_Kinder_Gesamtgewebe/all_in_1folder/SamplesTable_premiRNA.txt"
  # x = i
  #20/7/2012  28.5.13 ergaenyt, dass nichtvorhandene spalten aus my_columns ergaenzt und aufgefuellt mit NA werden
  #funktion, die die PRobenattribute aus den Dateien HT12v4 sausliest

  my_columns <- c(c("Index", "Sample.ID", "Sample_ID", "Pool_ID",
                    "Sample_Well", "Sample_Plate", "Sample.Group",
                    "Sentrix.Barcode" , "Sample.Section","Detected.Genes..0.01.",
                    "Detected.Genes..0.05.", "Signal.Average", "Signal.P05",
                    "Signal.P25", "Signal.P50", "Signal.P75", "Signal.P95",
                    "BIOTIN", "CY3_HYB", "HOUSEKEEPING", "LABELING",
                    "LOW_STRINGENCY_HYB", "NEGATIVE..background.", "Noise"))
  #testen, bei erstellen der formel
  #        names(test)[(!names(test)  %in% my_columns )]
  #        my_columns[!my_columns %in% names(test)]
  filename <- x
  vortest <- readLines(filename, n=1)
  lines2skip = ifelse(vortest == "[Header]", 8,  0)
  all_info <- read.delim(filename, skip = lines2skip, dec = mydec, as.is = T
  )
  # ergaenzen
  missingcol = setdiff(my_columns, names(all_info))
  if(length(missingcol) >0) message("\nin file", filename, " missing cols are: ", missingcol  )
  for(mi in missingcol){
    all_info[mi] = NA
  }
  all_info = all_info[,my_columns]
  if(docheck==T & identical(all_info$Sample_ID, all_info$Sample.ID)==F) {
    if(all(is.na(all_info$Sample_ID)==T)) message (paste("keine eintraege in spalte  Sample_ID zum vergleichen mit Spalte  Sample ID im info file von  file ", x))
    if(all(identical(all_info$Sample.ID, stringr::str_trim(all_info$Sample.ID) )==F)) warning("Leading/tailing spaces in 'Sample.ID' field found, ignoring them in consistency check")
    # if(all(is.na(all_info$Sample_ID)!=T)) stop (paste("spalte Sample_ID und
    #Sample ID matcht nicht im Sample-info file im file ", x))
  }
  if(docheck==F & identical(all_info$Sample_ID, stringr::str_trim(all_info$Sample.ID)==F)) {
    if(unique(all_info$Sample_ID) =="") message (paste("keine eintraege in spalte  Sample_ID zum vergleichen mit Spalte  Sample ID im info file von  file ", x))
    if(unique(all_info$Sample_ID) !="") message (paste("spalte Sample_ID und  Sample ID matcht nicht im Sample-info file im file ", x))
  }
  # venndiagram(x=all_info$Sample_ID, y=all_info$Sample.ID, unique=T, type="2")
  all_info$sample_f <- x
  all_info
}

checkKnutsBkgd <- function(filename_1,filename_2){
  #20/7/2012
  ##Anschauen, dass sich die Mittelwerte der Backgroundkorrigierten Proben auch unterscheiden von den nicht-backgroundkorrigierten Files
  ##dabei einlesen und Auswerten der ersten 100 Zeilen
  x=filename_1
  file1 <- readLines(x, n=101)
  file1 <- (stringr::str_split(file1, "\t"))
  file1 <- data.frame(file1,stringsAsFactors=F)
  file1 <- data.frame(t(file1),stringsAsFactors=F)
  names(file1) <- as.character(file1[1,])
  file1 <- file1[-1,]
  y=filename_2
  file2 <- readLines(y, n=101)
  file2 <- (stringr::str_split(file2, "\t"))
  file2 <- data.frame(file2,stringsAsFactors=F)
  file2 <- data.frame(t(file2),stringsAsFactors=F)
  names(file2) <- as.character(file2[1,])
  file2 <- file2[-1,]
  file1 <- file1[file1$PROBE_ID %in% intersect(file1$PROBE_ID,file2$PROBE_ID), ]
  file2 <- file2[file2$PROBE_ID %in% intersect(file1$PROBE_ID,file2$PROBE_ID), ]
  if(dim(file1)[1]<50) stop(paste(dim(file1)[1],"Vergleiche sind zu wenig vergleiche (Bedingung: mind 50 g?ltige von 100 potentiellen Vergleichen"))
  # print(paste(x,y,dim(file1)[1],identical(file1,file2)))
  return(identical(file1,file2))
}

## ----load----------------------------------------------------------------


#laden parameter
if(is.null(paramfile)==F) param <- read.delim(paramfile, as.is = T)
# str(param)
#beispiel aller ilmn Genexpresssionsids eines vollst. datensatzes aus  "r:/genstat/01_daten/1111_ht12v4_lifeb3_hss/Einzelanalyse_nonorm_nobgkd_HSS38-AMI457.txt"

laengeilmn = dim(allilmn)[1]
if(laengeilmn != 47323) stop("47323 ilmn-IDs not correctly loaded")
# ilmn controllannotation


## ----read----------------------------------------------------------------
if(file_names_of_files_and_folders== "from_paramfile") file_names_of_files_and_folders_used = getParam2("file_names_of_files_and_folders", myparam = param) else file_names_of_files_and_folders_used = file_names_of_files_and_folders

message("Using files as input found in ", file_names_of_files_and_folders_used)
file_overview = read.delim(file_names_of_files_and_folders_used, as.is = T)

showClassDF(file_overview)
## ----checkvol------------------------------------------------------------
checkvollstaendig <- file_overview[is.na(file_overview$sample)|is.na(file_overview$nobkgd)|is.na(file_overview$con),]
checkvollstaendig #keine NA filenames
if(dim(checkvollstaendig)[1]    != 0) stop("Problem found: Incomplete filesets - three files of type 'nobkgd', 'con' and 'sample' NOT allways found).\n")

status01_allelfilesda = "Allways complete filesets identified, i.e. three files of type 'nobkgd', 'con' and 'sample', no problem identified.\n";
message(status01_allelfilesda)

## ----FiletattributeBaue--------------------------------------------------
# Informationen aus  bkgd, nobkgd, control  auslesen
sample_overview <- c()
file_types <-  c("nobkgd_f", "con_f")
for (i in 1:length(file_types)){
  my_filetype <- file_types[i]
  for (j in na.omit(unlist(file_overview[,my_filetype]))) {
    # j = na.omit(unlist(file_overview[,my_filetype]))[1]
    j = stringr::str_replace_all(string=j, pattern='\\\\ ', replacement=" ") # workaround fuer utf-8 kodierung falls leerzeichen auftauchen
    message('reading ', j)
    # j = "SamplesTable_premiRNA.txt"
    samples <- (getIDsFromKnut(j))
    if(class(samples)=="try-error") samples = NA
    filenames <- j
    file_class <- i
    sample_overview_neu <- data.frame(samples,filenames, file_class, stringsAsFactors = F)
    sample_overview <- rbind(sample_overview, sample_overview_neu)
  }
}
sample_overview
showClassDF(sample_overview)
## ----addsample-----------------------------------------------------------
next_fileclass <- max(sample_overview$file_class)+1
for (i in na.exclude(unlist(file_overview[,c("sample_f")]))){
  #   i = na.exclude(unlist(file_overview[,c("sample_f")]))[1]
  # i =  "SamplesTable_premiRNA.txt"
  message("reading ", i)
  samples <- getIDsFromKnutSamplefiles(i,  file_overview = file_overview)
  filenames <- i
  file_class <- next_fileclass
  # print(paste(samples,filenames, file_class))
  sample_overview_neu <- data.frame(samples,filenames, file_class)
  sample_overview <- rbind(sample_overview, sample_overview_neu)
}
# cleanup
rm(sample_overview_neu)

## ----ergi----------------------------------------------------------------
dim(sample_overview)

## ----checi11-------------------------------------------------------------
checkalle = showNA(sample_overview)
checkalle
if(sum(checkalle$NAs) != 0) {
  status01b_allelPersoneninfosda = "\nProblem found: Some sample IDs not present in all three files of type 'nobkgd', 'con' and 'sample'\n"
  stop(status01b_allelPersoneninfosda)
} else  {status01b_allelPersoneninfosda = "\nAll sample IDs present in all three files of type 'nobkgd', 'con' and 'sample', no problem found.\n"; message(status01b_allelPersoneninfosda)}

## ----iterncode-----------------------------------------------------------
newprefix <-  param[param$variable =="prefix_new_samplename", "value"]
if(stringr::str_length(newprefix) < 1) stop("please provide `prefix_new_samplename` with a string length of at least 1")
newprefix
id_table <- data.frame(oriname=unique(sample_overview$samples))
id_table$genstatid <- paste0(newprefix, seq(along=id_table$oriname))
sample_overview$newID <- id_table[match(sample_overview$samples, id_table$oriname), "genstatid"]
ht(sample_overview)

## ----unique--------------------------------------------------------------
ht(sample_overview)
# Check ob proben doppelt vorkommen--> die wuerden bei reshape unter den Tisch fallen
# Hier muss man noch einen rename-wrapper bauen
check <-table(sample_overview$samples) #[sample_overview$file_classsamples)
table(check)

#welche betrifft das
mehrfach_probes <- table(sample_overview$samples)[as.numeric(check) > 3]
mehrfach_probes
status02_mehrfachIDs = ifelse(length(mehrfach_probes) ==0, "No duplicated IDs across files detected, no problem found.\n", "Problem found and resolved: Duplicated IDs across files detected, treated according parameter renameDublettes.\n")
status02_mehrfachIDs
ht(sample_overview[ sample_overview$samples %in% names(mehrfach_probes),])
sample_overview[ sample_overview$samples %in% names(mehrfach_probes)[1],]
if(any(as.numeric(names(table(check))) < 3)) stop("Probleme mit ID import  - checke Skript code sdlfjaklsdff")
if(any(as.numeric(names(table(check)))/3 != round(as.numeric(names(table(check)))/3,0))) stop("Probleme mit ID import  - checke Skript code sdlfjaklsdff2")

## ----doublette-----------------------------------------------------------
sample_overview$samplename_count <- check[match(sample_overview$samples, names(check))]/3
if(max(as.numeric(names(table(check))) > 3)) {
  message(paste0("CAVE - Proben mit Nummer(n): \n", names(mehrfach_probes),
                 "\n kommt/kommen doppelt vor, modifiziere code fuerr diesen Fall"))
  message("\n\n")
  message(paste0("betrifft \n", length(names(mehrfach_probes)), " von ", length(unique(sample_overview$samples)), " Proben"))
  message("\n\n")
  # print(sample_overview[ sample_overview$samples %in% names(mehrfach_probes),])


  if(renameDublettes== "from_paramfile") param2 = getParam2("renameDublettes", myparam = param) else param2 = renameDublettes

  if(param2 == "F" | param2 == "FALSE") stop("Found samples with duplicated name. As parameter 'renameDublettes'==F , stopping here. Consider setting parameter renameDublettes to TRUE)")
  if(param2 == "T" | param2 == "TRUE") {
    if(any(ceiling(sample_overview$samplename_count) != sample_overview$samplename_count)) stop("manche proben kommen nicht in allen files nobkgd, con und sample vor")
    sample_overview <- sample_overview[order(sample_overview$samples, sample_overview$file_class),]
    to_rename <- unique(sample_overview[ sample_overview$samplename_count>1 , "samples"])
    for(i in to_rename ){   #  i="4556239"
      message("Working on ", i)
      max_i <- unique(sample_overview[ sample_overview$samplename_count>1 , "samplename_count"])
      index <- 1
      for(j in 1:max_i){  #j=1
        sample_overview[sample_overview$samples == i  , "samplename_count"][c(j, j+max_i, j+max_i*2)] <- index
        index <- index+1
      }

    }
    stopifnot(all(is.na(sample_overview$samplename_count))==F)
    check <- table(table(paste(sample_overview$samples, sample_overview$samplename_count)))
    stopifnot(identical(names(check), "3"))
    #suffix .1 .2 anfuegen
    sample_overview[sample_overview$samples %in% to_rename, "newID"] <- paste0(
      sample_overview[sample_overview$samples %in% to_rename, "newID"], ".",                                                                  sample_overview[sample_overview$samples %in% to_rename, "samplename_count"])
    message("\n Renamed samples: \n")
    message(sample_overview[ sample_overview$samples %in% to_rename,])
  }
  else stop("Only  TRUE or FALSE is a valid entry in parameter renameDublettes")
}
ht(id_table)
showNA(id_table)

## ----wide----------------------------------------------------------------
showNA(sample_overview)
sample_overview_l <- reshape(sample_overview, idvar= c("samples", "newID",  "samplename_count"), timevar="file_class",direction="wide")
sample_overview_l$samplename_count <- NULL
head(sample_overview_l)
showNA(sample_overview_l)
sample_overview_l = plyr::rename(sample_overview_l, c(samples = "old_ID",  newID = "samples", filenames.1 ="nobkgd_f",filenames.2 ="con_f", filenames.3="sample_f"))
ht(sample_overview_l,2)
check <- table(table(sample_overview_l$samples))
check
if(identical(names(check), "1")==F) stop("code 12341.2 proben kommen immer noch doppelt vor, kann nicht sein")
showNA(sample_overview_l)
##Liste mit weiteren Probenannotationen hinzuf?gen
#unique identifyer bauen
sample_overview_l$fileset_id <- file_overview[ match_hk(sample_overview_l$sample_f, file_overview$sample_f),"fileset_id"]

## ----umbennennen---------------------------------------------------------
ht(sample_overview_l[ grep("_neu", sample_overview_l$fileset_id),])
sample_overview_l$fileset_id = stringr::str_replace( sample_overview_l$fileset_id, "_neu","")
ht(sample_overview_l,3)

## ----weiter4-------------------------------------------------------------
sample_overview_l$paste_id <- with(sample_overview_l, paste(old_ID, fileset_id, sep="."))
#weitere annotationen hinzuf?gen aus sample - files
#Da manche Dateien nicht alle Attribute enthalten,werden diese getrennt ausgelesen
#proben mit allen Attributen
#proben mit allen Spalten
head(sample_overview_l)
showNA(sample_overview_l)
all_infofiles <- unique(sample_overview_l$sample_f)
all_infofiles <- all_infofiles[!all_infofiles %in% grep("/NA$",all_infofiles, value=T )] #hier wird de facto nichts gefiltert
getwd()

if(colseparator== "from_paramfile") colseparator_used = getParam2("colseparator", myparam = param) else colseparator_used = colseparator


all_sample_infos <- c()
for (i in setdiff(all_infofiles,c(""))){

  # message(paste("processing...",i))
  try(sample_infos <- getAllInfosFromKnutSamplefiles3(i))
  try(all_sample_infos <- rbind(all_sample_infos, sample_infos))
}
showNA(all_sample_infos)
unique(all_sample_infos$sample_f)
dim(all_sample_infos)
ht(all_sample_infos)
rm(sample_infos)

## ------------------------------------------------------------------------
showNA(all_sample_infos)
na_sentrix_filter = is.na(all_sample_infos$Sentrix.Barcode)
if(any(na_sentrix_filter)) {
  all_sample_infos$sentr_help = sapply(stringr::str_split(all_sample_infos$Sample.ID, "_"), "[", 1)
  all_sample_infos$section_help = sapply(stringr::str_split(all_sample_infos$Sample.ID, "_"), "[", 2)
  all_sample_infos[ na_sentrix_filter, 'Sentrix.Barcode'] =  all_sample_infos[ na_sentrix_filter, 'sentr_help']
  all_sample_infos[ na_sentrix_filter, 'Sample.Section'] =  all_sample_infos[ na_sentrix_filter, 'section_help']
  showNA(all_sample_infos)
}

## ----proceeding----------------------------------------------------------
#Was steht in den  manchmal fehlenden Attributen drinne, wenn die da sind?
sapply(all_sample_infos[,c("Pool_ID", "Sample_Well", "Sample_Plate")], function (x) table(x, useNA = "always"))
if(identical(unique(unlist(all_sample_infos[,c("Pool_ID", "Sample_Well", "Sample_Plate")])), NA) == F) stop("code 39867632 unerwartete eintr?ge,wo NA erwartet")
all_sample_infos = all_sample_infos[ ,names(all_sample_infos) %nin% c("Pool_ID", "Sample_Well", "Sample_Plate")]
###Filter: ist der SEntrix ID richtig angegeben
bad_sentrix <- unique(all_sample_infos$Sentrix.Barcode)
bad_sentrix
# str(bad_sentrix)
check <- table(stringr::str_length(bad_sentrix))
# noch NA checken, wurde aber ausch schon in ExtraktionsFunktion gemacht
sm3 = all_sample_infos[ is.na(all_sample_infos$Sample_ID),]
ht(sm3,2)
qlist33 = venn2(all_sample_infos$Sample.ID, all_sample_infos$Sample_ID, plotte = showVennplots)
#identifizierung, wo der chip-identifyer falsch ist
bad_sentrix_file <- unique(all_sample_infos[stringr::str_length(all_sample_infos$Sentrix.Barcode)  %nin% c(10, 11),"sample_f"])
bad_sentrix_file
sm_badsentrix = all_sample_infos[stringr::str_length(all_sample_infos$Sentrix.Barcode) %nin% c(10, 11),]
sm_badsentrix
# ## ersetze haendisch sentrix id im glauben dass reihenfolge passt
# all_sample_infos[ all_sample_infos$sample_f ==  bad_sentrix_file,"Sentrix.Barcode"] = c(rep("misSentri1", 12), rep("misSentri2", 12), rep("misSentri3", 12), rep("misSentri3", 12))
# bad_sentrix_file <- unique(all_sample_infos[stringr::str_length(all_sample_infos$Sentrix.Barcode) != 10,"sample_f"])
# bad_sentrix_file
if(length(bad_sentrix_file) != 0) status03_badSentrix = paste0("Problem found: no 10 or 11-digit Sentrix ID found in file(s): \n", bad_sentrix_file) else status03_badSentrix ="Sentrix ID allways found with 10 or 11 digits, no problem identified."
if(length(bad_sentrix_file) != 0) warning(status03_badSentrix) else message(status03_badSentrix)

## ----ProcessInfosEtBAckroundcheck----------------------------------------
#zusammenmergen der infos zu den Dateien und zu den probenattributen
dim(sample_overview_l)
dim(all_sample_infos)
head(sample_overview_l)
head(all_sample_infos)
showNA(sample_overview_l)
#bauen eines identifyers aus filename und evtl. doppelten
matcher = unique(sample_overview_l[,c('sample_f',"fileset_id")])
all_sample_infos$fileset_id <- matcher[match_hk(all_sample_infos$sample_f,
                                                matcher$sample_f),"fileset_id"]
all_sample_infos$paste_id <- paste(all_sample_infos$Sample.ID, all_sample_infos$fileset_id, sep=".")
sample_overview_l <- merge(sample_overview_l, all_sample_infos[,!names(all_sample_infos) %in% c("sample_f","fileset_id")], by="paste_id", sort=F, all.x=T,  incomparables=c(NA, NaN))
ht(sample_overview_l)
rm(all_sample_infos)
showNA(sample_overview)
#auf inkorrekte nomenklatur pr?fen
allnames <-   names(sample_overview_l)
check <- allnames[grep("bgkd", allnames)]
if(length(check) != 0) stop("darf nicht sein code laklafklafljfalj")

## ----ShowExpressedGenes--------------------------------------------------
# exprimierte Gene visualisieren
par(mfrow=c(1,1))
titele = "Detected HT12-Probes pval > 0.05"
# hist(x=sample_overview_l$Detected.Genes..0.05., breaks=200, col="red", main = titele)
# boxplot(sample_overview_l$Detected.Genes..0.05., col="red", main=titele)

## ----allelIlmnDA---------------------------------------------------------
# pruefen, ob alle ilmn imputiert worden
# str(allilmn)
for(i in na.omit(unique(sample_overview_l$nobkgd_f))){
  #   i = na.omit(unique(sample_overview_l$nobkgd_f))[1]
  message("\ncheck following file for presence of 47,323 ilmn-probes:      ", i, "\n")
  pfad <- i
  fileset_id <- unique(na.omit(sample_overview_l[ sample_overview_l$nobkgd_f ==i, "fileset_id"]))
  ilmn_read <- fread(pfad, select = "PROBE_ID")
  ilmn_read = ilmn_read$PROBE_ID
  # print(sum(ilmn_read %in% allilmn$ilmn))
  allilmn$neu <- allilmn$ilmn %in% ilmn_read
  names(allilmn)[names(allilmn)=="neu"] <- fileset_id
}

# str(allilmn)
hh(allilmn,2)
uebersicht_ilmns <- sapply(allilmn[,2:ncol(allilmn), drop = F], function(x) sum(x))
checkalleilmndrinn <- table(uebersicht_ilmns == laengeilmn, useNA="ifany")
checkalleilmndrinn
uebersicht_ilmns[ uebersicht_ilmns != laengeilmn]
missing_ilmns <- sapply(allilmn[,2:ncol(allilmn), drop = F], function(x) paste(allilmn$ilmn[x ==F], collapse = ", "))
missing_ilmns = missing_ilmns[ stringr::str_length(missing_ilmns) !=0]
missing_ilmns = data.frame(missing_ilmns=missing_ilmns, fileset_id = names(missing_ilmns))
all(missing_ilmns$fileset_id %in% sample_overview_l$fileset_id)
missing_ilmns$nobkgd_f = file_overview[ match_hk(missing_ilmns$fileset_id, file_overview$fileset_id), "nobkgd_f"]
missing_ilmns$n_ilmn = uebersicht_ilmns[ match_hk(missing_ilmns$fileset_id, names(uebersicht_ilmns))]
missing_ilmns

# nachschauen, wass das ist
missing_ilmns_IDs = stringr::str_trim(unlist(stringr::str_split(missing_ilmns$missing_ilmns, ", ")))
relcols = c("ilmn", "nuid", "Source", "ILMN_Gene",
            "Entrez_Gene_ID", "Accession", "Symbol",
            "Probe_Type", "Cytoband", "Definition",
            "Synonyms", "is_control", "EntrezReannotated",
            "ProbeQuality",  "ProbeSequence", "SymbolReannotated",
            "SymbolReannotated_orgHsEg", "GenomicLocation_chr",
            "GenomicLocation_start", "GenomicLocation_ende",
            "badprobe_dunning")

if(length(missing_ilmns_IDs)>0) missingilmns_det = ilmnAnnot014allgInfos[ paste0("ILMN_", ilmnAnnot014allgInfos$ilmn) %in% missing_ilmns_IDs,relcols] else missingilmns_det = structure(list(ilmn = integer(0), nuid = character(0), Source = character(0),                                                                                      ILMN_Gene = character(0), Entrez_Gene_ID = integer(0), Accession = character(0),                                                                                   Symbol = character(0), Probe_Type = character(0), Cytoband = character(0),                                                                                  Definition = character(0), Synonyms = character(0), is_control = logical(0),                                                                                EntrezReannotated = integer(0), ProbeQuality = character(0),                                                                                ProbeSequence = character(0), SymbolReannotated = character(0),                                                                                                SymbolReannotated_orgHsEg = character(0), GenomicLocation_chr = character(0),                                                                                          GenomicLocation_start = integer(0), GenomicLocation_ende = integer(0),                                                                                      badprobe_dunning = logical(0)), .Names = c("ilmn", "nuid", "Source", "ILMN_Gene", "Entrez_Gene_ID", "Accession", "Symbol", "Probe_Type", "Cytoband", "Definition", "Synonyms", "is_control", "EntrezReannotated", "ProbeQuality", "ProbeSequence", "SymbolReannotated", "SymbolReannotated_orgHsEg", "GenomicLocation_chr", "GenomicLocation_start", "GenomicLocation_ende", "badprobe_dunning"), row.names = integer(0), class = "data.frame")



status04_missingilmns = ifelse(length(checkalleilmndrinn)!= 1, paste("Problem found: Following expression files did not include all",
                                                                     laengeilmn, "Illumina probes:\n",
                                                                     paste0(missing_ilmns$nobkgd_f, "\n")),
                               paste("All expression files included ", laengeilmn,
                                     "Illumina probes, no problem identified.\n"))
message(status04_missingilmns)

## ----allelIlmnDA2--------------------------------------------------------
head(all_con)

# all_con = read.delim("_archive//illumina_annot_control_beads_110407hk.txt")

for(i in na.omit(unique(sample_overview_l$con_f))){
  # i = "ControlProbeProfile_nonorm_nobkgd_A1-Kohorte.txt"
  message("Check following file for presence of 883 Illumina-control probes:      ", i, "\n")
  pfad <-i
  fileset_id <- unique(na.omit(sample_overview_l[ sample_overview_l$con_f ==i, "fileset_id"]))
  ilmn_read <- fread(pfad, select = 1:2)
  # str(ilmn_read)
  ilmn_read = ilmn_read$ProbeID
  # print(sum(ilmn_read %in% all_con$Array_Address_Id))
  all_con$neu <- all_con$Array_Address_Id %in% ilmn_read
  names(all_con)[names(all_con)=="neu"] <- fileset_id
}

ht(all_con,1)
uebersicht_ilmnCons <- sapply(all_con[ duplicated(all_con$Probe_Id)==F,7:ncol(all_con), drop = F], function(x) sum(x))
checkall_condrinn <- table(uebersicht_ilmnCons == length(all_con$Array_Address_Id), useNA="ifany")
checkall_condrinn
uebersicht_ilmnCons[ uebersicht_ilmnCons != length(all_con$Array_Address_Id)]
missing_ilmnCons <- sapply(all_con[,7:ncol(all_con), drop = F], function(x) paste(all_con$Probe_Id[x ==F], collapse = ", "))
missing_ilmnCons = missing_ilmnCons[ stringr::str_length(missing_ilmnCons) !=0]
missing_ilmnCons = data.frame(missing_ilmnCons=missing_ilmnCons, fileset_id = names(missing_ilmnCons))
all(missing_ilmnCons$fileset_id %in% sample_overview_l$fileset_id)
missing_ilmnCons$nobkgd_f = file_overview[ match_hk(missing_ilmnCons$fileset_id, file_overview$nobkgd), "nobkgd_f"]
missing_ilmnCons$n_ilmn = uebersicht_ilmnCons[ match_hk(missing_ilmnCons$fileset_id, names(uebersicht_ilmnCons))]
laengeilmncon = unique(uebersicht_ilmnCons)
status04_missingilmnCons = ifelse(length(checkall_condrinn)!= 1,
                                  paste("Problem found: Following control files did not include all",
                                        length(all_con$Array_Address_Id), "Illumina control probes:\n",
                                        paste0(missing_ilmnCons$nobkgd_f, "\n")),
                                  paste("All control files included all",
                                        laengeilmncon, "Illumina control probes, no problem identified.\n"))


## ----anhuebschen---------------------------------------------------------
# gleiche Ids
sample_overview_l$Sample.IDclean = stringr::str_trim(sample_overview_l$Sample.ID) # oben, in funktion getIDsFromKnutSamplefiles wird ein warning generiert fuer diesen fall
sample_overview_l$old_ID.IDclean = stringr::str_trim(sample_overview_l$old_ID) # oben, in funktion  getIDsFromKnutSamplefiles wird ein warning generiert fuer diesen fall
checkedIDs = apply(sample_overview_l[, c('old_ID.IDclean', "Sample.IDclean", "Sample_ID")], 1,function (x) length(na.omit(unique(stringr::str_trim(x)))) >1)
check = sample_overview_l[checkedIDs,c(  'old_ID.IDclean', "Sample.IDclean", "Sample_ID")]
check
checkdetail = sample_overview_l[checkedIDs,]
status05checkedIDsConsistency = ifelse(nrow(check) != 0,('Problem found: Different IDs for identical samples used within columns "Sample.ID", "Sample_ID", and rownames of Expressionfiles. Please clarify reasons!'),('No different IDs for identical samples used within columns "Sample.ID", "Sample_ID", and rownames of Expressionfiles, no problem identified.\n'))
status05checkedIDsConsistency_data = check
if(nrow(check) != 0) warning(status05checkedIDsConsistency) else message(status05checkedIDsConsistency)
sample_overview_l$Sample.IDclean = NULL
sample_overview_l$old_ID.IDclean = NULL

## add on 25.6.18 old ID darf trailing leerzeichen haben, wird implizit ignoriert durch str_trim. Die trailing spacesim expression file scheinen ohnehin ignoriert zu werden,also sogar notwendig, wenn einheitlich trailing spaces vergeben wurden
sample_overview_l$old_ID = stringr::str_trim(sample_overview_l$old_ID)


## ----releva--------------------------------------------------------------
names(sample_overview_l)[names(sample_overview_l) == "samples"] = "new_ID"
good_columns = c("paste_id", "old_ID",  "new_ID",  "nobkgd_f", "con_f", "sample_f", "fileset_id", "Index", "Sample.Group", "Sentrix.Barcode", "Sample.Section", "Detected.Genes..0.01.", "Detected.Genes..0.05.", "Signal.Average", "Signal.P05", "Signal.P25", "Signal.P50", "Signal.P75", "Signal.P95", "BIOTIN", "CY3_HYB", "HOUSEKEEPING", "LABELING", "LOW_STRINGENCY_HYB", "NEGATIVE..background.", "Noise")  # , "bkdg_prblm" entfernt, weil nicht mehr daten dafuer da
qlist44 = venn2(good_columns, names(sample_overview_l), plotte = showVennplots)
# str(qlist44)
sample_overview_l2 = sample_overview_l[,good_columns]

## ----sentrihaeurf und umwandeln--------------------------------------------------------
sample_overview_l2$Sentrix.Barcode = as.character(sample_overview_l2$Sentrix.Barcode) # neu 21.8.18 to avoid int64 class
check_haeufigkeitSentrixIDs = table(table(sample_overview_l2$Sentrix.Barcode))
check_haeufigkeitSentrixIDs
if(as.numeric(max(names(check_haeufigkeitSentrixIDs), na.rm = T)) >12) warning("sentrix IDs should only have 12 indis")

## ----speichern-----------------------------------------------------------


# for doku
alle_statusse = ls()[grep("^status", ls())]
alle_statusse
status01_allelfilesda
status02_mehrfachIDs
status03_badSentrix
status05checkedIDsConsistency

qlist332 = venn2(sample_overview$samples, sample_overview_l2$old_ID, plotte = showVennplots)


sample_overview_l2$in_study = TRUE
sample_overview_l2$reason4exclusion = NA_character_
sample_overview_l2$subgroup = 'all'
sample_overview_l2$processingbatch = sample_overview_l2$fileset_id
sample_overview_l2$strangebatch = 0

res = vector(mode = "list", length = 1)

class(res) = unique(c(class(res), "HT12prepro"))

res$chipsamples = sample_overview_l2
res$dokuobjects_checkExtractChipsamples = list(missing_ilmns, missing_ilmnCons,missingilmns_det,sample_overview_l2, check_haeufigkeitSentrixIDs, mget(alle_statusse))
names(res$dokuobjects_checkExtractChipsamples)  = c('missing_ilmns', 'missing_ilmnCons','missingilmns_det','sample_overview_l2', 'check_haeufigkeitSentrixIDs', 'alle_statusse')
res$history = data.frame(calls = paste(Sys.time(), deparse(myparameters)))
return(res)

}
