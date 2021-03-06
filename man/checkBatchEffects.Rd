% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkBatchEffects.R
\name{checkBatchEffects}
\alias{checkBatchEffects}
\title{Identify remaining batch effects in expression data}
\usage{
checkBatchEffects(
  ht12object,
  paramfile = NULL,
  showPlots = T,
  round4ANOVAcheck = 5
)
}
\arguments{
\item{ht12object}{A list object of class HT12prepro created with function removeBatchEffects()}

\item{paramfile}{Path to the file specifying parameters}

\item{round4ANOVAcheck}{parameter for control how identical Results of MATRIXEQTL and standard ANova have to be. Provided number rounds the F statistic value. Applied in function runMAtrixEQTLAnova called by calcAnovaSentrixRunSpecialbatchViaMatrixEQTL2() . Note that slightly different results may occur  when e.g. using two completely nested covariates}
}
\value{
A list object of class HT12prepro where  the slot with  probe-related attributes of the current processing-stage named `$genesdetail` is updated as well as the  slot with the history of the commands named `$history`.  QQ plots of the association are shown.
}
\description{
Identification of transcripts where association with batches is still stronger than expected by chance after batch-adjustment. Stronger than chance (i.e. overinflation) is defined as association stronger p-value after Bonferroni correction. All defined Batches will be tested, i.e. `Sentrix.Barcode`,  `processingbatch`, and `strangebatch` in the table from slot `$chipsamples` if more than one category of the batch is found. Probes are not filtered but annotated when found to be overinflated. This function makes use of package MatrixEQTL in order to allow fast ANOVA on matrices. Therefore, some messages from these functions report doing an "eQTL analysis" allthough in fact an ANOVA is done on batch effects
}
