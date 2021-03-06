% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualizePreprocessing.R
\name{visualizePreprocessing}
\alias{visualizePreprocessing}
\title{Visualize expression data before and after preprocessing}
\usage{
visualizePreprocessing(
  ht12object,
  paramfile = NULL,
  show_qcplots_individSamplelevel = "from_paramfile",
  showPlots = T
)
}
\arguments{
\item{ht12object}{A list object of class HT12prepro created with function filterAtypicalExpressed()}

\item{paramfile}{Path to the file specifying parameters}

\item{show_qcplots_individSamplelevel}{Show QC plots comparing expression values before / after preprocessing either binned on SentrixID ('FALSE') or binned on Sample ID ('TRUE')? If "from_paramfile", than the parameter will be read from the paramfile with the location of this file given in parameter paramfile.}

\item{showPlots}{Show plots. If FALSE, plots are only stored in slot $dokuobjects_visualizePreprocessing}
}
\value{
A list object of class HT12prepro including additional slots  $dokuobjects_visualizePreprocessing
}
\description{
Visualisation of control features and expression data before and after preprocessing in several plots. Thereby, raw data before preprocessing is log2 transformed for a more meaningful comparison
}
