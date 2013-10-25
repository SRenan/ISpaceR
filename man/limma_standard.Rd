\docType{methods}
\name{limma_standard}
\alias{limma_standard}
\title{Standard limma analysis}
\usage{
  limma_standard(eset, contrast = "study_time_reported",
    FDRthresh = 0.1)
}
\arguments{
  \item{eset}{An \code{ExpressionSet}. The expression set
  to analyse.}

  \item{contrast}{A \code{character}. The characteristic
  used to make contrasts. It should be a valid column name
  in the phenoData of eset.}

  \item{FDRthresh}{A \code{numeric}. The threshold used for
  filtering the False Discovery Rate.}
}
\value{
  A named \code{list} of \code{data.frame} of the same
  length as the number of contrasts.
}
\description{
  A wrapper around limma functions to do a standard gene
  expression analysis. Fit a linear model for each gene,
  computes some statistics and returns the result in a
  table (or list of tables). By default the function tries
  a time course experiment.
}
\seealso{
  \code{\link{read_exprs_mat}} \code{\link{write_gea}}
  \code{\link{lmFit}} \code{\link{topTable}}
}

