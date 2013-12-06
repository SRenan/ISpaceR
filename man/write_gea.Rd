\name{write_gea}
\alias{write_gea}
\title{Write gene expression analysis result to the database}
\arguments{
  \item{topTable}{A \code{data.frame}, as returned by the
  \code{limma_standard} function, or a \code{list} of
  \code{data.frame}.}

  \item{analysis_accession}{A \code{character} vector. The
  accession associated with the gene expression analysis.
  The vector should be of the same length as topTable if
  topTable is a list.}

  \item{description}{A \code{character} vector. The
  description of the analysis.}
}
\value{
  A \code{boolean}. TRUE if eveything went as expected.
}
\description{
  This function writes the output of a standard limma
  analysis to the database and fill the appropriate tables.
}
\seealso{
  \code{\link{read_exprs_mat}} \code{\link{limma_standard}}
  \code{\link{labkey.importRows}}
}

