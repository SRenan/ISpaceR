\docType{methods}
\name{read_exprs_mat}
\alias{read_exprs_mat}
\title{Read expression matrix}
\usage{
  read_exprs_mat(e_m_a)
}
\arguments{
  \item{e_m_a}{A \code{character}. The accession number of
  the matrix or matrices to read.}
}
\value{
  A single \code{ExpressionSet} if a single
  expression_matrix_accession has been given. A \code{list}
  of \code{ExprssionSet} otherwise, with one element for
  each given matrix.
}
\description{
  Reads the gene_expression_matrix table and creates an
  ExpressionSet for the selected matrix
}
\note{
  This function uses \code{RCurl} to query ImmuneSpace.
  Therefore, a valid .netrc file is needed when executed
  from a remote location.
}
\examples{
labkey.url.base <- "https://www.immunespace.org/"
labkey.url.path <- "/Yale/SDY61/"
read_exprs_mat("EXPM0001")
}
\seealso{
  \code{\link{limma_standard}} \code{\link{write_gea}}
}

