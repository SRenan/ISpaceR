\docType{methods}
\name{read_exprs_mat}
\alias{read_exprs_mat}
\title{Read expression matrix}
\usage{
  read_exprs_mat(expression_mat_accession)
}
\arguments{
  \item{expression_mat_accession}{A \code{character}. The
  accession number of the matrix or matrices to read.}
}
\value{
  An \code{ExpressionSet} for each matrix.
}
\description{
  Reads the gene_expression_matrix table and creates an
  ExpressionSet for the selected matrix
}
\examples{
labkey.url.base <- "https://www.immunespace.org/"
labkey.url.path <- "/Yale/SDY61/"
read_exprs_mat("EXPM0001")
}
\seealso{
  \code{\link{limma_standard}} \code{\link{write_gea}}
}

