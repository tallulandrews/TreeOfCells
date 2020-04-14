\name{Input Checks}
\alias{check_row_col_names}
\title{Check a matrix had row and column names}
\description{
	Checks the input matrix has row and column names and if it doesn't, it creates some.
}
\usage{
	check_row_col_names(mat)
}
\arguments{
  \item{mat}{a matrix}
}
\details{
	Checks the input matrix has row and column names and if it doesn't, it creates some.
	Created row names are of the form: gene_#
	Created column names are of the form: cell_type_#

	This is primarily used for formating simulated data.
}
\value{
	the input matrix with new row or column namese if needed.
}
\examples{
	vals <- matrix(rnorm(70), ncol=7)
	mat <- check_row_col_names(vals);
}
