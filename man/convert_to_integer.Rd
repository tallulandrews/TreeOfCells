\name{convert_to_integer}
\alias{convert_to_integer}
\title{Turn a matrix of expression values into integer counts}
\description{
	Reformats a provided expression matrix to be compatible with NBumi modelling. 
}
\usage{
	convert_to_integer(mat)
}
\arguments{
  \item{mat}{a numeric matrix of expression values (ideally raw UMI counts or raw read counts), columns = samples, rows = genes.}
}
\details{
	Coerces the provided data to a matrix then rounds all values up (ceiling) to integers and removes all rows where all values are zero.	
}
\value{
	Rounded, integer matrix of the original data.
}
\examples{
mat <- matrix(rgamma(1000, shape=10, scale=20), ncol=10)
mat_int <- convert_to_integer(mat)
is.integer(mat_int)
}
\keyword{negative binomial}
