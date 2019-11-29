\name{calculate_summary_values}
\alias{calculate_summary_values}
\title{Calculate Summary Statistics}
\description{
Calculates summary statistics to facilitate the fitting of a ZINB model, performs various data input checks.
}
\usage{
	calculate_summary_values(counts)
	
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI or read counts, columns = samples, rows = genes.}
}
\details{
Calculates various summary statistics: total counts by genes and cells, total dropouts by genes and cells and dimensions; which are necessary for fitting library-size adjusted negative binomial models or library-size adjusted zero-inflated negative binomial models.

Compatible with both regular matrices and sparse matrices as defined in the Matrix package.

}
\value{
	A named list of: 
		tis = total molecules per cell, 
		tjs = total molecules per gene, 
		dis = total dropouts per cell, 
		djs = total dropouts per genes, 
		total = total molecules in dataset, 
		nc = number of cells, 
		ng = number of genes
}
\examples{
	counts <- matrix(rnbinom(10000, mu=10, size=0.5), ncol=20);
	stats <- calculate_summary_values(counts);
}
\keyword{single cell, negative binomial}
