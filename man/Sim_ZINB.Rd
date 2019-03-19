\name{ZINB.sim}
\alias{ZINB.sim}
\title{Simulate counts from a zero-inflated negative binomial.}
\description{
	Simulates a matrix of counts from a zero-inflated negative binomial, and provided parameters.
}
\usage{
	ZINB.sim(gene_means, lib.size, mean.size.fun, K=NULL)
}
\arguments{
  \item{gene_means}{vector of mean expression values for each gene.}
  \item{lib.size}{vector of relative library sizes for each cell.}
  \item{mean.size.fun}{a function that generates a the "size" parameter of the negative binomial distribution (1/dispersion) for a given mean expression level.}
  \item{K}{Michaelis-Menten parameter for zero-inflation.}
}
\details{
	1) Multiplies the gene_means and lib.size to get a matrix of gene x cell means for the negative binomial. 
	2) Generates gene-specific size parameters using the provided function
	3) Simulates counts from a negative binomial for each gene x cell entry in the matrix
	4) Adds dropout to each row according to a Michaelis-Menten equation as specified by K using \code{add_dropouts}.
}
\value{
	gene x cell matrix of simulated counts.
}
\examples{
constant_disp <- function(x) {return(1/0.2)}
g_means <- runif(10, min=1, max=100)
c_size <- runif(5)
zinb_counts <- ZINB.sim(g_means, c_size, constant_disp, K=50)
nb_counts <- ZINB.sim(g_means, c_size, constant_disp, K=NULL)
}
\keyword{negative binomial, dropouts, zero-inflation, simulation}
