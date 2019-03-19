\name{child.means}
\alias{child.means}
\title{Calculate Child Mean Expression}
\description{
	Calculates the mean expression of each gene in a child cell-population from the respective parental means and a relative distance between them.
}
\usage{
	child.means(parent.means, relative.dist, sigma=0.5, prop_genes=0.2)
}
\arguments{
  \item{parent.means}{a vector of mean gene expression values from parental cell-type.}
  \item{relative.dist}{a scalar value of the relative distance between child and parent cell populations.}
  \item{sigma}{standard deviation of normal distribution of log2-fold changes.}
  \item{prop_genes}{proportion of genes differentially expressed between parent/child.}
}
\details{
	This function calculates mean expression of a child cell-type from a progenitor (parent) cell-type as follows:
	log2(child) = log2(parent) + dist*rnorm(mean=0, sd=sigma)
	
	genes are selected at random with probability equal to \code{prop_genes} to be adjusted in this way.
}
\value{
	vector of mean expression values for the child cell-population.
}
\examples{
	parent <- 10^rnorm(100, mean=0, sd=2);
	child <- child.means(parent, 1)
	child <- child.means(parent, 1, sigma=1, prop_genes=0.9)
}
\keyword{simulation}
