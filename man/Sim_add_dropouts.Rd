\name{add_dropouts}
\alias{add_dropouts}
\title{Adds dropouts to a vector of simulated counts}
\description{
	Adds dropouts to a vector of simulated counts in accordance with the M3Drop Michaelis-Menten dropout model.
}
\usage{
	add_dropouts(x,mu,K)
}
\arguments{
  \item{x}{a vector of simulated read/UMI counts.}
  \item{mu}{known mean expression for this gene.}
  \item{K}{Michaelis-Menten parameter.}
}
\details{
	Uses the Michaelis-Menten equation to add dropouts to a vector of countsrelative to the mean expression of the gene.
	\deqn{P_{drop} = 1 - \frac{mu}{mu + K}}{P_drop = 1 - mu /(mu + K)}
	
	This is useful for generating zero-inflated negative binomial simulations that resemble Smartseq2 single-cell RNA-seq data.
}
\value{
	vector of counts with added dropouts.
}
\examples{
counts <- rnbinom(100, size=0.5, mu=200)
counts <- add_dropouts(counts, 200, 50)
}
\keyword{negative binomial, dropouts, zero-inflation, simulation}
