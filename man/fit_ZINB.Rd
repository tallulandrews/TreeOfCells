\name{fit_ZINB_to_matrix}
\alias{fit_ZINB_to_matrix}
\alias{fit_NB_to_matrix}
\alias{fit_zero_inflated_negative_binomial}
\title{Calculate Summary Statistics}
\description{
Calculates summary statistics to facilitate the fitting of a ZINB model, performs various data input checks.
}
\usage{
	fit_ZINB_to_matrix(counts)
	fit_NB_to_matrix(counts)
	fit_zero_inflated_negative_binomial(obs, g_row, vals, e=0.00001)
	
}
\arguments{
  \item{counts}{a numeric matrix of raw UMI or read counts, columns = samples, rows = genes.}
  \item{obs}{a single row from a count matrix}
  \item{g_row}{index/row number of the "obs" row}
  \item{vals}{summary statistics from \code{calculate_summary_values}}
  \item{e}{acceptable error for convergence of the EM algorithm}
}
\details{
\code{fit_ZINB_to_matrix} calculates summary statistics and runs \code{fit_zero_inflated_negative_binomial} on each row of the count matrix.

\code{fit_NB_to_matrix} calculates summary statistics and fits a library size adjusted negative binomial to each row of the count matrix.

\code{fit_zero_inflated_negative_binomial} fit a library-size adjusted zero-inflated negative binomial model to the expression of one gene using an expectation-maximization (EM) algorithm.

Library size adjusted zero-inflated negative binomial model:
\deqn{X~NBinom(mu_{ij}, r_j)*(1-d_j)+0*d_j}{X ~ NBinom(muij, rj) *( 1-dj) + 0*dj}

mu_ij = li*muj = mean expression of gene j in cell i
li = relative library size of cell i
muj = mean expression of gene j
rj = dispersion (size parameter of R's nbinom functions)
dj = dropout rate of gene j


}
\value{
\code{fit_zero_inflated_negative_binomial} a vector of mu, r, and d of the ZINB.
\code{fit_ZINB_to_matrix} a table of mu, r, and d for each gene in the matrix.
\code{fit_NB_to_matrix} a table of mu and r for each gene in the matrix.
}
\examples{
	counts <- matrix(rnbinom(10000, mean=10, size=0.5), ncol=20);
	fits <- fit_NB_to_matrix(counts);

	z_inflate <- matrix(runif(10000), ncol=20) < 0.2
	counts[z_inflate] <- 0;
	fits <- fit_ZINB_to_matrix(counts);

	vals <- calculate_summary_values(counts);
	fit <- fit_zero_inflated_negative_binomial(counts[3,], 3, vals)
}
\keyword{single cell, negative binomial}
\seealso{calculate_summary_values}
