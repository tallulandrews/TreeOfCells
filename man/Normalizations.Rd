\name{Normalization}
\alias{quantile_norm}
\alias{scnorm}
\alias{scaled_regression}
\alias{sctrans}
\alias{log_transform}
\alias{total_norm}
\alias{gene_length_norm}
\title{Profile Normalization Methods}
\description{
Methods for normalizing cell-type profiles.
}
\usage{
   quantile_norm(merged_profiles, conditions=NULL)
   scnorm(merged_profiles, conditions=rep(1, ncol(merged_profiles$mus)))
   scaled_regression(merged_profiles, conditions=NULL)
   sctrans(merged_profiles, conditions=NULL)
   log_transform(merged_profiles, conditions=NULL)
   total_norm(merged_profiles, min_rs = 10^-10, conditions=NULL)
   gene_length_norm(merged_profiles, gene_lengths, conditions=NULL)
}
\arguments{
  \item{merged_profiles}{Output from: merge_profiles function.}
  \item{conditions}{vector of different condition (only used in scnorm).}
  \item{gene_lengths}{a data.frame with two columns: symbol, length. The former is the gene IDs that matches your cell-type profiles.}
  \item{min_rs}{scalar value, minimum allowable value for the dispersion parameter "r"}
}
\details{
	Normalizes the various parameters of the cell-type profiles using different methods:

	scaled_regression: uses the \code{nlme} package, to perform weighted linear regression to eliminate the effect of total expression from each gene independently. Then scales each gene to have a mean of 0 and sd of 1.

	quantile_norm: quantile normalizes all parameters using \code{preprocessCore::normalize.quantiles}
	scnorm: normalizes all parameters using \code{SCnorm}.
	sctrans: normalizes means using \code{sctransform}.
	log_transform: log (base e) transforms both means and dispersion parameters.
	total_norm: normalizes both means and dispersion parameters using the median total expression.
	gene_length_norm: rescales each gene by the provided gene lengths.
}
\value{
	merged_profiles object with corrected matrices.
}
\examples{
obj <- list(
	mus = matrix(rnorm(200), ncol=10),
	rs  = matrix(rgamma(200, shape=0.75, scale=1), ncol=10),
	ds  = matrix(runif(200), ncol=10),
	Ns  = matrix(rpois(200, lambda=30), ncol=10))
colnames(obj$mus) <- paste("celltype", 1:ncol(obj$mus), sep="")
rownames(obj$mus) <- paste("gene", 1:nrow(obj$mus), sep="")
n_obj1 <- quantile_norm(obj)
n_obj2 <- scaled_regression(obj)
n_obj3 <- total_norm(obj)
n_obj4 <- log_transform(obj)
n_obj5 <- gene_length_norm(obj, gene_lengths=list(symbol=rownames(obj$mus), length=1:20))
}
\keyword{distance}
