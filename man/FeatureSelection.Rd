\name{Feature Selection}
\alias{fold_change_fs}
\alias{distribution_olap_fs}
\alias{dropout_fs}
\alias{top_diff_fs}
\title{Profile Feature Selection Methods}
\description{
Methods for performing feature selection on cell-type profiles.
}
\usage{
	top_diff_fs(merged_profiles, ntop=3000)
	dropout_fs(merged_profiles, diff=0.75)
	distribution_olap_fs(merged_profiles, quantile=0.5)
	fold_change_fs(merged_profiles, diff_t=2, max_t=1)
}
\arguments{
  \item{merged_profiles}{Output from: merge_profiles function.}
  \item{ntop}{number of features to select.}
  \item{diff}{difference in proportion of zeros}
  \item{quantile}{quantile threshold on overlap}
  \item{diff_t}{log2 fold-change minimum threshold.}
  \item{max_t}{lower threshold on maximum expression level.}
}
\details{
	Performs feature selection on the cell-type profiles using different methods:

	top_diff_fs : select genes with the biggest difference in mean expression across cell-types.
	dropout_fs : select genes with significantly higher difference in dropout proportion than \code{diff} (requires unscaled matrices)
	distribution_olap_fs : the distributions of the cell-type with the lowest expression and highest expression must by no more than proportion than \code{quantile}. (requires unscaled matrices)
	fold_change_fs : finds all genes passing a log2 fold-change threshold and expressed at at least max_t level in one cell-type.
}
\value{
	a vector of gene names of the selected features.
}
\examples{
obj <- list(
	mus = matrix(rgamma(200, shape=2.5, scale=2), ncol=10),
	rs  = matrix(rgamma(200, shape=0.75, scale=1), ncol=10),
	ds  = matrix(runif(200), ncol=10),
	Ns  = matrix(rpois(200, lambda=30), ncol=10))
fs1 <- top_diff_fs(obj, 10)
fs2 <- dropout_fs(obj)
fs3 <- distribution_olap_fs(obj)
fs4 <- fold_change_fs(obj)
}
\keyword{feature selection}
