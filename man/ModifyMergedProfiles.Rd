\name{apply_feature_selection}
\alias{apply_feature_selection}
\title{Apply Feature Selection}
\description{
	Subset all matrices of a merged_profiles object to contain only the specified features.
}
\usage{
	apply_feature_selection(merged_profiles, features)
}
\arguments{
  \item{merged_profiles}{Output from: merge_profiles function.}
  \item{features}{a vector of gene names.}
}
\details{
	apply_feature_selection subsets all the matrices of a merged_profiles object to contain only the specified features.
}
\value{
	a merged_profiles object.
}
\examples{
obj <- list(
	mus = matrix(rnorm(200), ncol=10),
	rs  = matrix(rgamma(200, shape=0.75, scale=1), ncol=10),
	ds  = matrix(runif(200), ncol=10),
	Ns  = matrix(rpois(200, lambda=30), ncol=10))
rownames(obj$mus) <- paste("gene", 1:nrow(obj$mus), sep="")
rownames(obj$rs) <- paste("gene", 1:nrow(obj$mus), sep="")
rownames(obj$ds) <- paste("gene", 1:nrow(obj$mus), sep="")
rownames(obj$Ns) <- paste("gene", 1:nrow(obj$mus), sep="")

obj <- apply_feature_selection(obj, c("gene1", "gene3", "gene5", "gene10"))
}
\keyword{feature selection}
