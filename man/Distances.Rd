\name{Distances}
\alias{base_dist}
\alias{pearson_dist}
\alias{spearman_dist}
\alias{phis_dist}
\alias{rhop_dist}
\alias{cosine_dist}
\title{Calculate Summary Statistics}
\description{
Calculates summary statistics to facilitate the fitting of a ZINB model, performs various data input checks.
}
\usage{
	base_dist(merged_profiles, type="euclidean")
	pearson_dist(merged_profiles, type="pearson")
	spearman_dist(merged_profiles, type="spearman")
	phis_dist(merged_profiles)
	rhop_dist(merged_profiles)
	cosine_dist(merged_profiles)
}
\arguments{
  \item{merged_profiles}{Output from: merge_profiles function.}
  \item{type}{type of distance to calcuate for methods with multiple options.}
}
\details{
	Calculates all pairwise distances between cell-types, using different distance measures.
	base_dist = wrapper of \code{stats::dist} methods
	pearson_dist/spearman_dist = 1 - correlation, wrapper of \code{stats::cor} methods.
	phis_dist = wrapper of \code{propr::phis}.
	rhop_dist = 1 - \code{propr::perb}.
	cosine_dist = 1 - \code{lsa::cosine}.
	
	For phis and rhop proportion methods, profiles are shifted up so that none of the values are negative.
}
\value{
	Distance matrix.
}
\examples{
obj <- list(
	mus = matrix(rnorm(200), ncol=10),
	rs  = matrix(rgamma(200, shape=0.75, scale=1), ncol=10),
	ds  = matrix(runif(200), ncol=10),
	Ns  = matrix(rpois(200, lambda=30), ncol=10))
d1 <- base_dist(obj, type="manhattan")
d2 <- spearman_dist(obj)
d3 <- pearson_dist(obj)
d4 <- phis_dist(obj)
d5 <- rhop_dist(obj)
d6 <- cosine_dist(obj)
}
\keyword{distance}
