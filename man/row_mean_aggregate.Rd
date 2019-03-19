\name{row_mean_aggregate}
\alias{row_mean_aggregate}
\alias{row_var_aggregate}
\title{Row Means/Variances by Factor}
\description{
	Efficiently calculates the row means or row variances of a matrix split by a factor.
}
\usage{
	row_mean_aggregate(mat, groups)
	row_var_aggregate(mat, groups)
}
\arguments{
  \item{mat}{a matrix}
  \item{groups}{a factor of length ncol(mat) to split the matrix by.}
}
\details{
	Calculates the row means[variances] for columns belonging to each level of the provided groups. Uses an approach much more efficient than \code{aggregate}. Requires at least 2 samples from each group.
}
\value{
	a matrix of dimension nrow(mat) x length(levels(groups)) of meanas[variances] for each row for each group.
}
\examples{
	groups <- c("A", "A", "B", "B", "C","C","C")
	vals <- matrix(rnorm(70), ncol=7)
	group_means = row_mean_aggregate(vals, groups)
	group_vars = row_var_aggregate(vals, groups)
}
