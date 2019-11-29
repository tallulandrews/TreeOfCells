\name{Cell Type Stats}
\alias{row_mean_aggregate}
\alias{row_var_aggregate}
\alias{col_mean_aggregate}
\alias{col_var_aggregate}
\title{Row/Column Means/Variances by Factor}
\description{
	Efficiently calculates the row means or row variances or column means or columns variances of a matrix split by a factor.
}
\usage{
	row_mean_aggregate(mat, groups)
	row_var_aggregate(mat, groups)
	col_mean_aggregate(mat, groups)
	col_var_aggregate(mat, groups)
}
\arguments{
  \item{mat}{a matrix}
  \item{groups}{a factor of length ncol(mat) to split the matrix by.}
}
\details{
	Calculates the row[column] means[variances] for columns belonging to each level of the provided groups. Uses an approach much more efficient than \code{aggregate}. Requires at least 2 samples from each group.
}
\value{
	a matrix of dimension nrow(mat) x length(levels(groups)) of means[variances] for each row for each group.
	a matrix of dimension ncol(mat) x length(levels(groups)) of means[variances] for each column for each group.
}
\examples{
	groups <- c("A", "A", "B", "B", "C","C","C")
	vals <- matrix(rnorm(70), ncol=7)
	group_means = row_mean_aggregate(vals, groups)
	group_vars = row_var_aggregate(vals, groups)
	vals <- t(vals)
	group_means = col_mean_aggregate(vals, groups)
	group_vars = col_var_aggregate(vals, groups)
}
