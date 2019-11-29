\name{tree_of_cells}
\alias{tree_of_cells}
\title{Build Tree of Cells}
\description{
Build Tree of Cells.
}
\usage{
tree_of_cells( merged_profiles, gene_features=rownames(merged_profiles$mus), dist_func=pearson_dist, plot_type="radial")
}
\arguments{
  \item{merged_profiles}{Output from: merge_profiles function.}
  \item{gene_features}{genes to use}
  \item{dist_func}{distance function for cell-type distances.}
  \item{plot_type}{how to plot the final tree.}
}
\details{
	uses the neighbour joining algorithm (from \code{phangorn} package) to build a tree of cell-types from a pair-wise distance matrix between cell-types. Cell-type distances are calculated using the provided function and the specificed genes.
}
\value{
	list containing the pair-wise distance matrix and the tree object.
}
\examples{
obj <- list(
	mus = matrix(rnorm(200), ncol=10),
	rs  = matrix(rgamma(200, shape=0.75, scale=1), ncol=10),
	ds  = matrix(runif(200), ncol=10),
	Ns  = matrix(rpois(200, lambda=30), ncol=10))
colnames(obj$mus) <- paste("celltype", 1:ncol(obj$mus), sep="")
rownames(obj$mus) <- paste("gene", 1:nrow(obj$mus), sep="")
colnames(obj$rs) <- paste("celltype", 1:ncol(obj$mus), sep="")
colnames(obj$ds) <- paste("celltype", 1:ncol(obj$mus), sep="")
colnames(obj$Ns) <- paste("celltype", 1:ncol(obj$mus), sep="")
out <- tree_of_cells(obj);
}
\keyword{distances, feature selection, tree}
