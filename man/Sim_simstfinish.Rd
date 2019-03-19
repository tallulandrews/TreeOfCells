\name{sim.start.finish.defaults}
\alias{sim.start.finish.defaults}
\title{Simulate data from a tree.}
\description{
	Simulates data from a random tree using default values.
}
\usage{
	sim.start.finish.defaults()
}
\arguments{
}
\details{
	Simulates a random tree, add random mean expression values, and randomly generated numbers of cells for each popluation. Then generates a count matrix and SingleCellExperiment object of the result, using default values.
}
\value{
	A SingleCellExperiment object containing:
		assays: a count matrix
		colData: containing "lib.size", "Group", "is.Leaf"
		rowData: containing true mean expression in each group.
		metadata$true_tree: a \code{data.tree} object of the input tree.
}
\examples{
sce <- sim.start.finish.defaults();
}
\keyword{simulation, tree, zero-inflation, negative binomial}
\seealso{data.tree, create.rnd.tree, tree.means, tree.n_cells, ZINB.sim, sim.tree}
