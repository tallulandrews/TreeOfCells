\name{sim.tree}
\alias{sim.tree}
\title{Simulate data from a tree.}
\description{
	Simulates a count matrix from a given tree, and stores it along with appropriate metadata in a SingleCellExperiment object.
}
\usage{
	sim.tree(tree, mean.size.fun=function(a){return(rep(0.2,length(a)))}, lib.size.gamma.params=list(shape=2, scale=2), K=NULL)
}
\arguments{
  \item{tree}{a \code{data.tree} object of the root of the tree.}
  \item{mean.size.fun}{a function to generate the size parameter (1/dispersion) of a gene's negative binomial distribution given its mean.}
  \item{lib.size.gamma.params}{named list of "shape" and "scale" parameters for a gamma distribution to draw cell-specific library sizes from.}
  \item{K}{Michaelis-Menten parameter for the zero-inflation added by \code{add_dropouts}.}
}
\details{
	Takes a \code{data.tree} tree and uses \code{ZINB.sim} to generate counts for the population described by each node of the tree. These are combined into a single count matrix. This is stored in a SingleCellExperiment object along with a cell-specific metadata table that contains the randomly generated library sizes, the name of the node the cell belongs to and whether it is a leaf population (colData). In addition, the true mean expression of each gene in each cell population is stored in a second metadata table (rowData). 

	The true tree relationships between cell-populations is stored in metadata slot of the SingleCellExperiment object as "true_tree" a \code{data.tree} object.

	The input tree must contain "name", "gene_means", and "n_cells" attributes at each node. These are added by \code{create.rnd.tree}, \code{tree.means} and \code{tree.n_cells} respectively. Or they can be added by hand using the functions provided in the \code{data.tree} package.
}
\value{
	A SingleCellExperiment object containing:
		assays: a count matrix
		colData: containing "lib.size", "Group", "is.Leaf"
		rowData: containing true mean expression in each group.
		metadata$true_tree: a \code{data.tree} object of the input tree.
}
\examples{
tree <- create.rnd.tree(depth=2);
tree.means(tree);
tree.n_cells(tree);
sce <- sim.tree(tree);
}
\keyword{simulation, tree, zero-inflation, negative binomial}
\seealso{data.tree, create.rnd.tree, tree.means, tree.n_cells, ZINB.sim}
