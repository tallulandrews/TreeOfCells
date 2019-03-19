\name{tree.means}
\alias{tree.means}
\title{Generate the mean expression of genes through a tree}
\description{
	Generates mean expression values genes at each node through a tree.
}
\usage{
	tree.means(tree, de.sigma=0.5, prop_genes=0.2, root_means=2^rnorm(20000, mean=6, sd=2))
}
\arguments{
  \item{tree}{a \code{data.tree} object of the root of the tree.}
  \item{de.sigma}{standard deviation of log2 fold-change for differential expression between cell-types.}
  \item{prop_genes}{proportion of genes differentially expressed between each layer of the tree.}
  \item{root_means}{mean expression of genes at the root of the tree.}

}
\details{
	Takes a tree created by \code{create.rnd.tree} and uses \code{child.means} to calculate mean gene expression values at each node which are added as the "gene_means" attribute to each node.

	Custom created trees can be used as well as long as each node contains a "dist" attribute. 
}
\value{
	nothing, trees are passed as references so modifications will be visible to the existing tree variable.
}
\examples{
tree <- create.rnd.tree(dept=3);
tree.means(tree, root_means=2^rnorm(20, mean=6, sd=2))
tree$Get("gene_means")
}
\keyword{simulation, tree}
\seealso{data.tree, create.rnd.tree, tree.means, child.means}
