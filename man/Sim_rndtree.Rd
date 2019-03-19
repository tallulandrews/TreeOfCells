\name{create.rnd.tree}
\alias{create.rnd.tree}
\title{Create a random tree}
\description{
	Creates a random tree using specified parameters.
}

\usage{
	create.rnd.tree(depth, dist.gamma.params=list(shape=40, scale=0.05), split.prob=0.5)
}
\arguments{
  \item{depth}{Number of layers of cell-types in the tree.}
  \item{dist.gamma.params}{named list of parameters of a gamma distribution for distances between layers.}
  \item{split.prob}{Probability of a lineage splitting into two at each layer of the tree.}
}
\details{
	Uses the \code{data.tree} package to build a randomly generated tree-structure. Starts with a single root and builds down each layer adding child nodes each with a unique "name" attribute and a "dist" attribute specifying its distance from its parent. The \code{depth} arguement specifies the number of additional layers beyond the root. 
}
\value{
	A \code{data.tree} object of the root of the randomly created tree.
}
\examples{
sim_tree <- create.rnd.tree()
sim_tree <- create.rnd.tree(depth=5, split.prob=0.8)
}
\keyword{simulation, tree}
