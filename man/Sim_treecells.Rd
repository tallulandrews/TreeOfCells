\name{tree.n_cells}
\alias{tree.n_cells}
\title{Generate the mean expression of genes through a tree}
\description{
	Generates mean expression values genes at each node through a tree.
}
\usage{
	tree.n_cells(tree, min=10, max=1000, rnd.fun=function(n){runif(n, min=min, max=max)}, ...)
}
\arguments{
  \item{tree}{a \code{data.tree} object of the root of the tree.}
  \item{min}{minimum number of cells per population.}
  \item{max}{maximum number of cells per population.}
  \item{rnd.fun}{random number generator function for number of cells per population}
  \item{...}{additional arguments for rnd.fun}
}

\details{
	Takes a tree created by \code{create.rnd.tree} and randomly generates the number of cells at each node which are added as the "n_cells" attribute to each node.

	Custom created trees can be used as well as long as they are a \code{data.tree} object. 
}
\value{
	nothing, trees are passed as references so modifications will be visible to the existing tree variable.
}
\examples{
tree <- create.rnd.tree(dept=3);
tree.n_cells(tree) # default uniform
tree.n_cells(tree, min=5, rnd.fun=function(n){rnbinom(n, mean=50, size=0.1)}) # custom negative binomial;
tree.n_cells(tree, min=5, rnd.fun=function(n){rep(50, n)}) # 50 cells of each node.
tree$Get("n_cells")
}
\keyword{simulation, tree}
\seealso{data.tree, create.rnd.tree, tree.n_cells}
