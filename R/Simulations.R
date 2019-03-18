# Tree-like relationships between cell-types

# Zero inflation according to M3Drop model (bg) #
add_dropouts <- function(x,mu,K){
	# p_drop = 1-mu/(mu+K) = K/(mu+K)
	# p_no_drop = mu/(mu+K) = sum(x>0)/length(x)
	# mu*length(x)/sum(x>0) = mu+K
	expect_pos <- mu*length(x)/sum(x>0);
	p_drop <- K/expect_pos
	toss <- runif(n=length(x));
	x[toss < p_drop] <- 0;
	return(x);
}

# Fit gamma parameters to observed values #
fit_gamma <- function(x) {
	s = var(x)/mean(x)
	a = mean(x)/s
	return(list(shape=a, scale=s))
}


# Simulate from a ZINB (bg) #
ZINB.sim <- function(gene_means, lib.size, mean.size.fun, K=NULL) {
	lib.size <- lib.size/mean(lib.size);
	mu.mat <- gene_means %*% t(lib.size);
	size.vec <- mean.size.fun(gene_means);
	sim <- sapply(1:nrow(mu.mat), function(i) {sapply(mu.mat[i,], function(a){rnbinom(1, mu=a, size=size.vec[i])})})
	if (!is.null(K)) {
		sim <- sapply(1:nrow(sim), function(i){add_dropouts(sim[i,], gene_means[i], K)})
	}
	return(sim);
}

# Change in expression along branches (bg) #
child.means <- function(parent.means, relative.dist, sigma=0.5, prop_genes=0.2) {
	child <- log2(parent.means) + relative.dist*rnorm(length(parent.means), mean=0, sd=sigma)
	de.genes <- runif(length(child)) < prop_genes
	child[!de.genes] <- parent.means[!de.genes]
	return(2^child);	
}

# Simulate a random tree #
create.rnd.tree <- function(depth=4, dist.gamma.params=list(shape=40, scale=0.05), split.prob=0.5) {
	require(data.tree);
	type_tree <- Node$new("Root");
	layer <- c(type_tree);
	for (i in 1:depth) {
		new_layer = vector();
		for (n in 1:length(layer)) {
			if (runif(1) < split.prob) {
				# two children
				n_child <- 2
			} else {
				# one child
				n_child <- 1
			}
			for (c in 1:n_child) {
				child.name <- (paste(layer[[n]]$name, c, sep="_"))
				child.dist <- rgamma(1, shape=dist.gamma.params$shape, scale=dist.gamma.params$scale)
				child <- layer[[n]]$AddChild(child.name, de_dist=child.dist);
				new_layer <- c(new_layer, child)
			}
		}
		layer <- new_layer
	}
	return(type_tree)
}

# Set the mean expression of genes for each node of the tree #
tree.means <- function(tree, de.sigma=0.5, prop_genes=0.2, root_means=2^rnorm(20000, mean=6, sd=2)) {
	require(data.tree);
	tree$Do(function(n){
		if (n$isRoot) {
			n$gene_means <- root_means
		} else {
			n$gene_means <- child.means(n$parent$gene_means, n$de_dist, sigma=de.sigma, prop_genes=prop_genes)
		}
	})
}

# Randomly set the number of cells for each node of the tree #
tree.n_cells <- function(tree, min=10, max=1000, rnd.fun=function(n){runif(n, min=min, max=max)}, ...) {
	require(data.tree);
	tree$Do(function(n){
		cells <- round(rnd.fun(1, ...));
		if (cells > max) {cells <- max};
		if (cells < min) {cells <- min};
		n$n_cells <- cells;
	})
}

# Simulate the dataset from a given tree #
sim.tree <- function(tree, mean.size.fun=function(a){return(rep(0.2,length(a)))}, lib.size.gamma.params=list(shape=2, scale=2), K=NULL) {
	require(data.tree);
	all_nodes <- tree$Get("name")
	leaf_names <- tree$Get("name", filterFun=isLeaf)
	OUT <- vector()
	ANN <- vector()
	for (node_name in all_nodes) {
		node <- FindNode(tree, node_name);	
		lib.size <- rgamma(node$n_cells, shape=lib.size.gamma.params$shape, scale=lib.size.gamma.params$scale)
		mat <- t( ZINB.sim(node$gene_means, lib.size, mean.size.fun, K) )
		colnames(mat) <- paste(node$name, "_cell", 1:ncol(mat),sep="")
		ann <- cbind(lib.size, rep(node_name, node$n_cells), rep(node$isLeaf, node$n_cells))
		rownames(ann) <- colnames(mat);
		ANN <- rbind(ANN, ann)
		OUT<- cbind(OUT,mat)
	}
	require("SingleCellExperiment")
	colnames(ANN) <- c("lib.size", "Group", "is.Leaf")
	rownames(ANN) <- colnames(OUT);
	true_means <- tree$Get("gene_means")
	gene_names <- paste("Gene", 1:nrow(OUT), sep="")
	rownames(OUT) <- gene_names;
	rownames(true_means) <- gene_names;
	ANN <- data.frame(lib.size=as.numeric(ANN[,1]), Group=factor(ANN[,2]), is.Leaf=as.logical(ANN[,3]))
	SCE <- SingleCellExperiment(assays=list(counts=OUT), colData=ANN, rowData=true_means);
	SCE@metadata$group_tree <- Clone(tree);
	return(SCE);
}

# Do full simulation using default values #
sim.start.finish.defaults <- function() {
	tree <- create.rnd.tree()
	tree.means(tree);
	tree.n_cells(tree);
	sce <- sim.tree(tree)
	return(sce);
}
