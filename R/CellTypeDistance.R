# Utility Functions
factor_counts <- function(vec) {
        vec <- as.factor(vec)
        x <- split(seq(length(vec)), vec)
        result <- sapply(x, function(a) length(a))
        return(result);
}

row_mean_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[,a]))

        return(result);
}

row_var_aggregate <- function(mat, groups) {
        # Much faster version
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) matrixStats::rowVars(MAT[,a]))

        return(result);
}

merge_lists <- function(list1, list2) {
	if ( sum( names(list1) %in% names(list2) ) > 0) {
		warning("overlaping names in lists")
	}
	return(c(list1, list2))
}

cosine_dist <- function(x,y) {x %*% y / sqrt(x%*%x * y%*%y)}

celltype_stats <- function(expr_mat, celltypes) {
	# I think rowMeans and rowSums work for sparse matrices too.
	# need to catch cases with only one cell in a cell-type
	celltypes <- as.factor(celltypes);
	x <- split(seq(ncol(expr_mat)), celltypes)
	
	ns <- sapply(x, length)
	if (min(ns) < 2) {
		remove_types <- names(ns)[ns < 2]
		remove_cells <- celltypes %in% remove_types
		expr_mat <- expr_mat[,!remove_cells]
		celltypes <- celltypes[!remove_cells]
		celltypes <- factor(celltypes)
		warning(paste("Warning: removing", sum(remove_cells), "cells from", length(remove_types), "celltypes because less than two cells in the celltype."))
	}
	means <- sapply(x, function(a) rowMeans(expr_mat[,a]))
	detect <- sapply(x, function(a) rowSums(expr_mat[,a] > 0))
	expr_mat <- expr_mat-means[,as.numeric(celltypes)];
	vars <- sapply(x, function(a) rowSums(expr_mat[,a]^2))/(ns-1);
	return(list(m=means, v=vars, d=detect, n=ns));
}

# Distance Functions
celltype_diff_pairs <- function(expr_mat, celltypes, mt.method="none", q.value.threshold=0.05, filter_pct=0.5) { # Tested
	celltypes <- as.factor(celltypes)
	type_stats <- celltype_stats(expr_mat, celltypes);

	out <- list()

	require("proxy")
	for (t in levels(celltypes)){
		mat <- cbind(type_stats$m[,t], type_stats$v[,t],
			     rep(type_stats$n[t], times=nrow(type_stats$m)), 
			     type_stats$d[,t]/type_stats$n[t])
		rownames(mat) <- rownames(expr_mat);
		mat <- mat[mat[,2] > 0,]
		filter <- quantile(mat[,4], probs=c(filter_pct/2,1-filter_pct/2))
		mat <- mat[mat[,4] > filter[2] | mat[,4] < filter[1],]

		my_dist_fun <- function(x,y) { pnorm( abs(x[1]-y[1])/(sqrt(x[2]+y[2])/sqrt(x[3])), lower.tail=FALSE ) }
		D <- proxy::dist(mat, my_dist_fun)
		D <- as.matrix(D);
		diag(D) <- 1;
		if (mt.method == "none") { # Is this going to work?
			Q <- D
		} else {
			Q <- matrix(p.adjust(D, method=mt.method), ncol=ncol(D))
		}
		pairs <- which(Q < q.value.threshold, arr.ind=TRUE);
		q.value <- Q[pairs[,1:2]]
		p.value <- D[pairs[,1:2]]

		# Order genes column1 = low, column2 = high
		flip <- mat[pairs[,1], 1] > mat[pairs[,2],1] # column1 = low, column2=high
		tmp <- cbind(pairs[flip,2], pairs[flip,1])
		pairs[flip,] <- tmp;
		
		res <- cbind(rownames(D)[pairs[,1]], colnames(D)[pairs[,2]], p.value, q.value );
		colnames(res) <- c("lowGene", "hiGene", "p.value", "q.value")
		out[[t]] <- res;
	}
	return(out);
}

congruent_pair_distance <- function(list_of_gene_pair_matrices) {
	# output from celltype_diff_pairs
	my_pair_dist<-function(x,y){
		agree <- sum( paste(x$lowGene, x$hiGene) %in% paste(y$lowGene, y$hiGene) )
		disagree <- sum( paste(x$lowGene, x$hiGene) %in% paste(y$hiGene, y$lowGene) )
		if (agree == 0) {agree <- 10^-10}
		return(disagree/agree);
	}
	D <- proxy::dist(list_of_gene_pair_matrices, method=my_pair_dist)
	return(D)
}

summarize_pairs <- function(gene_pair_matrix, suppress.plot=TRUE) {
	all <- nrow(gene_pair_matrix);
	uni <- unique( c(gene_pair_matrix[,1], gene_pair_matrix[,2]))
	hist_data <- sapply(uni, function(g){sum(gene_pair_matrix[,1]==g)+sum(gene_pair_matrix[,2]==g)})
	if (!suppress.plot) {
		hist(hist_data, xlab="Number of pairs", ylab="Number of genes", col="grey80");	
	}
	return(list(all=all, unique=uni, dist=hist_data)) # check that hist_data keeps gene names
}

correct_gene_length <- function(gene_lengths, expression_vec, is.log=FALSE, suppress.plot=TRUE){
	if (!is.log) {
		expression_vec <- log2(expression_vec+1)
	}
	expression <- expression_vec[expression_vec > 0]
	gene_lengths <- log2(gene_lengths)
	gene_lengths <- gene_lengths[names(gene_lengths) %in% names(expression)]
	reg <- lm(expression ~ gene_lengths)

	if (!suppress.plot) {
		plot(gene_lengths, expression, xlab="log gene length", ylab="log expression")
		abline(reg, col="red");
		abline(h=median(expression), col="red", lty=2)
	}

	corrected <- median(expression)+reg$residuals;
	corrected[corrected < 0] <- 0;
	expression_vec[expression_vec > 0] <- corrected;

	if (is.log) {
		return(expression_vec);
	} else {
		return(2^expression_vec - 1)
	}
}

get_gene_length <- function(species=c("Mmus", "Hsap")){
	require("biomaRt")

}

# Testing Functions

down_sample_counts <- function(expr_mat, depth=min(colSums(expr_mat))) {
    down_sample <- function(x) {
        prob <- depth/sum(x)
        return(
            unlist(
                lapply(
                    x,
                    function(y) {
                        rbinom(1, y, prob)
                    }
                )
            )
        )
    }
    down_sampled_mat <- apply(expr_mat, 2, down_sample)
    return(down_sampled_mat)
}

down_sample_cells <- function(expr_mat, ncells) {
	subset <- sample(1:ncol(expr_mat), ncells)
	return(expr_mat[,subset])
}

down_sample_genes <- function(expr_mat, ngenes) {
	subset <- sample(1:nrow(expr_mat), ngenes)
	return(expr_mat[subset,])
}
