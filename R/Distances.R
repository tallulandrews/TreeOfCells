get_val_mat_for_dist <- function(merged_profiles) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	if ("scaled" %in% names(merged_profiles)) {
		return(merged_profiles$scaled);
	}
	else if ("norm" %in% names(merged_profiles)) {
		return(merged_profiles$norm);
	}
	else if ("mus" %in% names(merged_profiles)) {
		return(merged_profiles$mus);
	} else {
		stop("No recognized matrix for distance calculation.")
	}
}

base_dist <- function(merged_profiles, type="euclidean") {
	as.matrix(dist(t(get_val_mat_for_dist(merged_profiles)), method=type));
}
cosine_dist <- function(merged_profiles) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	1-lsa::cosine(get_val_mat_for_dist(merged_profiles))

}
pearson_dist <- function(merged_profiles, type="pearson") {
	merged_profiles <- check_merged_profiles(merged_profiles)
	1-as.matrix(cor(get_val_mat_for_dist(merged_profiles), method=type));
}
spearman_dist <- function(merged_profiles, type="spearman") {
	merged_profiles <- check_merged_profiles(merged_profiles)
	1-as.matrix(cor(get_val_mat_for_dist(merged_profiles), method=type));
}

KL_Poisson <- function(params1, params2) {
	# Acceptable speed.
	 gene_scores <- sapply(1:length(params1$mus), 
		function(i){
			lambda1 <- params1$mus[i];
			lambda2 <- params2$mus[i];
			KL_1_2 <- lambda1*log(lambda1/lambda2) - lambda1 + lambda2
			KL_2_1 <- lambda2*log(lambda2/lambda1) - lambda2 + lambda1
			return(KL_1_2+KL_2_1)
		}
	)
	return(mean(gene_scores));
}

dmat_wrapper <- function(merged_profiles, d_func) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	out <- matrix(0, ncol=ncol(merged_profiles$mus), nrow=nrow(merged_profiles$mus));
	for (i in 1:ncol(out)) {
		for (j in i:ncol(out)) {
			d <- d_func(list(mus=merged_profiles$mus[,i], rs=merged_profiles$rs[,i]),
				    list(mus=merged_profiles$mus[,j], rs=merged_profiles$rs[,j]))
			out[i,j] <-d
			out[j,i] <-d
		}
	}
	return(out)
}

get_bins <- function(params1, params2, q=0.01, nbins=100) {
	low<-min(qnbinom(q, mu=params1$mu, size=params1$r), qnbinom(q, mu=params2$mu, size=params2$r))
	high<-max(qnbinom(1-q, mu=params1$mu, size=params1$r), qnbinom(1-q, mu=params2$mu, size=params2$r))
	bins <- seq(from=low, to=high, length=nbins);
	bins <- unique(round(bins));
	return(bins);
}

phis_dist <- function(merged_profiles){
	merged_profiles <- check_merged_profiles(merged_profiles)
	mat <- get_val_mat_for_dist(merged_profiles)
	if (length(colnames(mat)) != ncol(mat)) {
		colnames(mat) <- paste("cell", 1:ncol(mat), sep="")
	}
	lowest <- min(mat)
	if (lowest < 0) {
		mat <- mat+abs(lowest)
	}
	propr::phis(mat, select=colnames(mat))@matrix
}
rhop_dist <- function(merged_profiles) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	mat <- get_val_mat_for_dist(merged_profiles)
	if (length(colnames(mat)) != ncol(mat)) {
		colnames(mat) <- paste("cell", 1:ncol(mat), sep="")
	}
	lowest <- min(mat)
	if (lowest < 0) {
		mat <- mat+abs(lowest)
	}
	1-propr::perb(mat, select=colnames(mat))@matrix
}
