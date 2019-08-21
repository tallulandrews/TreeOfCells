base_dist <- function(merged_profiles, type="euclidean") {
	as.matrix(dist(t(merged_profiles$mus), method=type));
}
cosine_dist <- function(merged_profiles) {
	require("lsa")
	lsa::cosine(merged_profiles$mus)

}
cor_dist <- function(merged_profiles, type="pearson") {
	1-as.matrix(cor(merged_profiles$mus, method=type));
}
KL_NegBinom <- function(merged_profiles) {
	log_factorial <- function(x) {
		# Srinivasa Ramanujan approximation: https://en.wikipedia.org/wiki/Stirling%27s_approximation#cite_note-Mortici2011-2-16
		# cite Cristinel Martici 2011?
		return( x*log(x)-x+1/6*log(8*x^3+4*x^2+x+1/30) + 1/2*log(pi) )
	}
	# Even with this I find no closed-form solution!
	

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


approx_EMD <- function(params1, params2, ntop=1000) {
	# still way too slow :(
	genes1 <- names(sort(params1$mus, decreasing=TRUE)[1:ntop])
	genes2 <- names(sort(params2$mus, decreasing=TRUE)[1:ntop])
	genes <- union(as.character(genes1), as.character(genes2));
	
	gene_scores <- sapply(which(names(params1$mus) %in% genes), function(i){Wasserstein_EMD(list(mu=params1$mus[i], r=params1$rs[i]), list(mu=params2$mus[i], r=params2$rs[i]))})
	return(mean(gene_scores));
}

get_bins <- function(params1, params2, q=0.01, nbins=100) {
	low<-min(qnbinom(q, mu=params1$mu, size=params1$r), qnbinom(q, mu=params2$mu, size=params2$r))
	high<-max(qnbinom(1-q, mu=params1$mu, size=params1$r), qnbinom(1-q, mu=params2$mu, size=params2$r))
	bins <- seq(from=low, to=high, length=nbins);
	bins <- unique(round(bins));
	return(bins);
}

Wasserstein_EMD <- function(params1, params2, err  = 10^-3) {
	# From: https://en.wikipedia.org/wiki/Earth_mover%27s_distance
	# no work different ranges used for different comparisons
	
	emd <- 0
	diff <- 0
	c <- 0
	#while(c < params1$mu | c < params2$mu | (emd[length(emd)] > err & diff > err)) {
	#	diff <- pnbinom(c, mu=params1$mu, size=params1$r) - pnbinom(c, mu=params2$mu, size=params2$r)	
	#	if (diff < err){c<-c+1; next;}
	#	emd_next = diff + emd[length(emd)];
	#	emd <- c(emd, emd_next);
	#	c <- c+1;
	#}

	bins <- get_bins(params1, params2);
	for (b_i in 2:length(bins)) {
		p1 <- qnbinom(bins[b_i], mu=params1$mu, size=params1$r) - qnbinom(bins[b_i-1], mu=params1$mu, size=params1$r)
		p2 <- qnbinom(bins[b_i], mu=params2$mu, size=params2$r) - qnbinom(bins[b_i-1], mu=params2$mu, size=params2$r)
		diff <- p1-p2
		emd_next = diff + emd[length(emd)];
		emd <- c(emd, emd_next);
	}
	return(sum(abs(emd))/length(emd))
}

Hellinger_distance <- function(params1, param2, err= 10^-3) {
	# https://en.wikipedia.org/wiki/Hellinger_distance
	bins <- get_bins(params1, params2);
	for (b_i in 2:length(bins)) {
		p1 <- qnbinom(bins[b_i], mu=params1$mu, size=params1$r) - qnbinom(bins[b_i-1], mu=params1$mu, size=params1$r)
		p2 <- qnbinom(bins[b_i], mu=params2$mu, size=params2$r) - qnbinom(bins[b_i-1], mu=params2$mu, size=params2$r)
		tot <- tot+(sqrt(p1)-sqrt(p2))^2
	}
	return(1/sqrt(2)*sqrt(tot))
	# for poisson: sqrt(1-e^(-1/2*(sqrt(lambda1)-sqrt(lambda2))^2))
}
Hellinger_poisson <- function(params1, params2) {
	# Acceptable speed.
	 gene_scores <- sapply(1:length(params1$mus), 
		function(i){
			lambda1 <- params1$mus[i];
			lambda2 <- params2$mus[i];
			return(sqrt(1-exp(-1/2*(sqrt(lambda1) - sqrt(lambda2))^2)))
		}
	)
	return(mean(gene_scores));
}

phi_s <- function(merged_profiles){
	require(propr)
	lowest <- min(merged_profiles$mus)
	if (lowest < 0) {
		merged_profiles$mus <- merged_profiles$mus+abs(lowest)
	}
	propr::phis(merged_profiles$mus, select = colnames(merged_profiles$mus))@matrix
}
rho_p <- function(merged_profiles) {
	require(propr)
	lowest <- min(merged_profiles$mus)
	if (lowest < 0) {
		merged_profiles$mus <- merged_profiles$mus+abs(lowest)
	}
	1-propr::perb(merged_profiles$mus, select=colnames(merged_profiles$mus))@matrix
}
