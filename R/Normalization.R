quantile_norm <- function (merged_profiles, conditions=NULL) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	for (i in 1:length(merged_profiles)) {
		tmp <- preprocessCore::normalize.quantiles(merged_profiles[[i]])
		colnames(tmp) <- colnames(merged_profiles[[i]])
		rownames(tmp) <- rownames(merged_profiles[[i]])
		merged_profiles[[i]] <- tmp
	}
	merged_profiles$norm <- merged_profiles$mus;
	return(merged_profiles)
}

downsample_norm <- function (merged_profiles) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	#require("mclust")
	zeros <- Matrix::colSums(merged_profiles$mus == 0);
	zero_lim <- quantile(zeros, 0.90);
	zero_perc <- zero_lim/nrow(merged_profiles$mus)
	downsamp <- apply(merged_profiles$mus, 2, function(x){
		if (sum(x==0) < zero_lim) {
			thresh <- quantile(x, zero_perc);
			x[x <= thresh] <- 0;
		}
		return(x);
	})
	sf <- Matrix::colSums(downsamp);
	downsamp <- t(t(downsamp)/sf*median(sf))
	merged_profiles$norm <- downsamp;
	return(merged_profiles)
}
	

scnorm <- function (merged_profiles, conditions=rep(1, ncol(merged_profiles$mus))) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	for (i in 1:length(merged_profiles)) {
		tmp <- SCnorm::SCnorm(merged_profiles[[i]], Conditions=conditions, FilterCellNum=10, NCores=1, FilterExpression=0.1)
		colnames(tmp) <- colnames(merged_profiles[[i]])
		rownames(tmp) <- rownames(merged_profiles[[i]])
		merged_profiles[[i]] <- tmp
	}
	merged_profiles$norm <- merged_profiles$mus;
	return(merged_profiles)
}

scaled_regression <- function(merged_profiles, conditions=NULL) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	# general least squares regression 
	# with custom variances for each datapoint
	# based on the number of cells for that cell-type
	# (i.e. variance of sample mean)
	ns <- apply(merged_profiles$Ns,2,max)
	var_structure <- nlme::varFixed(~sqrt(ns));
	sf <- Matrix::colSums(merged_profiles$mus);
	my_reg <- function(x) {
		if (var(x)<=1^-10){return(x)}
		else {nlme::gls(x~sf, weights=var_structure)$residuals}
	}

	norm <- t(apply(merged_profiles$mus, 1, my_reg))
	rownames(norm) <- rownames(merged_profiles$mus)
	colnames(norm) <- colnames(merged_profiles$mus)
	merged_profiles$norm <- norm;
	norm <- t(scale(t(norm)))
	exclude <- is.na(Matrix::rowMeans(norm))
	merged_profiles$scaled <-norm
	for (i in 1:length(merged_profiles)) {
		merged_profiles[[i]] <- merged_profiles[[i]][!exclude,]
	}
	return(merged_profiles);
}

sctrans <- function (merged_profiles, conditions=NULL) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	merged_profiles$norm <- suppressWarnings(sctransform::vst(merged_profiles$mus)$y)
	return(merged_profiles);
}

log_transform <- function(merged_profiles, conditions=NULL) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	v <- merged_profiles$mus + merged_profiles$mus^2/merged_profiles$rs
	logmus <- log(merged_profiles$mus+1)
	logrs <- logmus^2/(log(v)/logmus - logmus)
	logrs[!is.finite(logrs)] <- 0;
	colnames(logmus) <- colnames(merged_profiles$mus)
	rownames(logmus) <- rownames(merged_profiles$mus)
	colnames(logrs) <- colnames(merged_profiles$rs)
	rownames(logrs) <- rownames(merged_profiles$rs)
	logprofiles <- merged_profiles;
	logprofiles$mus <- logmus;
	logprofiles$norm <- logmus;
	logprofiles$rs <- logrs;
	return(logprofiles);
}

total_norm <- function(merged_profiles, min_rs = 10^-10, conditions=NULL) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	SF <- Matrix::colSums(merged_profiles$mus)
	SF <- SF/median(SF);
	v <- merged_profiles$mus + merged_profiles$mus^2/merged_profiles$rs
	normed_mus <- t(t(merged_profiles$mus)/SF)
	normed_v <- t(t(merged_profiles$mus)/(SF^2))
	normed_rs <- (normed_mus^2)/(normed_v - normed_mus)
	normed_rs[normed_rs < min_rs] <- min_rs
	normed_profiles <- merged_profiles
	normed_profiles$mus <- normed_mus
	normed_profiles$norm <- normed_mus
	#normed_profiles$rs <- normed_rs
	return(normed_profiles);
}

gene_length_norm <- function(merged_profiles, gene_lengths, conditions=NULL){
	merged_profiles <- check_merged_profiles(merged_profiles)
	GF <- gene_lengths$length[match(rownames(merged_profiles$mus), gene_lengths$symbol)];
	GF[is.na(GF)] <- mean(GF[!is.na(GF)]);
	merged_profiles$norm <- merged_profiles$mus/GF*median(GF);
	return(merged_profiles);
}
