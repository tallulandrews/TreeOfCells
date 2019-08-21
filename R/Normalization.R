quantile_norm <- function (merged_profiles, conditions=NULL) {
	require(preprocessCore)
	for (i in 1:length(merged_profiles)) {
		tmp <- normalize.quantiles(merged_profiles[[i]])
		colnames(tmp) <- colnames(merged_profiles[[i]])
		rownames(tmp) <- rownames(merged_profiles[[i]])
		merged_profiles[[i]] <- tmp
	}
	return(merged_profiles)
}

scnorm <- function (merged_profiles, conditions=rep(1, ncol(merged_profiles$mus))) {
	require(SCnorm)
	for (i in 1:length(merged_profiles)) {
		tmp <- SCnorm(merged_profiles[[i]], Conditions=conditions, FilterCellNum=10, NCores=1, FilterExpression=0.1)
		colnames(tmp) <- colnames(merged_profiles[[i]])
		rownames(tmp) <- rownames(merged_profiles[[i]])
		merged_profiles[[i]] <- tmp
	}
	return(merged_profiles)
}

sctrans <- function (merged_profiles, conditions=NULL) {
	require("sctransform")
	merged_profiles$mus <- suppressWarnings(sctransform::vst(merged_profiles$mus)$y)
	return(merged_profiles);
}

log_transform <- function(merged_profiles, conditions=NULL) {
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
	logprofiles$rs <- logrs;
	return(logprofiles);
}

total_norm <- function(merged_profiles, min_rs = 10^-10, conditions=NULL) {
	SF <- Matrix::colSums(merged_profiles$mus)
	SF <- SF/median(SF);
	v <- merged_profiles$mus + merged_profiles$mus^2/merged_profiles$rs
	normed_mus <- t(t(merged_profiles$mus)/SF)
	normed_v <- t(t(merged_profiles$mus)/(SF^2))
	normed_rs <- (normed_mus^2)/(normed_v - normed_mus)
	normed_rs[normed_rs < min_rs] <- min_rs
	normed_profiles <- merged_profiles
	normed_profiles$mus <- normed_mus
	normed_profiles$rs <- normed_rs
	return(normed_profiles);
}

gene_length_norm <- function(merged_profiles, gene_lengths, conditions=NULL){
	GF <- gene_lengths$length[match(rownames(merged_profiles$mus), gene_lengths$symbol)];
	GF[is.na(GF)] <- mean(GF[!is.na(GF)]);
	merged_profiles$mus <- merged_profiles$mus/GF*median(GF);
	return(merged_profiles);
}
