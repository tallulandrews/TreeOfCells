run_TM_benchmark <- function(test_fun, norm_facs=total_norm, norm_10x=norm_facs, norm_final=norm_facs) {
	dat10x <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/ForTreeOfCells/All_profiles_10X.rds")
	datFACs <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/ForTreeOfCells/All_profiles_FACS.rds")

	# Norm
	dat10x <- norm_10x(dat10x);
	datFACs <- norm_facs(datFACs);

	# Match columns
	common <- colnames(dat10x$mus)[colnames(dat10x$mus) %in% colnames(datFACs$mus)]
	n_common_col <- length(common)
	reorder <- match(common, colnames(dat10x$mus));
	for (i in 1:length(dat10x)) {
		dat10x[[i]] <- dat10x[[i]][,reorder]
		colnames(dat10x[[i]]) <- paste("x10", colnames(dat10x[[i]]), sep="_")
	}
	reorder <- match(common, colnames(datFACs$mus));
	for (i in 1:length(datFACs)) {
		datFACs[[i]] <- datFACs[[i]][,reorder]
		colnames(datFACs[[i]]) <- paste("FACS", colnames(datFACs[[i]]), sep="_")
	}
	# Match rows
	common <- rownames(dat10x$mus)[rownames(dat10x$mus) %in% rownames(datFACs$mus)]
	n_common_row <- length(common)
	reorder <- match(common, rownames(dat10x$mus));
	for (i in 1:length(dat10x)) {
		dat10x[[i]] <- dat10x[[i]][reorder,]
	}
	reorder <- match(common, rownames(datFACs$mus));
	for (i in 1:length(datFACs)) {
		datFACs[[i]] <- datFACs[[i]][reorder,]
	}
	# Combine
	dat_comb <- list();
	for (i in 1:length(datFACs)) {
		dat_comb[[i]] <- cbind(datFACs[[i]], dat10x[[i]]);
	}
	names(dat_comb) <- names(datFACs)
	dat_comb <- norm_final(dat_comb, conditions=c(rep(1, ncol(datFACs[[1]])), rep(2, ncol(dat10x[[1]]))));

	
	# Calculate distances
	dist_all <- test_fun(dat_comb)
	dist_cross <- dist_all[1:n_common_col, n_common_col+(1:n_common_col)]
	dist_facs <- dist_all[1:n_common_col, 1:n_common_col]
	dist_10x <- dist_all[n_common_col+(1:n_common_col), n_common_col+(1:n_common_col)]

	# Calculate statistics
	a <- as.vector(dist_10x)
	b <- as.vector(dist_facs)

	diag(dist_facs) <- Inf
	diag(dist_10x) <- Inf

	min_facs <-apply(dist_facs,1,min)
	min_10x <-apply(dist_10x,1,min)
	n_best <- sum(diag(dist_cross) < min_facs & diag(dist_cross) < min_10x)
	n_facs <- sum(diag(dist_cross) < min_facs)
	n_10x <- sum(diag(dist_cross) < min_10x)
	min_cross <- apply(dist_cross, 1, min);
	n_match <- sum(min_cross == diag(dist_cross))
	
	reg <- lm(a~b+0)
	return(list(
#		dist_all=dist_all,
#		dist_10x = dist_10x, 
#		dist_facs=dist_facs, 
		sim=cor(a, b), 
		diff=sum(abs(rank(a) - rank(b))),
		var=sum(((a-b)^2))/length(a)^2,
		rsq=summary(reg)$r.squared,
		prop_recip=n_best/length(common), 
		prop_facs=n_facs/length(common), 
		prop_10x=n_10x/length(common),
		prop_match=n_match/length(common)
		));
}

run_TM_benchmark2 <- function(test_fun) {
	dat10x <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/ForTreeOfCells/All_profiles_10X.rds")
	datFACs <- readRDS("/lustre/scratch117/cellgen/team218/TA/scRNASeqDatasets/TabulaMuris/ForTreeOfCells/All_profiles_FACS.rds")
	common <- colnames(dat10x$mus)[colnames(dat10x$mus) %in% colnames(datFACs$mus)]
	reorder <- match(common, colnames(dat10x$mus));
	for (i in 1:length(dat10x)) {
		dat10x[[i]] <- dat10x[[i]][,reorder]
	}
	reorder <- match(common, colnames(datFACs$mus));
	for (i in 1:length(datFACs)) {
		datFACs[[i]] <- datFACs[[i]][,reorder]
	}
	dat_comb <- list();
	for (i in 1:length(datFACs)) {
		dat_comb[[i]] <- cbind(datFACs[[i]], dat10x[[i]]);
	}

	
	
	dist_all <- test_fun(dat_comb)
	dist_cross <- dist_all[1:length(common), length(common)+(1:length(common))]
	dist_facs <- dist_all[1:length(common), 1:length(common)]
	dist_10x <- dist_all[length(common)+(1:length(common)), length(common)+(1:length(common))]
	diag(dist_facs) <- Inf
	diag(dist_10x) <- Inf

	min_facs <-apply(dist_facs,1,min)
	min_10x <-apply(dist_10x,1,min)
	n_best <- sum(dist_all < min_facs & dis_all < min_10x)
	n_facs <- sum(dist_all < min_facs)
	n_10x <- sum(dist_all < min_10x)
	return(list(
		#dist_all=dist_all,
		prop_recip=n_best/length(common), 
		prop_facs=n_facs/length(common), 
		prop_10x=n_10x/length(common)));
}

