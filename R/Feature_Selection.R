fold_change_fs <- function(merged_profiles, diff_t=2, max_t=1) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	fc <- log2(apply(merged_profiles$mus+1, 1, max)) - log2(apply(merged_profiles$mus+1, 1, min))
	max_x <- log2(apply(merged_profiles$mus+1, 1, max))
	fs <- rownames(merged_profiles$mus)[fc > diff_t & max_x > max_t]
	return(fs)
}

distribution_olap_fs <- function(merged_profiles, quantile=0.5){
	merged_profiles <- check_merged_profiles(merged_profiles)
	min_r <- 10^-5;
	max_r <- 10^10;
	lowest <- apply(merged_profiles$mus, 1, function(x) {which(x==min(x))[1]})
	highest <- apply(merged_profiles$mus, 1, function(x) {which(x==max(x))[1]})
	olap_quantile <- function(i) {
		# upper threshold of lowest value
		m_l <- merged_profiles$mus[i,lowest[i]];
		r_l <- merged_profiles$rs[i,lowest[i]];
		if (is.na(r_l) | r_l < min_r | r_l >= max_r) {
			bottom_up <- qpois(quantile, lambda=m_l);
		} else {	
			bottom_up <- qnbinom(quantile, mu=m_l, size=r_l)
		}
		if (is.na(bottom_up)) {bottom_up <- 0}
		# lower threshold of highest value
		m_u <- merged_profiles$mus[i,highest[i]];
		r_u <- merged_profiles$rs[i,highest[i]];
		if (is.na(r_u) | r_u < min_r | r_u >= max_r) {
			top_down <- qpois(1-quantile, lambda=m_u)
		} else {	
			top_down <- qnbinom(1-quantile, mu=m_u, size=r_u)
		}
		if (is.na(top_down)) {top_down <- -1}
		return(top_down-bottom_up)
	}
	test <- sapply(1:nrow(merged_profiles$mus), olap_quantile)
	return(rownames(merged_profiles$mus)[test > 0]);
}

dropout_fs <- function(merged_profiles, diff=0.75) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	p_zero_d <- merged_profiles$ds
	p_zero_nb <-(1+merged_profiles$mus/merged_profiles$rs)^(-merged_profiles$rs);
	p_zero <- p_zero_d+p_zero_nb
	Ns <- apply(merged_profiles$Ns, 2, max);
	stdev <- t(t(p_zero*(1-p_zero))/Ns)
	lowest <- apply(p_zero, 1, function(x) {which(x==min(x))[1]})
        highest <- apply(p_zero, 1, function(x) {which(x==max(x))[1]})

	prop_diff <- function(i) {
		p_l <- p_zero[i,lowest[i]]
		p_h <- p_zero[i, highest[i]]
		Z <- ((p_h-p_l)-diff)/sqrt(stdev[i,lowest[i]]+stdev[i,highest[i]])
		return(Z);
	}
	test <- sapply(1:nrow(merged_profiles$mus), prop_diff)
	test[is.na(test)] <- 1;
	pvals <- pnorm(test, lower.tail=F)
	return(rownames(merged_profiles$mus)[pvals <  0.05/(nrow(merged_profiles$mus) * choose(ncol(merged_profiles$mus),2))])
}


FS_fold_change_sig <- function(merged_profiles, diff_t=2, max_t=1, quantile=0.75) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	# What is this? Appoximating NB with Normal???
	top_down_lim <- merged_profiles$mus - sqrt(merged_profiles$mus+merged_profiles$mus^2/merged_profiles$rs)/sqrt(merged_profiles$Ns)*2
	top <- unlist( apply(top_down_lim, 1, function(x) {max(x, na.rm=TRUE)}) )
	top[top < 0] <- 0
	bottom_up_lim <- merged_profiles$mus + sqrt(merged_profiles$mus+merged_profiles$mus^2/merged_profiles$rs)/sqrt(merged_profiles$Ns)*2
	bottom <- unlist( apply(bottom_up_lim, 1, function(x) {min(x, na.rm=TRUE)}) )
	bottom[!is.finite(bottom)] <- 1

	sig <- top > bottom


	top_mus <- merged_profiles$mus[cbind(1:length(top), top)] 
	bottom_mus <- merged_profiles$mus[cbind(1:length(bottom), bottom)] 
	top_v <- top_mus + top_mus^2/merged_profiles$rs[cbind(1:length(top), top)]
	bottom_v <- bottom_mus + bottom_mus^2/merged_profiles$rs[cbind(1:length(bottom), bottom)]
	top_N <- merged_profiles$Ns[cbind(1:length(top), top)]
	
	fc <- log2(apply(merged_profiles$mus+1, 1, max)) - log2(apply(merged_profiles$mus+1, 1, min))
	max_x <- log2(apply(merged_profiles$mus+1, 1, max))
	return(fc > diff_t & max_x > max_t)
}

top_diff_fs <- function(merged_profiles, ntop=3000) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	mat <- get_val_mat_for_fs(merged_profiles);
	diff <- apply(mat, 1, max) - apply(mat,1,min);
	keep <- rownames(mat)[diff > quantile(diff, probs=1-ntop/length(diff))]
	return(keep);
}
	
apply_feature_selection <- function(merged_profiles, features) {
	merged_profiles <- check_merged_profiles(merged_profiles)
	for (i in 1:length(merged_profiles)) {
		merged_profiles[[i]] <- merged_profiles[[i]][rownames(merged_profiles[[i]]) %in% features,]
	}
	return(merged_profiles)
}

get_val_mat_for_fs <- function(merged_profiles) {
	merged_profiles <- check_merged_profiles(merged_profiles)
        if ("norm" %in% names(merged_profiles)) {
                return(merged_profiles$norm);
        }
        else if ("mus" %in% names(merged_profiles)) {
                return(merged_profiles$mus);
        } else {
                stop("No recognized matrix for distance calculation.")
        }
}

