FS_fold_change <- function(merged_profiles, diff_t=2, max_t=1) {
	fc <- log2(apply(merged_profiles$mus+1, 1, max)) - log2(apply(merged_profiles$mus+1, 1, min))
	max_x <- log2(apply(merged_profiles$mus+1, 1, max))
	return(fc > diff_t & max_x > max_t)
}

FS_fold_change_sig <- function(merged_profiles, diff_t=2, max_t=1) {
	low_lim <- merged_profiles$mus - sqrt(merged_profiles$mus+merged_profiles$mus^2/merged_profiles$rs)/sqrt(merged_profiles$Ns)*2
	top <- unlist( apply(low_lim, 1, function(x) {max(x, na.rm=TRUE)}) )
	top[top < 0] <- 0
	bottom_up_lim <- merged_profiles$mus + sqrt(merged_profiles$mus+merged_profiles$mus^2/merged_profiles$rs)/sqrt(merged_profiles$Ns)*2
	bottom <- unlist( apply(bottom_up_lim, 1, function(x) {min(x, na.rm=TRUE)}) )
	bottom[!is.finite(bottom)] <- 0

	top_mus <- merged_profiles$mus[cbind(1:length(top), top)] 
	bottom_mus <- merged_profiles$mus[cbind(1:length(bottom), bottom)] 
	top_v <- top_mus + top_mus^2/merged_profiles$rs[cbind(1:length(top), top)]
	bottom_v <- bottom_mus + bottom_mus^2/merged_profiles$rs[cbind(1:length(bottom), bottom)]
	top_N <- merged_profiles$Ns[cbind(1:length(top), top)]
	
	fc <- log2(apply(merged_profiles$mus+1, 1, max)) - log2(apply(merged_profiles$mus+1, 1, min))
	max_x <- log2(apply(merged_profiles$mus+1, 1, max))
	return(fc > diff_t & max_x > max_t)
}

olap_fs <- function(merged_profiles) {
	lowest_params <-
	highest_params  <-
	for (g in genes) {
		bins <- get_bins(lowest_params, highest_params)
		for (i in bins) {
			ps <- sapply(1:ncol(merged_profiles$mus), function(x){pnbinom(i, merged_profiels$mus[g,x], merged_profiles$rs[g,x])})
			ps/sum(ps)
####IMCOMPLETE IMCOMPLETE#####
		}
	}
}
