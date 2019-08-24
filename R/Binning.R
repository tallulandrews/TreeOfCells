mclust_bins <- function(normed) {
	# Tested
	require(mclust);
	get_bins <- function(x) {
		if (var(x)==0) {return(rep(1, length(x)))}
		out <- Mclust(x, G=1:(ncol(normed$mus)/10), modelName="V")
		lims <- sapply(1:out$G, function(y){c(min(x[out$classification==y]), max(x[out$classification==y]) )})
		if (ncol(lims)==1){return(rep(1, length(x)))}
		check <- which(lims[1,2:ncol(lims)] < lims[2,1:(ncol(lims)-1)])
		while(length(check) > 0) {
			change <- out$classification %in% check;
			out$classification[change] <- out$classification[change]+1
			lims <- sapply(min(out$classification):max(out$classification), 
					function(y){c(min(x[out$classification==y]), 
						max(x[out$classification==y]) )})
			if (ncol(lims)==1){return(rep(1, length(x)))}
			check <- which(lims[1,2:ncol(lims)] < lims[2,1:(ncol(lims)-1)])+min(out$classification)-1
		}
		return(out$classification-min(out$classification)+1)
	}

	binned <- t(apply(normed$mus, 1, get_bins))
}
