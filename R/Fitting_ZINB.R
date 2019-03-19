#### Fitting #####
calculate_summary_values <- function(counts) {
        if (sum(counts < 0) >0) {stop("Expression matrix contains negative values! Please provide raw UMI counts!")}
	if ( sum(counts >= 1) != sum(counts > 0) ) {stop("Error: Expression matrix is not integers! Please provide raw UMI counts.")}
#        if (sum(!is.integer(counts)) >0) {stop("Expression matrix is not integers! Please provide a matrix (not data.frame) raw UMI counts!")}

        tjs <- rowSums(counts, na.rm=T) # Total molecules/gene
	if (sum(tjs <= 0) > 0) {stop("Error: all genes must have at least one detected molecule.")}
        tis <- colSums(counts, na.rm=T) # Total molecules/cell
	if (sum(tis <= 0) > 0) {stop("Error: all cells must have at least one detected molecule.")}
        djs <- ncol(counts)-rowSums(counts > 0, na.rm=T) # Observed Dropouts per gene
        dis <- nrow(counts)-colSums(counts > 0, na.rm=T) # Observed Dropouts per cell
        nc <- length(counts[1,]) # Number of cells
        ng <- length(counts[,1]) # Number of genes
        total <- sum(tis, na.rm=T) # Total molecules sampled
        return(list(tis=tis, tjs=tjs, dis=dis, djs=djs, total=total,nc=nc,ng=ng));
}

convert_to_integer <- function(mat) {
        mat <- ceiling(as.matrix(mat))
        storage.mode(mat) <- "integer"
        mat <- mat[rowSums(mat) > 0,]
        return(mat)
}


fit_ZINB_to_matrix <- function(counts) {
	vals <- calculate_summary_values(counts);
	out <- lapply(1:vals$ng, function(g) {
		fit_zero_inflated_negative_binomial(counts[g,], g, vals)
		})
	return(out);
}

fit_zero_inflated_negative_binomial <- function(obs, g_row, vals, e=0.00001) {
	l_is <- vals$tis*vals$nc/vals$total # cell-specific weights due to difference in library size
	d0 <- vals$djs[g_row]/vals$nc # dropout-rate, initially all zeros due to dropouts
	max_r <- 10^10 # r = size parameter for R's nbinom

	d_curr <- d0; 
	d_prev <- -100;
	while( abs(d_curr-d_prev) > e ) {

		# Fit the NB #
		mu_j <- vals$tjs[g_row]/(vals$nc-d_curr*vals$nc) # dropouts only affect zeros to just change the denominator of the mean
		mu_ijs <- mu_j*l_is; #account for different library sizes
		weights <- rep(1, times=length(mu_ijs)) # rather than remove some columns 
							
		weights[obs == 0] <- (1-d_curr*vals$nc/vals$djs[g_row]) # lets down-weight the contribution of 
									# zeros for fitting the NB
		# Note: error = n*variance 
		# weight-adjusted observed error correcting for cell-specific library sizes
		obs_err <- sum( (obs - mu_ijs)^2*weights )

		# fit the dispersion as:
		# Observed error = sum of variances of cell-specific NB distributions
		# variance of NB = mean+mean^2/dispersion
		# rearrange to solve for dispersion (r)
		r_j <- sum( mu_ijs^2*weights )/(obs_err - sum(mu_ijs*weights))
		if (r_j <= 0) {r_j <- max_r}

		# probability of a zero in each cell given the NB
		p_ijs <- (1 + mu_ijs/r_j)^(-r_j)
		d_exp <- mean(p_ijs)

		# Update dropout rate (d) estimate
		d_prev <- d_curr
		d_curr <- (vals$djs[g_row] - d_exp*vals$nc)/vals$nc
		if (d_curr <= 0) {d_curr <- 0}
	}
	return(c(mu_j, r_j, d_prev));
}


fit_NB_to_matrix <- function(counts) {
	vals <- calculate_summary_values(counts)
	min_size <- 10^-10
	my_rowvar <- sapply(1:nrow(counts), function(i) {
		mu_is <- vals$tjs[i]*vals$tis/vals$total
		var(as.vector(unlist(counts[i,]))-mu_is)
		})
	size <- vals$tjs^2*(sum(vals$tis^2)/vals$total^2)/((vals$nc-1)*my_rowvar-vals$tjs)
	max_size <- 10*max(size);
	size[size < 0] <- max_size;
	#return(list(var_obs=my_rowvar, mus=vals$tjs/vals$nc, sizes=size, vals=vals))
	return(data.frame(mus=vals$tjs/vals$nc, sizes=size))
}
