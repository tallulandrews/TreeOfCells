#### Fitting #####
calculate_summary_values <- function(counts) {
        if (sum(counts < 0) >0) {stop("Expression matrix contains negative values! Please provide raw UMI counts!")}
	if ( sum(counts >= 1) != sum(counts > 0) ) {stop("Error: Expression matrix is not integers! Please provide raw UMI counts.")}
#        if (sum(!is.integer(counts)) >0) {stop("Expression matrix is not integers! Please provide a matrix (not data.frame) raw UMI counts!")}

        tjs <- Matrix::rowSums(counts, na.rm=T) # Total molecules/gene
	if (sum(tjs <= 0) > 0) {stop("Error: all genes must have at least one detected molecule.")}
        tis <- Matrix::colSums(counts, na.rm=T) # Total molecules/cell
	if (sum(tis <= 0) > 0) {stop("Error: all cells must have at least one detected molecule.")}
        djs <- ncol(counts)-Matrix::rowSums(counts > 0, na.rm=T) # Observed Dropouts per gene
        dis <- nrow(counts)-Matrix::colSums(counts > 0, na.rm=T) # Observed Dropouts per cell
        nc <- length(counts[1,]) # Number of cells
        ng <- length(counts[,1]) # Number of genes
        total <- sum(tis, na.rm=T) # Total molecules sampled
        return(list(tis=tis, tjs=tjs, dis=dis, djs=djs, total=total,nc=nc,ng=ng));
}

convert_to_integer <- function(mat) {
        mat <- ceiling(as.matrix(mat))
        storage.mode(mat) <- "integer"
        mat <- mat[Matrix::rowSums(mat) > 0,]
        return(mat)
}


fit_ZINB_to_matrix <- function(counts) {
	zeros <- which(Matrix::rowSums(counts) == 0);
	if (length(zeros) > 0) {
		counts <- counts[-zeros,]
	}

	vals <- calculate_summary_values(counts);
	out <- sapply(1:vals$ng, function(g) {
		fit_zero_inflated_negative_binomial(counts[g,], g, vals)
		})
	counts <- check_row_col_names(counts);
	colnames(out) <- rownames(counts);
	out <- t(out);
	if (length(zeros) > 0) {
		out <- rbind(out, matrix(0, nrow=length(zeros), ncol=ncol(out)))
		rownames(out) <- c(rownames(counts), names(zeros))
	}
#	out <- out[order(rownames(out)),]
	colnames(out) <- c("mu", "r", "d", "N")
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
	return(c(mu_j, r_j, d_prev, vals$nc));
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
	size[size < min_size] <- min_size;
	
	#return(list(var_obs=my_rowvar, mus=vals$tjs/vals$nc, sizes=size, vals=vals))
	return(data.frame(mu=vals$tjs/vals$nc, r=size, d=rep(0, length(size)), N=vals$nc))
}

fit_ZINB_to_SCE <- function(sce, lab_column="cell_type1") {
	if (class(sce) != "SingleCellExperiment") {stop("Error: Input must be an SCE object")}
	if (! ("counts" %in% names(SummarizedExperiment::assays(sce)))) {stop("Error: Input must contain a counts matrix")}
	if (! (lab_column %in% colnames(SingleCellExperiment::colData(sce)))) {stop("Error: cannot find cell-type labels")}

	type_labs <- SingleCellExperiment::colData(sce)[,lab_column]
	OUT_list <- list()
	for (type in levels(factor(type_labs))) {
	        out <- fit_ZINB_to_matrix(SummarizedExperiment::assays(sce)[["counts"]][,type_labs == type & !is.na(type_labs)]);
        	OUT_list[[type]] <- out;
	}
	return(OUT_list)
}
