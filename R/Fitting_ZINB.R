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

fit_negative_binomial <- function(counts) {
	vals <- calculate_summary_values(counts);
#	mus <- (vals$tjs) %*% t(vals$tis/vals$total)
	
	min_size <- 10^-10;
	my_rowvar <- sapply(1:nrow(counts), function(i){
				mu_is <- vals$tjs[i]*vals$tis/vals$total
				var(as.vector(unlist(counts[i,]))-mu_is)
			})
		
	size <- vals$tjs^2*(sum(vals$tis^2)/vals$total^2)/((vals$nc-1)*my_rowvar-vals$tjs) # for this to work with sparse matrices might need to implement in C
	max_size <- 10*max(size);
	size[size < 0] <- max_size;
	size[size < min_size] <- min_size;
#	size[size > max_size] <- max_size;

	return(list(var_obs=my_rowvar, sizes=size, vals=vals))
}

calc_expected_zeros <- function(fit) {
	vals <- fit$vals;
	droprate_exp <- vector(length=vals$ng)
	exp_size <- fit$sizes

	for (i in 1:vals$ng) {
                mu_is <- vals$tjs[i]*vals$tis/vals$total
                p_is <- (1+mu_is/exp_size[i])^(-exp_size[i]);
                #p_var_is <- p_is*(1-p_is);
                droprate_exp[i] <- sum(p_is)/vals$nc;
		#droprate_exp_err[i] <- sqrt(sum(p_var_is)/(vals$nc^2));
        }
	return(droprate_exp);
}

fit_zero_inflated_negative_binomial <- function() {

}

