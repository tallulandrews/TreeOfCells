check_profiles <- function(profiles) {
}

check_merged_profiles <- function(merged_profiles) {
	bits <- c("mus", "rs", "ds", "Ns");
	if (sum(names(merged_profiles) %in% bits) != 4) {
		stop("Merged Profiles does not contain correct matrices")
	}  
	dimension <- dim(merged_profiles$mus);
	if (! (identical(dim(merged_profiles$rs), dimension) & 
		identical(dim(merged_profiles$ds), dimension) &
		identical(dim(merged_profiles$Ns), dimension))) {
		stop("Merged profile matrices are not all the same size.")
	}
	# Column names
	c_names <- as.character(colnames(merged_profiles$mus))
	if (length(c_names) != ncol(merged_profiles$mus)){
		warning("No cell-type names, I will create some!")
		c_names <- paste("cell_type_", 1:ncol(merged_profiles$mus), sep="")
		colnames(merged_profiles$mus) <- c_names;
	}
	if (! identical(as.character(colnames(merged_profiles$rs)), c_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		colnames(merged_profiles$rs) <- c_names
	}
	if (! identical(as.character(colnames(merged_profiles$ds)), c_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		colnames(merged_profiles$ds) <- c_names
	}
	if (! identical(as.character(colnames(merged_profiles$Ns)), c_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		colnames(merged_profiles$Ns) <- c_names
	}

	# row names
	r_names <- as.character(rownames(merged_profiles$mus))
	if (length(r_names) != nrow(merged_profiles$mus)){
		warning("No gene names, I will create some!")
		r_names <- paste("gene_", 1:nrow(merged_profiles$mus), sep="")
		rownames(merged_profiles$mus) <- r_names;
	}
	if (! identical(as.character(rownames(merged_profiles$rs)), r_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		rownames(merged_profiles$rs) <- r_names
	}
	if (! identical(as.character(rownames(merged_profiles$ds)), r_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		rownames(merged_profiles$ds) <- r_names
	}
	if (! identical(as.character(rownames(merged_profiles$Ns)), r_names)) {
		warning("Disagreeing cell-type names, defaulting to those of mus")
		rownames(merged_profiles$Ns) <- r_names
	}

	return(merged_profiles);

}

check_row_col_names <- function(mat) {
	r_names <- as.character(rownames(mat))
        if (length(r_names) != nrow(mat)){
                warning("No gene names, I will create some!")
                r_names <- paste("gene_", 1:nrow(mat), sep="")
                rownames(mat) <- r_names;
        }
	
	c_names <- as.character(colnames(mat))
        if (length(c_names) != ncol(mat)){
                warning("No cell-type names, I will create some!")
                c_names <- paste("cell_type_", 1:ncol(mat), sep="")
                colnames(mat) <- c_names;
        }
	return(mat);
}
