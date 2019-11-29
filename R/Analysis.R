normalize_feature_selection <- function(merged_profiles, norm_func=scaled_regression, fs_func="none"){
	if (class(norm_func) != "character") {
		merged_profiles <- norm_func(merged_profiles)
	}
	if (class(fs_func) != "character") {
		features <- fs_func(merged_profiles)
	} else {
		features <- rownames(merged_profiles$mus);
	}
	return(list(merged_profiles=merged_profiles, FS=features));
}

merge_across_protocols <- function(merged_profile_list, norm_func=scaled_regression, fs_func="none") {
	# Get common genes across all & features IDed in at least 2
	features <- vector();
	common_genes <- vector();
	for (i in 1:length(merged_profile_list)) {
		norm_fs <- normalize_feature_selection(merged_profile_list[[i]], norm_func=norm_func, fs_func=fs_func)
		merged_profile_list[[i]] <- norm_fs$merged_profiles;
		this_genes <- as.character(rownames(norm_fs$merged_profiles$mus))
		if (i == 1) {
			common_genes <- this_genes;
			features <- as.character(norm_fs$FS)
		} else {
			common_genes <- intersect(as.character(common_genes), this_genes)
			features <- c(features, as.character(norm_fs$FS));
		}
	}
	features <- unique(features[duplicated(features)]);

	# Align datasets 
	common_genes <- sort(common_genes)
	MERGE <- list();
	for (i in 1:length(merged_profile_list)) {
		this_genes <- as.character(rownames(norm_fs$merged_profiles$mus))
		profiles <- merged_profile_list[[i]]
		reorder <- match(common_genes, rownames(profiles$mus));
		for (j in 1:length(profiles)) {
                	profiles[[j]] <- profiles[[j]][reorder,]
			colnames(profiles[[j]]) <- paste(names(merged_profile_list)[j],colnames(profiles[[j]]), sep="");

			if (i == 1) {
				MERGE[[j]] <- profiles[[j]];
			} else {
				MERGE[[j]] <- cbind(MERGE[[j]], profiles[[j]]);
			}
		}
		names(MERGE) <- names(profiles);
        }

	return(list(merged_profiles=MERGE, features=features));
}

tree_of_cells <- function( merged_profiles, gene_features=rownames(merged_profiles$mus), dist_func=pearson_dist, plot_type="radial") {

		
	merged_profiles <- apply_feature_selection(merged_profiles, gene_features)	
	dist_mat <- dist_func(merged_profiles)
	tree <- phangorn::NJ(dist_mat);
	plot(tree, cex=0.6, type=plot_type)
	return(list(distances=dist_mat, tree=tree))
}
