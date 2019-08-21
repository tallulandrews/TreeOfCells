get_path_dist <- function(node, field) {
	path <- vector();
	tmp <- node;
	for(i in 1:node$level){
		old_names <- names(path);
		if (is.null(tmp[[field]])) {
			tmp[[field]] <- 0;
			warning(paste("Warning: Missing distance for", tmp$name, "filling in with 0."))
		}
		path <- c(path, tmp[[field]])
		names(path) <- c(old_names, tmp$name);
		tmp <- tmp$parent;
	}
	return(rev(path));
}

tree_weighted_Distance <- function (node1, node2, field="dist"){
	# modified from data.tree package
    if (!identical(node1$root, node2$root))
        stop("node1 and node2 must be in same tree!")
    path1 <- get_path_dist(node1, field);
    path2 <- get_path_dist(node2, field);
	
    i <- 1
    maxi <- min(node1$level, node2$level)
    common_dist <- 0
    while (names(path1)[i] == names(path2)[i] && i <= maxi) {
	common_dist <- common_dist + path1[i]
	i <- i + 1; 
	}
    distance <- sum(path1) + sum(path2) - 2 * common_dist
    return(distance)
}

run_Sim_benchmark <- function(test_fun, sim_file, sim_profiles=NULL) {
	require(data.tree)
	dat <- readRDS(sim_file);
	if (is.null(sim_profiles)) {
		outfile <- paste("Profile", sim_file, sep="_");
        	sce <- readRDS(files);
        	profiles <- fit_ZINB_to_SCE(sce, "Group");
        	sim_profiles <- merge_profiles(profiles)
		saveRDS(sim_profiles, outfile);
	}
	
	dist1 <- test_fun(sim_profiles)
	#Truth
	node_dist <- dat@metadata$true_tree$Get("de_dist")
	tree <- dat@metadata$true_tree;
	distance_mat <- matrix(0, nrow=length(node_dist), ncol=length(node_dist));
	for (i in 1:length(node_dist)) {
		node1 <- FindNode(tree, names(node_dist)[i])
		for(j in (i+1):length(node_dist)) {
			node2 <- FindNode(tree, names(node_dist)[j]);
			d <- tree_weighted_Distance(node1, node2, field="de_dist")
			distance_mat[i,j] <- distance_mat[j,i] <- d;
		}
	}

	a <- as.vector(dist1)
	b <- as.vector(dist2)
	reg <- lm(a~b+0)
	return(list(
		disttest = dist1, 
		disttrue = dist2, 
		sim=cor(as.vector(dist1), as.vector(dist2)), 
		diff=sum(abs(rank(a) - rank(b))),
		var=sum(((a-b)^2))/length(a)^2,
		rsq=summary(reg)$r.squared
		))
}
