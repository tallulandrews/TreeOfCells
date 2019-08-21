merge_lists <- function(list1, list2, name1="A", name2="B") {
	if (sum(names(list1) %in% names(list2)) > 0) {
		names(list1) <- paste(name1, names(list1), sep="_")
		names(list2) <- paste(name2, names(list2), sep="_")
	}
	return(c(list1, list2))
}


merge_profiles <- function(list_of_profiles) {
	mus <- vector();
	ds <- vector();
	rs <- vector();
	Ns <- vector();
	for (p in list_of_profiles) {
		if (length(mus) == 0) {
			mus <- matrix(p[,"mu"], ncol=1)
			rs <- matrix(p[,"r"], ncol=1)
			ds <- matrix(p[,"d"], ncol=1)
			Ns <- matrix(p[,"N"], ncol=1)
			rownames(mus) <- rownames(p)
			rownames(ds) <- rownames(p)
			rownames(rs) <- rownames(p)
			rownames(Ns) <- rownames(p)
		} else {
			if (!identical(rownames(p) , rownames(mus))){
				all_g <- sort(union(rownames(p), rownames(mus)))
				mus <- mus[match(all_g, rownames(mus)),]
				rs <- rs[match(all_g, rownames(rs)),]
				ds <- ds[match(all_g, rownames(ds)),]
				Ns <- Ns[match(all_g, rownames(Ns)),]
				p <- p[match(all_g, rownames(p)),]
				rownames(mus) <- all_g
				rownames(ds) <- all_g
				rownames(rs) <- all_g
				rownames(Ns) <- all_g
				rownames(p) <- all_g
			} 
			mus <- cbind(mus, p[,"mu"])
			rs <- cbind(rs, p[,"r"])
			ds <- cbind(ds, p[,"d"])
			Ns <- cbind(Ns, p[,"N"])
		}
	}
	mus[is.na(mus)] <- 0
	rs[is.na(rs)] <- 0
	ds[is.na(ds)] <- 0
	Ns[is.na(Ns)] <- 0
	colnames(mus) <- names(list_of_profiles)
	colnames(Ns) <- names(list_of_profiles)
	colnames(rs) <- names(list_of_profiles)
	colnames(ds) <- names(list_of_profiles)
	return(list(mus=mus, rs=rs, ds=ds, Ns=Ns));
}
