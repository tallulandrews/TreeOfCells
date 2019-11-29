\name{merge_lists}
\alias{merge_lists}
\alias{merge_profiles}
\alias{merge_across_protocols}
\title{Merge Profiles}
\description{
	Methods for merging profiles across cell-types, and experimental protocols.
}
\usage{
	merge_across_protocols(merged_profile_list, norm_func=scaled_regression, fs_func="none")
	merge_profiles(list_of_profiles)
	merge_lists(list1, list2, name1="A", name2="B")
}
\arguments{
	\item{merged_profile_list}{a named list of output from merge_profiles, one item for each different experimental protocol (e.g. Chromium, Smartseq, Cellseq).}
	\item{list_of_profiles}{a named list of output from ZINB fitting, one item for each cell-type.}
	\item{norm_func}{normalization function to use}
	\item{fs_func}{feature selection function to use}
	\item{list1}{first list}
	\item{list2}{second list}
	\item{name1}{the name of the first list, used as a prefix to avoid duplicate items.}
	\item{name2}{the name of the second list, used as a prefix to avoid duplicate items.}
}
\details{
	\code{merge_lists} merges together two lists ensuring all contents of each list are present in the final list.
	\code{merge_profiles} reorganizes fitted ZINB parameters into a table of mus, rs, and ds where the columns are each cell-type, and rows are each gene.

	\code{merge_across_protocols} performs normalization and feature selection on profiles from each protocol separately then merges them together into a single set of mus, rs, and ds matrices. Reported features are detected in at least two of the protocols. 
}
\value{
	\code{merge_lists} a list of containing all contents of each input list.
	\code{merge_profiles} a list of mus, rs, ds and Ns, the parameters of the ZINB across all genes (rows) and cell-types (columns), Ns are the number of cells in each cell-type.
	\code{merge_across_protocols} a list of two items: 
		merged_profiles (list of parameters across all cell-types & protocols)
		features (consensus features from at least 2 protocols).
}
\examples{
	set.seed(1010)
	counts <- matrix(rnbinom(5000*100, mu=30, size=0.5), ncol=100);
	fits1 <- fit_NB_to_matrix(counts[,1:50]);
	fits2 <- fit_NB_to_matrix(counts[,51:75]);
	fits3 <- fit_NB_to_matrix(counts[,76:100]);
	merged_profiles <- merge_profiles(list(celltype1=fits1, celltype2=fits2, celltype3=fits3));
	cross_protocol <- merge_across_protocols(list(X10=merged_profiles, SS2=merged_profiles))

	l1 <- list("hat"=1, "cat"=2, "pat"=3)
	l2 <- list("dog"=3, "cat"=10, "horse"=4)
	merge_l <- merge_lists(l1, l2, name1="rhyme", name2="animal")
}
\keyword{profiles, merging}
