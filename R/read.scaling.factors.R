read.scaling.factors <- function(scaling.factor.files.dir, aberration.types) {
	scaling.factors <- list(center=list(),scale=list());
	for(aberration.type in aberration.types) {
		scaling.factors$center[[aberration.type]] <- read.table(
			paste0(scaling.factor.files.dir,aberration.type,'_gene_recenter.txt')
			);
		scaling.factors$scale[[aberration.type]] <- read.table(
			paste0(scaling.factor.files.dir,aberration.type,'_gene_rescale.txt')
			);
	}
	return(scaling.factors);
}
