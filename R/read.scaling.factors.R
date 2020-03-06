read.scaling.factors <- function(scaling.factor.files.dir, aberration.types) {
	scaling.factors <- list(center=list(),scale=list());
	for(aberration.type in aberration.types) {
		scaling.matrix <- read.table(
			paste0(scaling.factor.files.dir,aberration.type,'_gene_recenter.txt'),
			sep='\t',
			header=TRUE,
			check.names=FALSE
			);
		scaling.factors$center[[aberration.type]] <- as.numeric(scaling.matrix[1,]);
		names(scaling.factors$center[[aberration.type]]) <- colnames(scaling.matrix);

		scaling.matrix <- read.table(
			paste0(scaling.factor.files.dir,aberration.type,'_gene_rescale.txt'),
			sep='\t',
			header=TRUE,
			check.names=FALSE
			);
		scaling.factors$scale[[aberration.type]] <- as.numeric(scaling.matrix[1,]);
		names(scaling.factors$scale[[aberration.type]]) <- colnames(scaling.matrix);
	}
	return(scaling.factors);
}
