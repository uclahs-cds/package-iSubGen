write.scaling.factors <- function(scaling.factors,scaling.factor.files.dir) {
	for(x in c('center','scale')) {
		for(aberration.type in names(scaling.factors[[x]])) {
			write.table(
				scaling.factors[[x]][[aberration.type]],
				paste0(scaling.factor.files.dir,aberration.type,'_gene_re',x,'.txt'),
				quote=FALSE,
				sep='\t'
				);
		}
	}
}
