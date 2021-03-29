write.scaling.factors <- function(scaling.factors, scaling.factor.files.dir) {
	# write two files for each data type, one for the center values and one for the scale values
	for(x in c('center','scale')) {
		for(data.type in names(scaling.factors[[x]])) {
			write.table(
				t(as.data.frame(scaling.factors[[x]][[data.type]])),
				paste0(scaling.factor.files.dir,data.type,'_gene_re',x,'.txt'),
				quote = FALSE,
				sep = '\t',
				col.names = TRUE,
				row.names = FALSE
				);
			}
		}
	}
