read.scaling.factors <- function(scaling.factor.files.dir, data.types) {
	scaling.factors <- list(center = list(), scale = list());
	for(data.type in data.types) {
		scaling.matrix <- read.table(
			paste0(scaling.factor.files.dir, data.type, '_gene_recenter.txt'),
			sep = '\t',
			header = TRUE,
			check.names = FALSE
			);
		scaling.factors$center[[data.type]] <- as.numeric(scaling.matrix[1,]);
		names(scaling.factors$center[[data.type]]) <- colnames(scaling.matrix);

		scaling.matrix <- read.table(
			paste0(scaling.factor.files.dir, data.type, '_gene_rescale.txt'),
			sep = '\t',
			header = TRUE,
			check.names = FALSE
			);
		scaling.factors$scale[[data.type]] <- as.numeric(scaling.matrix[1,]);
		names(scaling.factors$scale[[data.type]]) <- colnames(scaling.matrix);
		}
	return(scaling.factors);
	}
