read.scaling.factors <- function(scaling.factor.files.dir, data.types) {
	# read two scaling files for each data type
	scaling.factors <- list();
	for(data.type in data.types) {
		scaling.factors[[data.type]] <- list();

		# re-centering values
		filename <- paste0(scaling.factor.files.dir, data.type, '_gene_recenter.txt');
		if(exists(filename)) {
			scaling.matrix <- read.table(
				filename,
				sep = '\t',
				header = TRUE,
				check.names = FALSE
				);
			scaling.factors[[data.type]]$center <- as.numeric(scaling.matrix[1,]);
			names(scaling.factors[[data.type]]$center) <- colnames(scaling.matrix);
			}
		else {
			warning(paste0("cannot open file '",filename,"': No such file or directory"));
			}

		# re-scaling values 
		filename <- paste0(scaling.factor.files.dir, data.type, '_gene_rescale.txt');
		if(exists(filename)) {
			scaling.matrix <- read.table(
				filename,
				sep = '\t',
				header = TRUE,
				check.names = FALSE
				);
			scaling.factors[[data.type]]$scale <- as.numeric(scaling.matrix[1,]);
			names(scaling.factors[[data.type]]$scale) <- colnames(scaling.matrix);
			}
		else {
			warning(paste0("cannot open file '",filename,"': No such file or directory"));
			}
		}
	return(scaling.factors);
	}
