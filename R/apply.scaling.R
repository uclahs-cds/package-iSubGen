apply.scaling <- function(data.matrices, scaling.factors) {

	# data.matrices can be a single matrix or a list of matrices
	# if the data is a single matrix then the class will be 'matrix'	
	if (class(data.matrices)[1] == 'matrix') {
	
		# check that scaling.factors is have elements with the names "center" and "scale"
		if (all(!c('center','scale') %in% names(scaling.factors))) {
			stop('for each data matrix, scaling.factor needs to be a list with center and scale ');
			}


		# check that scaling.factors are the correct format
		if (length(scaling.factors$center) != nrow(data.matrices)) {
			stop('the length of scaling.factors$center needs to match the number of rows in data.matrices');
			}
		if (length(scaling.factors$scale) != nrow(data.matrices)) {
			stop('the length of scaling.factors$scale needs to match the number of rows in data.matrices');
			}

		# if necessary adjust the format of the scaling factors for a single matrix
		if (class(scaling.factors$center) == 'list') {
			scaling.factors$center <- scaling.factors$center[[1]];
			warning('the first item from the scaling.factor$center list was used for scaling');
			}
		if (class(scaling.factors$scale) == 'list') {
			scaling.factors$scale <- scaling.factors$scale[[1]];
			warning('the first item from the scaling.factor$scale list was used for scaling');
			}

		# scale each row in the matrix by the corresponding scaling factors
		for(i in 1:nrow(data.matrices)) {
			center.adjustment <- scaling.factors$center[rownames(data.matrices)[i]];
			scale.adjustment <- 1;
			if(scaling.factors$scale[rownames(data.matrices)[i]] > 0) {
				scale.adjustment <- scaling.factors$scale[rownames(data.matrices)[i]];
				}
			data.matrices[i,] <- (data.matrices[i,] - center.adjustment) / scale.adjustment;
			}

		# return the scaled single matrix
		return(data.matrices);
		}
	else if (class(data.matrices) == 'list') {
		# if you make it to this point then data.matrices is a list
		# so check the format of the input and then recurse on each matrix

		# check that scaling.factors are the correct format
		if (any(sort(names(data.matrices)) != sort(names(scaling.factors)))) {
			stop('the scaling.factors list needs to have the same names as the data.matrices list');
			}

		# if you get to this point then data.matrices is a list of matrices
		for(data.type in names(data.matrices)) {

			# check that scaling.factors are the correct format
			if (length(scaling.factors[[data.type]]$center) != nrow(data.matrices[[data.type]])) {
				stop(paste0('scaling.factors$',data.type,'$center does not match the number of rows in data.matrices$',data.type));
				}
			if (length(scaling.factors[[data.type]]$scale) != nrow(data.matrices[[data.type]])) {
				stop(paste0('scaling.factors$',data.type,'$scale does not match the number of rows in data.matrices$',data.type));
				}
		
			# call the function for each data type
			data.matrices[[data.type]] <- apply.scaling(data.matrices[[data.type]],scaling.factors[[data.type]]);
			}

		# return the scaled list of matrices
		return(data.matrices);
		}
	stop('data.matrices needs to be a matrix or a list of matrices');
	}

