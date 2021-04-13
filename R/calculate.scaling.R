calculate.scaling <- function(data.matrices) {
	# if there is only one data type to scale
	if (class(data.matrices)[1] == 'matrix') {
		# return the mean and sd of each row
		return(list(
			center = apply(data.matrices,1,mean),
			scale = apply(data.matrices,1,sd)
			));
		}
	if (class(data.matrices) == 'list') {
		# if there are multiple data types to scale
		# return the mean and sd of each row for each data matrix
		scaling.factors <- list();
		for(data.type in names(data.matrices)) {
			scaling.factors[[data.type]] <- calculate.scaling(data.matrices[[data.type]]);
			}
		return(scaling.factors);
		}
	# if not a list or a matrix return an error message to let the user know how to correct the input
	stop('data.matrices needs to be a matrix or a list');
	}
