calculate.scaling <- function(data.matrices) {
	# if there is only one data type to scale
	if ('matrix' == class(data.matrices)[1] || 'data.frame' == class(data.matrices)[1]) {
		# return the mean and sd of each row
		return(list(
			center = apply(data.matrices, 1, mean),
			scale = apply(data.matrices, 1, sd)
			))
	if ('list' == class(data.matrices)) {
		# if there are multiple data types to scale
		# return the mean and sd of each row for each data matrix
		scaling.factors <- list();
		for(data.type in names(data.matrices)) {
			scaling.factors[[data.type]] <- calculate.scaling(data.matrices[[data.type]]);
			}
		return(scaling.factors);
		}
	# if `data.matrices` is not a list, a matrix, or a data.frame return an error message to let the user know how to correct the input
	stop('`data.matrices` must be a matrix, a data frame or a list');
	}
}