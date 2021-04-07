combine.integrative.features <- function(
	irf.matrix,
	cis.matrix,
	irf.rescale.recenter = NA,
	cis.rescale.recenter = NA,
	irf.rescale.denominator = NA,
	cis.rescale.denominator = NA,
	irf.weights = rep(1, ncol(irf.matrix)),
	cis.weights = rep(1, ncol(cis.matrix))
	) {

	# check that the matrices have the same rows
	if (nrow(irf.matrix) != nrow(cis.matrix)) {
		stop('irf.matrix and cis.matrix do not have the same number of patients');
		}
	if (all(rownames(irf.matrix) != rownames(cis.matrix))) {
		stop('irf.matrix and cis.matrix rownames do not match');
		}

	# check for the right length and adjust by each parameter
	if (length(irf.rescale.recenter) != 1 || !is.na(irf.rescale.recenter)) {
		if (length(irf.rescale.recenter) == 1 && irf.rescale.recenter == 'mean') {
			irf.rescale.recenter <- apply(irf.matrix, 2, mean);
			}
		if (length(irf.rescale.recenter) == 1) {
			irf.matrix <- irf.matrix - irf.rescale.recenter;
			}
		else if (length(irf.rescale.recenter) == ncol(irf.matrix)) {
			for(j in 1:ncol(irf.matrix)) {
				irf.matrix[,j] <- irf.matrix[,j] - irf.rescale.recenter[j];
				}
			}
		else {
			warning('irf.rescale.recenter is not a valid length so no recentering was completed for irf.matrix');
			}
		}
	if (length(cis.rescale.recenter) != 1 || !is.na(cis.rescale.recenter)) {
		if (length(cis.rescale.recenter) == 1 && cis.rescale.recenter == 'mean') {
			cis.rescale.recenter <- apply(cis.matrix, 2, mean);
			}
		if (length(cis.rescale.recenter) == 1) {
			cis.matrix <- cis.matrix - cis.rescale.recenter;
			}
		else if (length(cis.rescale.recenter) == ncol(cis.matrix)) {
			for(j in 1:ncol(cis.matrix)) {
				cis.matrix[,j] <- cis.matrix[,j] - cis.rescale.recenter[j];
				}
			}
		else {
			warning('cis.rescale.recenter is not a valid length so no recentering was completed for cis.matrix');
			}
		}
	if (length(irf.rescale.denominator) != 1 || !is.na(irf.rescale.denominator)) {
		if (length(irf.rescale.denominator) == 1 && irf.rescale.denominator == 'sd') {
			irf.rescale.denominator <- apply(irf.matrix, 2, sd);
			}
		if (length(irf.rescale.denominator) == 1) {
			irf.matrix <- irf.matrix / irf.rescale.denominator;
			}
		else if (length(irf.rescale.denominator) == ncol(irf.matrix)) {
			if (any(is.na(irf.rescale.denominator))) {
				warning('irf.rescale.denominator contains NAs');
				}
			if (any(irf.rescale.denominator == 0)) {
				warning('irf.rescale.denominator contains zeros');
				}
			for(j in 1:ncol(irf.matrix)) {
				irf.matrix[,j] <- irf.matrix[,j] / irf.rescale.denominator[j];
				}
			}
		else {
			warning('irf.rescale.denominator is not a valid length so no rescaling was completed for irf.matrix');
			}
		}
	if (length(cis.rescale.denominator) != 1 || !is.na(cis.rescale.denominator)) {
		if (length(cis.rescale.denominator) == 1 && cis.rescale.denominator == 'sd') {
			cis.rescale.denominator <- apply(cis.matrix, 2, sd);
			}
		if (length(cis.rescale.denominator) == 1) {
			cis.matrix <- cis.matrix / cis.rescale.denominator;
			}
		else if (length(cis.rescale.denominator) == ncol(cis.matrix)) {
			if (any(is.na(cis.rescale.denominator))) {
				warning('cis.rescale.denominator contains NAs');
				}
			if (any(cis.rescale.denominator == 0)) {
				warning('cis.rescale.denominator contains zeros');
				}
			for(j in 1:ncol(cis.matrix)) {
				cis.matrix[,j] <- cis.matrix[,j] / cis.rescale.denominator[j];
				}
			}
		else {
			warning('cis.rescale.denominator is not a valid length so no rescaling was completed for cis.matrix');
			}
		}
	if (length(irf.weights) != 1 || !is.na(irf.weights)) {
		if (length(irf.weights) == 1) {
			irf.weights <- rep(irf.weights, ncol(irf.matrix));
			}
		if (length(irf.weights) == ncol(irf.matrix)) {
			for(j in 1:ncol(irf.matrix)) {
				irf.matrix[,j] <- irf.matrix[,j] * irf.weights[j];
				}
			}
		else {
			warning('irf.weights is not a valid length so a weight of 1 was used compared to cis.weights');
			}
		}
	if (length(cis.weights) != 1 || !is.na(cis.weights)) {
		if (length(cis.weights) == 1) {
			cis.weights <- rep(cis.weights, ncol(cis.matrix));
			}
		if (length(cis.weights) == ncol(cis.matrix)) {
			for(j in 1:ncol(cis.matrix)) {
				cis.matrix[,j] <- cis.matrix[,j] * cis.weights[j];
				}
			}
		else {
			warning('cis.weights is not a valid length so a weight of 1 was used compared to irf.weights');
			}
		}

	# combine the matrices into one matrix
	integrative.feature.matrix <- cbind(irf.matrix, cis.matrix);

	# return the re-scaling and re-weighting parameter along with the integrative.feature.matrix for the matrix can be reproduced
	return(list(
		integrative.feature.matrix = integrative.feature.matrix,
		irf.rescale.recenter = irf.rescale.recenter,
		cis.rescale.recenter = cis.rescale.recenter,
		irf.rescale.denominator = irf.rescale.denominator,
		cis.rescale.denominator = cis.rescale.denominator,
		irf.weights = irf.weights,
		cis.weights = cis.weights
		));
	}
