combine.compressed.features <- function(
	irf,
	cis,
	irf.rescale.recenter = NA,
	cis.rescale.recenter = NA,
	irf.rescale.denominator = NA,
	cis.rescale.denominator = NA,
	irf.weights = rep(1, ncol(irf)),
	cis.weights = rep(1, ncol(cis))
	) { 

	# check that the matrices have the same rows
	if(nrow(irf) != nrow(cis)) {
		stop('irf and cis matrices do not have the same number of patients');
		}
	if(all(rownames(irf) != rownames(cis))) {
		stop('irf and cis matrices rownames do not match');
		}
	# check for the right length and adjust by each parameter
	if(length(irf.rescale.recenter) != 1 || !is.na(irf.rescale.recenter)) {
		if(length(irf.rescale.recenter) == 1 && irf.rescale.recenter == 'mean') {
			irf.rescale.recenter <- apply(irf, 2, mean);
			}
		if(length(irf.rescale.recenter) == 1) {
			irf <- irf - irf.rescale.recenter;
			}
		else if(length(irf.rescale.recenter) == ncol(irf)) {
			for(j in ncol(irf)) {
				irf[,j] <- irf[,j] - irf.rescale.recenter[j];
				}
			}
		else {
			warning('irf.rescale.recenter is not a valid length so no recentering was completed for irf');
			}
		}
	if(length(cis.rescale.recenter) != 1 || !is.na(cis.rescale.recenter)) {
		if(length(cis.rescale.recenter) == 1 && cis.rescale.recenter == 'mean') {
			cis.rescale.recenter <- apply(cis, 2, mean);
			}
		if(length(cis.rescale.recenter) == 1) {
			cis <- cis - cis.rescale.recenter;
			}
		else if(length(cis.rescale.recenter) == ncol(cis)) {
			for(j in ncol(cis)) {
				cis[,j] <- cis[,j] - cis.rescale.recenter[j];
				}
			}
		else {
			warning('cis.rescale.recenter is not a valid length so no recentering was completed for cis');
			}
		}
	if(length(irf.rescale.denominator) != 1 || !is.na(irf.rescale.denominator)) {
		if(length(irf.rescale.denominator) == 1 && irf.rescale.denominator == 'sd') {
			irf.rescale.denominator <- apply(irf, 2, sd);
			}
		if(length(irf.rescale.denominator) == 1) {
			irf <- irf/irf.rescale.denominator;
			}
		else if(length(irf.rescale.denominator) == ncol(irf)) {
			for(j in ncol(irf)) {
				irf[,j] <- irf[,j]/irf.rescale.denominator[j];
				}
			}
		else {
			warning('irf.rescale.denominator is not a valid length so no rescaling was completed for irf');
			}
		}
	if(length(cis.rescale.denominator) != 1 || !is.na(cis.rescale.denominator)) {
		if(length(cis.rescale.denominator) == 1 && cis.rescale.denominator == 'sd') {
			cis.rescale.denominator <- apply(cis, 2, sd);
			}
		if(length(cis.rescale.denominator) == 1) {
			cis <- cis/cis.rescale.denominator;
			}
		else if(length(cis.rescale.denominator) == ncol(cis)) {
			for(j in ncol(cis)) {
				cis[,j] <- cis[,j]/cis.rescale.denominator[j];
				}
			}
		else {
			warning('cis.rescale.denominator is not a valid length so no rescaling was completed for cis');
			}
		}
	if(length(irf.weights) != 1 || !is.na(irf.weights)) {
		if(length(irf.weights) == 1) {
			irf.weights <- rep(irf.weights, ncol(irf));
			}
		if(length(irf.weights) == ncol(irf)) {
			for(j in ncol(irf)) {
				irf[,j] <- irf[,j]*irf.weights[j];
				}
			}
		else {
			warning('irf.weights is not a valid length so a weight of 1 was used compared to cis.weights');
			}
		}
	if(length(cis.weights) != 1 || !is.na(cis.weights)) {
		if(length(cis.weights) == 1) {
			cis.weights <- rep(cis.weights, ncol(cis));
			}
		if(length(cis.weights) == ncol(cis)) {
			for(j in ncol(cis)) {
				cis[,j] <- cis[,j]*cis.weights[j];
				}
			}
		else {
			warning('cis.weights is not a valid length so a weight of 1 was used compared to irf.weights');
			}
		}

	# combine the matrices into one matrix
	compressed.feature.matrix <- cbind(irf, cis);

	# return the re-scaling and re-weighting parameter along with the compressed.feature.matrix for the matrix can be reproduced
	return(list(
		compressed.feature.matrix = compressed.feature.matrix,
		irf.rescale.recenter = irf.rescale.recenter,
		cis.rescale.recenter = cis.rescale.recenter,
		irf.rescale.denominator = irf.rescale.denominator,
		cis.rescale.denominator = cis.rescale.denominator,
		irf.weights = irf.weights,
		cis.weights = cis.weights
		));
	}
