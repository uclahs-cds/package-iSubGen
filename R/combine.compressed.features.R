combine.compressed.features <- function(
	reduced.features,
	pairwise.similiarity.features,
	reduced.features.rescale.recenter = NA,
	pairwise.similiarity.features.rescale.recenter = NA,
	reduced.features.rescale.denominator = NA,
	pairwise.similiarity.features.rescale.denominator = NA,
	reduced.features.weights = rep(1,ncol(reduced.features)),
	pairwise.similiarity.features.weights = rep(1,ncol(pairwise.similiarity.features))
	) { 

	# check that the matrices have the same rows
	if(nrow(reduced.features) != nrow(pairwise.similiarity.features)) {
		stop('reduced.features and pairwise.similiarity.features matrices do not have the same number of patients');
	}
	if(all(rownames(reduced.features) != rownames(pairwise.similiarity.features))) {
		stop('reduced.features and pairwise.similiarity.features matrices rownames do not match');
	}
	# check for the right length and adjust by each parameter
	if(!is.na(reduced.features.rescale.recenter)) {
		if(length(reduced.features.rescale.recenter) == 1 && reduced.features.rescale.recenter == 'mean') {
			reduced.features.rescale.recenter <- apply(reduced.features.rescale.recenter,2,mean);
		}
		if(length(reduced.features.rescale.recenter) == 1) {
			reduced.features <- reduced.features - reduced.features.rescale.recenter;
		} else if(length(reduced.features.rescale.recenter) == ncol(reduced.features)) {
			for(j in ncol(reduced.features)) {
				reduced.features[,j] <- reduced.features[,j] - reduced.features.rescale.recenter[j];
			}
		} else {
			warning('reduced.features.rescale.recenter is not a valid length so no recentering was completed for reduced.features');
		}
	}
	if(!is.na(pairwise.similiarity.features.rescale.recenter)) {
		if(length(pairwise.similiarity.features.rescale.recenter) == 1 && pairwise.similiarity.features.rescale.recenter == 'mean') {
			pairwise.similiarity.features.rescale.recenter <- apply(pairwise.similiarity.features.rescale.recenter,2,mean);
		}
		if(length(pairwise.similiarity.features.rescale.recenter) == 1) {
			pairwise.similiarity.features <- pairwise.similiarity.features - pairwise.similiarity.features.rescale.recenter;
		} else if(length(pairwise.similiarity.features.rescale.recenter) == ncol(pairwise.similiarity.features)) {
			for(j in ncol(pairwise.similiarity.features)) {
				pairwise.similiarity.features[,j] <- pairwise.similiarity.features[,j] - pairwise.similiarity.features.rescale.recenter[j];
			}
		} else {
			warning('pairwise.similiarity.features.rescale.recenter is not a valid length so no recentering was completed for pairwise.similiarity.features');
		}
	}
	if(!is.na(reduced.features.rescale.denominator)) {
		if(length(reduced.features.rescale.denominator) == 1 && reduced.features.rescale.denominator == 'sd') {
			reduced.features.rescale.recenter <- apply(reduced.features.rescale.recenter,2,sd);
		}
		if(length(reduced.features.rescale.denominator) == 1) {
			reduced.features <- reduced.features/reduced.features.rescale.denominator;
		} else if(length(reduced.features.rescale.denominator) == ncol(reduced.features)) {
			for(j in ncol(reduced.features)) {
				reduced.features[,j] <- reduced.features[,j]/reduced.features.rescale.denominator[j];
			}
		} else {
			warning('reduced.features.rescale.denominator is not a valid length so no rescaling was completed for reduced.features');
		}
	}
	if(!is.na(pairwise.similiarity.features.rescale.denominator)) {
		if(length(pairwise.similiarity.features.rescale.denominator) == 1 && pairwise.similiarity.features.rescale.denominator == 'sd') {
			pairwise.similiarity.features.rescale.recenter <- apply(pairwise.similiarity.features.rescale.recenter,2,sd);
		}
		if(length(pairwise.similiarity.features.rescale.denominator) == 1) {
			pairwise.similiarity.features <- pairwise.similiarity.features/pairwise.similiarity.features.rescale.denominator;
		} else if(length(pairwise.similiarity.features.rescale.denominator) == ncol(pairwise.similiarity.features)) {
			for(j in ncol(pairwise.similiarity.features)) {
				pairwise.similiarity.features[,j] <- pairwise.similiarity.features[,j]/pairwise.similiarity.features.rescale.denominator[j];
			}
		} else {
			warning('pairwise.similiarity.features.rescale.denominator is not a valid length so no rescaling was completed for pairwise.similiarity.features');
		}
	}
	if(!is.na(reduced.features.weights)) {
		if(length(reduced.features.weights) == 1) {
			reduced.features.weights <- rep(reduced.features.weights,ncol(reduced.features));
		}
		if(length(reduced.features.weights) == ncol(reduced.features)) {
			for(j in ncol(reduced.features)) {
				reduced.features[,j] <- reduced.features[,j]*reduced.features.weights[j];
			}
		} else {
			warning('reduced.features.weights is not a valid length so a weight of 1 was used compared to pairwise.similiarity.features.weights');
		}
	}
	if(!is.na(pairwise.similiarity.features.weights)) {
		if(length(pairwise.similiarity.features.weights) == 1) {
			pairwise.similiarity.features.weights <- rep(pairwise.similiarity.features.weights,ncol(pairwise.similiarity.features));
		}
		if(length(pairwise.similiarity.features.weights) == ncol(pairwise.similiarity.features)) {
			for(j in ncol(pairwise.similiarity.features)) {
				pairwise.similiarity.features[,j] <- pairwise.similiarity.features[,j]*pairwise.similiarity.features.weights[j];
			}
		} else {
			warning('pairwise.similiarity.features.weights is not a valid length so a weight of 1 was used compared to reduced.features.weights');
		}
	}

	# combine the matrices into one matrix
	compressed.feature.matrix <- cbind(reduced.features,per.patient.data.type.corr);

	# return the re-scaling and re-weighting parameter along with the compressed.feature.matrix for the matrix can be reproduced
	return(list(
		compressed.feature.matrix = compressed.feature.matrix
		reduced.features.rescale.recenter = reduced.features.rescale.recenter,
		pairwise.similiarity.features.rescale.recenter = pairwise.similiarity.features.rescale.recenter,
		reduced.features.rescale.denominator = reduced.features.rescale.denominator,
		pairwise.similiarity.features.rescale.denominator = pairwise.similiarity.features.rescale.denominator,
		reduced.features.weights = reduced.features.weights,
		pairwise.similiarity.features.weights = pairwise.similiarity.features.weights
		));
	}
