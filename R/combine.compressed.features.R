combine.compressed.features <- function(
	reduced.features,
	pairwise.similarity.features,
	reduced.features.rescale.recenter = NA,
	pairwise.similarity.features.rescale.recenter = NA,
	reduced.features.rescale.denominator = NA,
	pairwise.similarity.features.rescale.denominator = NA,
	reduced.features.weights = rep(1,ncol(reduced.features)),
	pairwise.similarity.features.weights = rep(1,ncol(pairwise.similarity.features))
	) { 

	# check that the matrices have the same rows
	if(nrow(reduced.features) != nrow(pairwise.similarity.features)) {
		stop('reduced.features and pairwise.similarity.features matrices do not have the same number of patients');
	}
	if(all(rownames(reduced.features) != rownames(pairwise.similarity.features))) {
		stop('reduced.features and pairwise.similarity.features matrices rownames do not match');
	}
	# check for the right length and adjust by each parameter
	if(length(reduced.features.rescale.recenter) != 1 || !is.na(reduced.features.rescale.recenter)) {
		if(length(reduced.features.rescale.recenter) == 1 && reduced.features.rescale.recenter == 'mean') {
			reduced.features.rescale.recenter <- apply(reduced.features,2,mean);
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
	if(length(pairwise.similarity.features.rescale.recenter) != 1 || !is.na(pairwise.similarity.features.rescale.recenter)) {
		if(length(pairwise.similarity.features.rescale.recenter) == 1 && pairwise.similarity.features.rescale.recenter == 'mean') {
			pairwise.similarity.features.rescale.recenter <- apply(pairwise.similarity.features,2,mean);
		}
		if(length(pairwise.similarity.features.rescale.recenter) == 1) {
			pairwise.similarity.features <- pairwise.similarity.features - pairwise.similarity.features.rescale.recenter;
		} else if(length(pairwise.similarity.features.rescale.recenter) == ncol(pairwise.similarity.features)) {
			for(j in ncol(pairwise.similarity.features)) {
				pairwise.similarity.features[,j] <- pairwise.similarity.features[,j] - pairwise.similarity.features.rescale.recenter[j];
			}
		} else {
			warning('pairwise.similarity.features.rescale.recenter is not a valid length so no recentering was completed for pairwise.similarity.features');
		}
	}
	if(length(reduced.features.rescale.denominator) != 1 || !is.na(reduced.features.rescale.denominator)) {
		if(length(reduced.features.rescale.denominator) == 1 && reduced.features.rescale.denominator == 'sd') {
			reduced.features.rescale.denominator <- apply(reduced.features,2,sd);
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
	if(length(pairwise.similarity.features.rescale.denominator) != 1 ||  !is.na(pairwise.similarity.features.rescale.denominator)) {
		if(length(pairwise.similarity.features.rescale.denominator) == 1 && pairwise.similarity.features.rescale.denominator == 'sd') {
			pairwise.similarity.features.rescale.denominator <- apply(pairwise.similarity.features,2,sd);
		}
		if(length(pairwise.similarity.features.rescale.denominator) == 1) {
			pairwise.similarity.features <- pairwise.similarity.features/pairwise.similarity.features.rescale.denominator;
		} else if(length(pairwise.similarity.features.rescale.denominator) == ncol(pairwise.similarity.features)) {
			for(j in ncol(pairwise.similarity.features)) {
				pairwise.similarity.features[,j] <- pairwise.similarity.features[,j]/pairwise.similarity.features.rescale.denominator[j];
			}
		} else {
			warning('pairwise.similarity.features.rescale.denominator is not a valid length so no rescaling was completed for pairwise.similarity.features');
		}
	}
	if(length(reduced.features.weights) != 1 || !is.na(reduced.features.weights)) {
		if(length(reduced.features.weights) == 1) {
			reduced.features.weights <- rep(reduced.features.weights,ncol(reduced.features));
		}
		if(length(reduced.features.weights) == ncol(reduced.features)) {
			for(j in ncol(reduced.features)) {
				reduced.features[,j] <- reduced.features[,j]*reduced.features.weights[j];
			}
		} else {
			warning('reduced.features.weights is not a valid length so a weight of 1 was used compared to pairwise.similarity.features.weights');
		}
	}
	if(length(pairwise.similarity.features.weights) != 1 || !is.na(pairwise.similarity.features.weights)) {
		if(length(pairwise.similarity.features.weights) == 1) {
			pairwise.similarity.features.weights <- rep(pairwise.similarity.features.weights,ncol(pairwise.similarity.features));
		}
		if(length(pairwise.similarity.features.weights) == ncol(pairwise.similarity.features)) {
			for(j in ncol(pairwise.similarity.features)) {
				pairwise.similarity.features[,j] <- pairwise.similarity.features[,j]*pairwise.similarity.features.weights[j];
			}
		} else {
			warning('pairwise.similarity.features.weights is not a valid length so a weight of 1 was used compared to reduced.features.weights');
		}
	}

	# combine the matrices into one matrix
	compressed.feature.matrix <- cbind(reduced.features,pairwise.similarity.features);

	# return the re-scaling and re-weighting parameter along with the compressed.feature.matrix for the matrix can be reproduced
	return(list(
		compressed.feature.matrix = compressed.feature.matrix,
		reduced.features.rescale.recenter = reduced.features.rescale.recenter,
		pairwise.similarity.features.rescale.recenter = pairwise.similarity.features.rescale.recenter,
		reduced.features.rescale.denominator = reduced.features.rescale.denominator,
		pairwise.similarity.features.rescale.denominator = pairwise.similarity.features.rescale.denominator,
		reduced.features.weights = reduced.features.weights,
		pairwise.similarity.features.weights = pairwise.similarity.features.weights
		));
}
