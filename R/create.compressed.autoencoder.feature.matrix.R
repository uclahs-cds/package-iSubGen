create.compressed.autoencoder.feature.matrix <- function(
	aberration.types,
	aberration.matrices,
	autoencoders,
	number.of.layers.in.decode=2,
	filter.to.common.patients=FALSE,
	patients.to.return=NULL
	) { 

	# pull out the patients to use
	patients <- NULL;
	for(aberration.type in aberration.types) {
		if(filter.to.common.patients) {
			if(is.null(patients)) {
				patients <- colnames(aberration.matrices[[aberration.type]])[grep('\\d',colnames(aberration.matrices[[aberration.type]]))];
			} else {
				patients <- intersect(patients,colnames(aberration.matrices[[aberration.type]])[grep('\\d',colnames(aberration.matrices[[aberration.type]]))]);
			}
		} else {
			# assume patient ids have at least one number in them and annotation columns don't
			patients <- union(patients,colnames(aberration.matrices[[aberration.type]]));
		}
	}
	patients <- sort(patients);
	if(is.null(patients.to.return)) {
		patients.to.return <- patients;
	} else {
		patients.to.return <- intersect(patients.to.return,patients);
	}

	# go through each aberration type and get the corresponding compressed autoencoder features
	autoencoder.results <- NULL;
	for(aberration.type in aberration.types) {
		if(aberration.type %in% names(autoencoders)) {
			# load the neural net for the aberration type
			model <- autoencoders[[aberration.type]];
			if(is.character(autoencoders[[aberration.type]]) && grep('hdf5$',autoencoders[[aberration.type]]) == 1) {
				model <- load_model_hdf5(
					autoencoders[[aberration.type]], 
					compile = FALSE);
			}

			# remove the decoding layers
			for(i in 1:number.of.layers.in.decode) {
				pop_layer(model);
			}

			autoencoder.aberration <- t(aberration.matrices[[aberration.type]]);
			reduced.features.autoencoder.predictions <- predict(model, x = autoencoder.aberration);
			rownames(reduced.features.autoencoder.predictions) <- rownames(autoencoder.aberration);
			colnames(reduced.features.autoencoder.predictions) <- paste0(aberration.type,1:ncol(reduced.features.autoencoder.predictions));

			autoencoder.results <- cbind(autoencoder.results,reduced.features.autoencoder.predictions);
		}
	}

	return(autoencoder.results);
}
