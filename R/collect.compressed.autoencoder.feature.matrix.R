collect.compressed.autoencoder.feature.matrix <- function(
	data.types,
	data.matrices,
	autoencoders
	number.of.layers.in.decode=2,
	filter.to.common.patients=FALSE,
	patients.to.return=NULL
	) { 

	# pull out the patients to use
	patients <- NULL;
	for(data.type in data.types) {
		if(filter.to.common.patients) {
			if(is.null(patients)) {
				patients <- colnames(data.matrices[[data.type]])[grep('\\d',colnames(data.matrices[[data.type]]))];
			} else {
				patients <- intersect(patients,colnames(data.matrices[[data.type]])[grep('\\d',colnames(data.matrices[[data.type]]))]);
			}
		} else {
			# assume patient ids have at least one number in them and annotation columns don't
			patients <- union(patients,colnames(data.matrices[[data.type]]));
		}
	}
	patients <- sort(patients);
	if(is.null(patients.to.return)) {
		patients.to.return <- patients;
	} else {
		patients.to.return <- intersect(patients.to.return,patients);
	}

	autoencoder.results <- NULL;
	for(data.type in data.types) {
		if(data.type %in% names(autoencoders)) {
			model <- autoencoders[[data.type]];
			if(is.character(autoencoders[[data.type]]) && grep('hdf5$',autoencoders[[data.type]]) == 1) {
				model <- load_model_hdf5(
					autoencoders[[data.type]], 
					compile = FALSE);
			}

			autoencoder.data <- data.matrices[[data.type]];

			if(!is.null(data.quantiles.files[[data.type]])) {
				data.quantiles <- read.table(
					data.quantiles.files[[data.type]],
					header=TRUE,
					sep='\t');
				autoencoder.data <- sapply(1:nrow(autoencoder.data),function(x) {as.numeric(autoencoder.data[x,] > data.quantiles[x,2]) - as.numeric(autoencoder.data[x,] < data.quantiles[x,1])});
				autoencoder.data <- t(as.matrix(autoencoder.data));
				rownames(autoencoder.data) <- rownames(data.matrices[[data.type]]);
				colnames(autoencoder.data) <- colnames(data.matrices[[data.type]]);
			}

			autoencoder.data <- t(autoencoder.data);

			for(i in 1:number.of.layers.in.decode) {
				pop_layer(model);
			}

			reduced.features.autoencoder.predictions <- predict(model, x = autoencoder.data);
			rownames(reduced.features.autoencoder.predictions) <- rownames(autoencoder.data);
			colnames(reduced.features.autoencoder.predictions) <- paste0(data.type,1:ncol(reduced.features.autoencoder.predictions));

			autoencoder.results <- cbind(autoencoder.results,reduced.features.autoencoder.predictions);
		}
	}

	return(autoencoder.results);
	}
