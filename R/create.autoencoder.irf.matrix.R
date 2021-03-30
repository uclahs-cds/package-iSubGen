create.autoencoder.irf.matrix <- function(
	data.types,
	data.matrices,
	autoencoders,
	filter.to.common.patients = FALSE,
	patients.to.return = NULL
	) {

	# pull out the patients to use
	patients <- NULL;
	for(data.type in data.types) {
		if (filter.to.common.patients) {
			# assume patient ids have at least one number in them and annotation columns don't
			if (is.null(patients)) {
				patients <- colnames(data.matrices[[data.type]])[grep('\\d', colnames(data.matrices[[data.type]]))];
				}
			else {
				patients <- intersect(patients, colnames(data.matrices[[data.type]])[grep('\\d', colnames(data.matrices[[data.type]]))]);
				}
			}
		else {
			patients <- union(patients, colnames(data.matrices[[data.type]]));
			}
		}
	patients <- sort(patients);
	if (is.null(patients.to.return)) {
		patients.to.return <- patients;
		}
	else {
		patients.to.return <- intersect(patients.to.return, patients);
		}

	# go through each data.type and get the corresponding compressed autoencoder features
	irf.matrix <- NULL;
	for(data.type in data.types) {
		if (data.type %in% names(autoencoders)) {
			# load the neural net for the data.type
			model <- autoencoders[[data.type]];
			if (is.character(autoencoders[[data.type]]) && grep('hdf5$',autoencoders[[data.type]]) == 1) {
				model <- load_model_hdf5(
					autoencoders[[data.type]], 
					compile = FALSE);
				}

			# create the autoencoder encoding layers from input to the bottleneck layer
			intermediate.layer.model <- keras_model(inputs = model$input, outputs = get_layer(model, "bottleneck")$output);

			# get the bottleneck values from the autoencoder
			bottleneck.values <- predict(intermediate.layer.model,x=t(data.matrices[[data.type]]));
			rownames(bottleneck.values) <- colnames(data.matrices[[data.type]]);
			colnames(bottleneck.values) <- paste0(data.type,1:ncol(bottleneck.values));

			# add the IRF features from this data type to the IRF matrix
			if (is.null(irf.matrix)) {
				irf.matrix <- bottleneck.values;
				}
			else {
				irf.matrix <- cbind(irf.matrix, bottleneck.values);
				}
			}
		}

	return(irf.matrix);
	}
