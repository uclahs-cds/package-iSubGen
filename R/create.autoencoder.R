create.autoencoder <- function(
	aberration.type,
	aberration.matrix,
	encoder.layers.node.nums = c(15,2),
	autoencoder.activation = 'tanh',
	optimization.loss.function = 'mean_squared_error',
	model.file.output.dir = '.'
	) { 

	# if the same activation function is going to be used for all layers expand it to a vector
	if(length(autoencoder.activation) == 1) {
		autoencoder.activation <- rep(autoencoder.activation,length(encoder.layers.node.nums));
	}

	# create an empty neural net and add the first/input layer
	model <- keras_model_sequential();
	model %>%
		layer_dense(units = encoder.layers.node.nums[1], activation = autoencoder.activation[1], input_shape = nrow(aberration.matrix));

	# add the rest of the neural net layers
	if(length(encoder.layers.node.nums) > 1) {
		# encoding layers
		for(i in 2:(length(encoder.layers.node.nums)-1)) {
			model %>%
				layer_dense(units = encoder.layers.node.nums[i], activation = autoencoder.activation[i]);
		}
		model %>% layer_dense(units = encoder.layers.node.nums[length(encoder.layers.node.nums)], activation = autoencoder.activation[length(encoder.layers.node.nums)],name='bottleneck');
		# decoding layers
		for(i in (length(encoder.layers.node.nums)-1):1) {
			model %>%
				layer_dense(units = encoder.layers.node.nums[i], activation = autoencoder.activation[i]);
		}
		
		# output layer
		model %>%
			layer_dense(units = nrow(aberration.matrix));
	}

	# set up training parameters
	model %>% compile(
		loss = optimization.loss.function, 
		optimizer = 'adam'
		);

	ae.output.file <- paste0(sub('/$','',model.file.output.dir),'/',aberration.type,'_model.hdf5');
	checkpoint <- callback_model_checkpoint(
		filepath = ae.output.file, 
		save_best_only = TRUE,
#		period = 1,
		verbose = 1
		);

	early.stopping <- callback_early_stopping(patience = 5);

	batch.size <- 50;
	if(batch.size > ncol(aberration.matrix)) {
		batch.size <- ncol(aberration.matrix) -1;
	}

	# train the neural net
	model %>% fit(
		x = t(aberration.matrix), 
		y = t(aberration.matrix), 
		epochs = 350, 
		batch_size = 50,
		validation_data = list(t(aberration.matrix), t(aberration.matrix)), 
		callbacks = list(checkpoint,early.stopping)
		);

	return(list(autoencoder=model,autoencoder.file=ae.output.file));
	}
