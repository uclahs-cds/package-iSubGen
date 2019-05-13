create.autoencoder <- function(
	data.type,
	data.matrix,
	decoding.layer.node.num=c(15,2),
	autoencoder.activation='tanh',
	optimization.loss.function='mean_squared_error',
	model.file.output.dir='.'
	) { 

	model <- keras_model_sequential();
	model %>%
		layer_dense(units = decoding.layer.node.num[1], activation = autoencoder.activation, input_shape = nrow(data.matrix));

	if(length(decoding.layer.node.num) > 1) {
		for(i in 2:length(decoding.layer.node.num)) {
			model %>%
				layer_dense(units = decoding.layer.node.num[i], activation = autoencoder.activation);
			}
		for(i in (length(decoding.layer.node.num)-1):1) {
			model %>%
				layer_dense(units = decoding.layer.node.num[i], activation = autoencoder.activation);
			}
								
		model %>%
			layer_dense(units = nrow(data.matrix));
		}

	model %>% compile(
		loss = optimization.loss.function, 
		optimizer = 'adam'
		);

	ae.output.file <- paste0(sub('/$','',model.file.output.dir),'/',data.type,'_model.hdf5');
	checkpoint <- callback_model_checkpoint(
		filepath = ae.output.file, 
		save_best_only = TRUE, 
		period = 1,
		verbose = 1
		);

	early.stopping <- callback_early_stopping(patience = 5);

	batch.size <- 50;
	if(batch.size > ncol(data.matrix)) {
		batch.size <- ncol(data.matrix) -1;
		}

	model %>% fit(
		x = t(data.matrix), 
		y = t(data.matrix), 
		epochs = 350, 
		batch_size = 50,
		validation_data = list(t(data.matrix), t(data.matrix)), 
		callbacks = list(checkpoint,early.stopping)
		);

	return(list(autoencoder=model,autoencoder.file=ae.output.file));
	}
