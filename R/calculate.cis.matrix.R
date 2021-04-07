calculate.cis.matrix <- function(
	data.types,
	data.matrices,
	dist.metrics,
	correlation.method = 'spearman',
	filter.to.common.patients = FALSE,
	patients.to.return = NULL,
	patients.for.correlations = NULL,
	patient.proportion = 0.8,
	feature.proportion = 1,
	num.iterations = 10,
	print.intermediary.similarity.matrices.to.file = TRUE,
	print.dir = '.',
	patient.proportion.seeds = seq(1,num.iterations),
	feature.proportion.seeds = seq(1,num.iterations)
	) {

	# pull out the patients to use
	patients <- NULL;
	for(data.type in data.types) {
		if (filter.to.common.patients) {
			if (is.null(patients)) {
				patients <- colnames(data.matrices[[data.type]])[grep('\\d', colnames(data.matrices[[data.type]]))];
				}
			else {
				patients <- intersect(patients, colnames(data.matrices[[data.type]])[grep('\\d', colnames(data.matrices[[data.type]]))]);
				}
			}
		else {
			# assume patient ids have at least one number in them and annotation columns don't
			patients <- union(patients, colnames(data.matrices[[data.type]]));
			}
		}
	if (is.null(patients.for.correlations)) {
		patients.for.correlations <- patients;
		}
	else {
		patients.for.correlations <- intersect(patients.for.correlations, patients);
		}

	# repeatly subsample the dataset and calculate integrative similarity
	per.patient.data.type.corr <- list();
	for(i in 1:num.iterations) {
		set.seed(patient.proportion.seeds[i]);
		selected.patients <- sample(patients.for.correlations,round(length(patients.for.correlations)*patient.proportion));
		data.matrices.subset <- data.matrices;
		# if the feature proportion is 1, then we don't need to filter the features
		# if its not 1, then the features need to be selected for the iteration
		if (feature.proportion != 1) {
			for(data.type in data.types) {
				set.seed(feature.proportion.seeds[i]);
				selected.features <- sample(rownames(data.matrices[[data.type]]), ceiling(nrow(data.matrices[[data.type]])*feature.proportion));
				data.matrices.subset[[data.type]] <- data.matrices.subset[[data.type]][selected.features,];
				}
			}
		per.patient.data.type.corr[[i]] <- calculate.integrative.similarity.matrix(
			data.types = data.types,
			data.matrices = data.matrices.subset,
			dist.metrics = dist.metrics,
			correlation.method = correlation.method,
			filter.to.common.patients = filter.to.common.patients,
			patients.to.return = patients.to.return,
			patients.for.correlations = selected.patients
			);
		if (print.intermediary.similarity.matrices.to.file) {
			write.table(
				per.patient.data.type.corr[[i]],
				file = paste0(print.dir,'/',Sys.Date(),'_correlation_matrix_seed_',i,'.txt'),
				col.names = TRUE,
				row.names = TRUE,
				sep = '\t',
				quote = FALSE
				);
			}
		}

	# combine iterations to calculate CIS by taking the median similarity
	median.per.patient.data.type.corr <- matrix(
		data = NA,
		nrow = nrow(per.patient.data.type.corr[[1]]),
		ncol = ncol(per.patient.data.type.corr[[1]])
		);
	for(i in 1:nrow(median.per.patient.data.type.corr)) {
		for(j in 1:ncol(median.per.patient.data.type.corr)) {
			median.per.patient.data.type.corr[i,j] <- median(sapply(per.patient.data.type.corr, function(x) {x[i,j]}));
			}
		}
	rownames(median.per.patient.data.type.corr) <- rownames(per.patient.data.type.corr[[1]]);
	colnames(median.per.patient.data.type.corr) <- colnames(per.patient.data.type.corr[[1]]);

	return(median.per.patient.data.type.corr);
	}
