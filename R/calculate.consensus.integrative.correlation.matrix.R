calculate.consensus.integrative.correlation.matrix <- function(
	aberration.types,
	aberration.matrices,
	dist.metrics,
	correlation.method = 'spearman',
	filter.to.common.patients = FALSE,
	patients.to.return = NULL,
	patients.for.correlations = NULL,
	patient.proportion = 0.8,
	feature.proportion = 0.8,
	num.iterations=10,
	print.correlation.matrices.to.file=TRUE,
	print.dir='.'
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
	if(is.null(patients.for.correlations)) {
		patients.for.correlations <- patients;
	} else {
		patients.for.correlations <- intersect(patients.for.correlations,patients);
	}

	per.patient.aberration.type.corr <- list(); 
	for(i in 1:num.iterations) {
		set.seed(i);
		selected.patients <- sample(patients.for.correlations,round(length(patients.for.correlations)*patient.proportion));
		per.patient.aberration.type.corr[[i]] <- calculate.integrative.correlation.matrix(
			aberration.types=aberration.types,
			aberration.matrices=aberration.matrices,
			dist.metrics=dist.metrics,
			correlation.method = correlation.method,
			filter.to.common.patients = filter.to.common.patients,
			patients.to.return = patients.to.return,
			patients.for.correlations = selected.patients
		);
		if(print.correlation.matrices.to.file) {
			write.table(
				per.patient.aberration.type.corr[[i]],
				file=paste0(print.dir,'/',Sys.Date(),'_correlation_matrix_seed_',i,'.txt'),
				col.names=TRUE,
				row.names=TRUE,
				sep='\t',
				quote=FALSE
			);
		}
	}
	median.per.patient.aberration.type.corr <- matrix(NA, nrow=nrow(per.patient.aberration.type.corr[[1]]), ncol=ncol(per.patient.aberration.type.corr[[1]]));
	for(i in 1:nrow(median.per.patient.aberration.type.corr)) {
		for(j in 1:ncol(median.per.patient.aberration.type.corr)) {
			median.per.patient.aberration.type.corr[i,j] <- median(sapply(per.patient.aberration.type.corr, function(x) {per.patient.aberration.type.corr[[x]][i,j]}));
		}
	}

	return(median.per.patient.aberration.type.corr);
}
