calculate.integrative.correlation.matrix <- function(
	aberration.types,
	aberration.matrices,
	dist.metrics,
	correlation.method = 'spearman',
	filter.to.common.patients = FALSE,
	patients.to.return = NULL,
	patients.for.correlations = NULL
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
	if(is.null(patients.for.correlations)) {
		patients.for.correlations <- patients;
	} else {
		patients.for.correlations <- intersect(patients.for.correlations,patients);
	}
	patient.pairs <- as.character(sapply(1:(length(patients.for.correlations)),function(x) {paste0(patients.for.correlations[x],':',patients.to.return)}));
	patients <- sort(unique(c(patients.to.return,patients.for.correlations)));

	# calculate pair-wise distances
	patient.paired.dists <- list();
	patient.paired.dists.matrix <- matrix(NA, ncol = length(aberration.types), nrow = length(patient.pairs));
	colnames(patient.paired.dists.matrix) <- aberration.types;
	rownames(patient.paired.dists.matrix) <- patient.pairs;
	for(aberration.type in aberration.types) {
		# filter to required patients 
		aberration.matrices[[aberration.type]] <- aberration.matrices[[aberration.type]][,sort(colnames(aberration.matrices[[aberration.type]])[colnames(aberration.matrices[[aberration.type]]) %in% patients])];

		# set up the matrix of distance pairs that are required
		dist.matrix <- matrix(NA, ncol = length(patients.for.correlations), nrow = length(patients.to.return));
		colnames(dist.matrix) <- patients.for.correlations;
		rownames(dist.matrix) <- patients.to.return;

		# determine the most efficient approach for calculating the distances
		dist.calc.operations <- list();
		pr <- length(patients.to.return);
		pc <- length(patients.for.correlations);
		opt.num.of.return.to.calc.at.once <- order(sapply(1:pr, function(k) {(pc + k)^2 * (ceiling(pr/k)-1) + (pc + ifelse((pr%%k) == 0, k, pr%%k))^2}))[1];
		pr.tracker <- 0;
		while(pr.tracker < pr) {
			if((pr.tracker + opt.num.of.return.to.calc.at.once) < pr) {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):(pr.tracker+opt.num.of.return.to.calc.at.once)];
				pr.tracker <- pr.tracker + opt.num.of.return.to.calc.at.once;
			} else {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):pr];
				pr.tracker <- pr;
			}
		}

		# calculate distances and fill in matrix
		for(dist.op in 1:length(dist.calc.operations)) {
			if(class(dist.metrics[[aberration.type]]) == 'character') {
				if(dist.metrics[[aberration.type]] %in% c('pearson','spearman')) {
					dist.result <- as.dist(
						1 - cor(aberration.matrices[[aberration.type]][,intersect(colnames(aberration.matrices[[aberration.type]]),unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)))],
						use = 'pairwise',
						method = dist.metrics[[aberration.type]])
						);
				}
				else {
					dist.result <- distance(
						t(aberration.matrices[[aberration.type]][,intersect(colnames(aberration.matrices[[aberration.type]]),unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)))]),
						method = dist.metrics[[aberration.type]],
						use.row.names=TRUE
						);
		#			colnames(dist.result) <- unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations));
		#			rownames(dist.result) <- unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations));
					
				}
			} else if(class(dist.metrics[[aberration.type]]) == 'function') {
				dist.result <- as.dist((dist.metrics[[aberration.type]])(t(aberration.matrices[[aberration.type]][,intersect(colnames(dist.metrics[[aberration.type]]),unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)))])));
			} else {
				stop(paste0('invalid option for ',aberration.type,' distance metric: ',dist.metrics[[aberration.type]]));
			}
			dist.matrix[
				intersect(colnames(aberration.matrices[[aberration.type]]),dist.calc.operations[dist.op][[1]]),
				intersect(colnames(aberration.matrices[[aberration.type]]),patients.for.correlations)
				] <- as.matrix(dist.result)[
					intersect(colnames(aberration.matrices[[aberration.type]]),dist.calc.operations[dist.op][[1]]),
					intersect(colnames(aberration.matrices[[aberration.type]]),patients.for.correlations)
					];
		}

		patient.paired.dists[[aberration.type]] <- dist.matrix;
		patient.pairs <- as.character(sapply(1:(ncol(dist.matrix)),function(x) {paste0(colnames(dist.matrix)[x],':',rownames(dist.matrix))}));
		patient.paired.dists.matrix[patient.pairs,aberration.type] <- as.numeric(dist.matrix);
	}
	patient.paired.dists.matrix <- patient.paired.dists.matrix[apply(patient.paired.dists.matrix,1,function(x) {sum(!is.na(x))}) >= 2,];

	split.rownames <- strsplit(rownames(patient.paired.dists.matrix),':');
	pair.patient1 <- sapply(1:length(split.rownames),function(i) {split.rownames[[i]][1]});
	pair.patient2 <- sapply(1:length(split.rownames),function(i) {split.rownames[[i]][2]});

	# calculate distance correlations between date types
	per.patient.aberration.type.corr <- matrix(NA, nrow = length(patients.to.return), ncol = length(aberration.types)*(length(aberration.types)-1)/2);
	rownames(per.patient.aberration.type.corr) <- patients.to.return;
	colnames(per.patient.aberration.type.corr) <- seq(1,ncol(per.patient.aberration.type.corr));
	aberration.type.pair.counter <- 0;
	for(i in 1:(length(aberration.types)-1)) {
		for(j in (i+1):length(aberration.types)) {

			aberration.type.pair.counter <- aberration.type.pair.counter +1;
			colnames(per.patient.aberration.type.corr)[aberration.type.pair.counter] <- paste0(aberration.types[i],':',aberration.types[j]);
			not.na.rows <- (!is.na(patient.paired.dists.matrix[,i])) & (!is.na(patient.paired.dists.matrix[,j]));

			for(patient in rownames(per.patient.aberration.type.corr)) {
				rows.to.use <- which(pair.patient2 == patient & not.na.rows);
				if(length(rows.to.use) > 1) {
					per.patient.aberration.type.corr[patient,aberration.type.pair.counter] <- cor(
						c(0,patient.paired.dists.matrix[rows.to.use,i]),
						c(0,patient.paired.dists.matrix[rows.to.use,j]),
						method = correlation.method
						);
				}
			}
		}
	}
	per.patient.aberration.type.corr <- per.patient.aberration.type.corr[which(apply(per.patient.aberration.type.corr,1,function(x) { sum(!is.na(x)) }) > 0),,drop=FALSE];

	return(per.patient.aberration.type.corr);
}
