calculate.integrative.similarity.matrix <- function(
	data.types,
	data.matrices,
	dist.metrics,
	correlation.method = 'spearman',
	filter.to.common.patients = FALSE,
	patients.to.return = NULL,
	patients.for.correlations = NULL
	) { 

	# pull out the patients to use
	patients <- NULL;
	for(data.type in data.types) {
		if (filter.to.common.patients) {
			# assume patient IDs have at least one number in them and annotation columns don't
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

	# to calculate integrative similarity values the set of patients that you want values for
	# does not need to be the same set of patients used for comparison when calculating the values
	# patients.to.return are the patients that you want similarity values for
	# patients.for.correlations are the patients to use for calculating similarity values
	# first we need to find the overlap with the patients in the data matrices
	if (is.null(patients.to.return)) {
		patients.to.return <- patients;
		}
	else {
		patients.to.return <- intersect(patients.to.return, patients);
		} 
	if (is.null(patients.for.correlations)) {
		patients.for.correlations <- patients;
		}
	else {
		patients.for.correlations <- intersect(patients.for.correlations, patients);
		}
	patient.pairs <- as.character(sapply(1:(length(patients.for.correlations)), function(x) {paste0(patients.for.correlations[x], ':', patients.to.return)}));
	patients <- sort(unique(c(patients.to.return, patients.for.correlations)));

	# calculate pair-wise distances
	patient.paired.dists <- list();
	patient.paired.dists.matrix <- matrix(NA, ncol = length(data.types), nrow = length(patient.pairs));
	colnames(patient.paired.dists.matrix) <- data.types;
	rownames(patient.paired.dists.matrix) <- patient.pairs;
	for(data.type in data.types) {
		# filter to required patients
		data.matrices[[data.type]] <- data.matrices[[data.type]][,sort(colnames(data.matrices[[data.type]])[colnames(data.matrices[[data.type]]) %in% patients])];

		# set up the matrix of distance pairs that are required
		dist.matrix <- matrix(NA, ncol = length(patients.for.correlations), nrow = length(patients.to.return));
		colnames(dist.matrix) <- patients.for.correlations;
		rownames(dist.matrix) <- patients.to.return;

		# determine the most efficient approach for calculating the distances
		# we need pr.num x pc.num distances but calculating distances creates 
		# matrix with the same columns and rows which would be (pr.num + pc.num)^2
		# calculating distances in sets can mean less unnecessary calculations are done
		# comparing the each of the patients.to.return to the patients.for.correlation
		# select the number of patients.to.return to compare to patients.for.correlation at a time
		# (pc.num + k) is the patients per comparison
		# then there will be (ceiling(pr.num/k)-1) comparisons
		# and then an additional (pc.num + j) where j is the remainder not calculated
		dist.calc.operations <- list();
		pr.num <- length(patients.to.return);
		pc.num <- length(patients.for.correlations);
		opt.num.of.return.to.calc.at.once <- order(sapply(
			1:pr.num,
			function(k) {(pc.num + k)^2 * (ceiling(pr.num/k)-1) + (pc.num + ifelse((pr.num%%k) == 0, k, pr.num%%k))^2}
			))[1];
		pr.tracker <- 0;
		while(pr.tracker < pr.num) {
			if ((pr.tracker + opt.num.of.return.to.calc.at.once) < pr.num) {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):(pr.tracker+opt.num.of.return.to.calc.at.once)];
				pr.tracker <- pr.tracker + opt.num.of.return.to.calc.at.once;
				}
			else {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):pr.num];
				pr.tracker <- pr.num;
				}
			}

		# calculate distances and fill in patient by patient distance matrix
		for(dist.op in 1:length(dist.calc.operations)) {
			if (class(dist.metrics[[data.type]]) == 'character') {

				if (dist.metrics[[data.type]] %in% c('pearson','spearman')) {
					# if the distance metric is a correlation, convert the correlation into a distance
					dist.result <- as.dist(
						1 - cor(data.matrices[[data.type]][,intersect(colnames(data.matrices[[data.type]]),unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)))],
						use = 'pairwise',
						method = dist.metrics[[data.type]])
						);
					}
				else {
					# for distances other than correlations use the distance function
					dist.result <- distance(
						t(data.matrices[[data.type]][,intersect(colnames(data.matrices[[data.type]]),unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)))]),
						method = dist.metrics[[data.type]],
						use.row.names = TRUE
						);
					}
				}
			else if (class(dist.metrics[[data.type]]) == 'function') {
				dist.result <- as.dist((dist.metrics[[data.type]])(t(data.matrices[[data.type]][,intersect(colnames(dist.metrics[[data.type]]), unique(c(dist.calc.operations[dist.op][[1]], patients.for.correlations)))])));
				}
			else {
				stop(paste0('invalid option for ', data.type, ' distance metric: ', dist.metrics[[data.type]]));
				}
			dist.matrix[
				intersect(colnames(data.matrices[[data.type]]), dist.calc.operations[dist.op][[1]]),
				intersect(colnames(data.matrices[[data.type]]), patients.for.correlations)
				] <- as.matrix(dist.result)[
					intersect(colnames(data.matrices[[data.type]]), dist.calc.operations[dist.op][[1]]),
					intersect(colnames(data.matrices[[data.type]]), patients.for.correlations)
					];
			}

		patient.paired.dists[[data.type]] <- dist.matrix;
		patient.pairs <- as.character(sapply(1:(ncol(dist.matrix)),function(x) {paste0(colnames(dist.matrix)[x], ':', rownames(dist.matrix))}));
		patient.paired.dists.matrix[patient.pairs,data.type] <- as.numeric(dist.matrix);
		}
	# find the rows with at least one value (not na) beyond the patient by patient comparison
	patient.paired.dists.matrix <- patient.paired.dists.matrix[apply(patient.paired.dists.matrix, 1, function(x) {sum(!is.na(x))}) > 1,];

	split.rownames <- strsplit(rownames(patient.paired.dists.matrix), ':');
	pair.patient1 <- sapply(1:length(split.rownames), function(i) {split.rownames[[i]][1]});
	pair.patient2 <- sapply(1:length(split.rownames), function(i) {split.rownames[[i]][2]});

	# calculate correlations (or integrative similarity) between data types
	per.patient.data.type.corr <- matrix(NA, nrow = length(patients.to.return), ncol = length(data.types)*(length(data.types)-1)/2);
	rownames(per.patient.data.type.corr) <- patients.to.return;
	colnames(per.patient.data.type.corr) <- seq(1, ncol(per.patient.data.type.corr));
	data.type.pair.counter <- 0;
	for(i in 1:(length(data.types)-1)) {
		for(j in (i+1):length(data.types)) {

			data.type.pair.counter <- data.type.pair.counter + 1;
			colnames(per.patient.data.type.corr)[data.type.pair.counter] <- paste0(data.types[i], ':', data.types[j]);
			not.na.rows <- (!is.na(patient.paired.dists.matrix[,i])) & (!is.na(patient.paired.dists.matrix[,j]));

			for(patient in rownames(per.patient.data.type.corr)) {
				rows.to.use <- which(pair.patient2 == patient & not.na.rows);
				if (length(rows.to.use) > 1) {
					per.patient.data.type.corr[patient,data.type.pair.counter] <- cor(
						c(0, patient.paired.dists.matrix[rows.to.use,i]),
						c(0, patient.paired.dists.matrix[rows.to.use,j]),
						method = correlation.method
						);
					}
				}
			}
		}
	per.patient.data.type.corr <- per.patient.data.type.corr[which(apply(per.patient.data.type.corr, 1, function(x) { sum(!is.na(x)) }) > 0), , drop=FALSE];

	return(per.patient.data.type.corr);
	}
