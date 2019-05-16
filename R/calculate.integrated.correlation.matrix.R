calculate.integrated.correlation.matrix <- function(
	data.types,
	data.matrices,
	dist.metrics,
	correlation.method='spearman',
	filter.to.common.patients=FALSE,
	patients.to.return=NULL,
	patients.for.correlations=NULL
	) { 

	# pull out the patients to use
	patients <- NULL;
	for(data.type in data.types) {
		if(filter.to.common.patients) {
			if(is.null(patients)) {
				patients <- colnames(data.matrices[[data.type]])[grep('\\d',colnames(data.matrices[[data.type]]))];
				}
			else {
				patients <- intersect(patients,colnames(data.matrices[[data.type]])[grep('\\d',colnames(data.matrices[[data.type]]))]);
				}
			}
		else {
			# assume patient ids have at least one number in them and annotation columns don't
			patients <- union(patients,colnames(data.matrices[[data.type]]));
			}
		}
	patients <- sort(patients);
	if(is.null(patients.to.return)) {
		patients.to.return <- patients;
		}
	else {
		patients.to.return <- intersect(patients.to.return,patients);
		}
	if(is.null(patients.for.correlations)) {
		patients.for.correlations <- patients;
		}
	else {
		patients.for.correlations <- intersect(patients.for.correlations,patients);
		}

	patient.pairs <- as.character(sapply(1:(length(patients.for.correlations)),function(x) {paste0(patients.for.correlations[x],':',patients.to.return)}));
	patients <- sort(unique(c(patients.to.return,patients.for.correlations)));

	# calculate pair-wise distances
	sample.paired.dists <- list();
	sample.paired.dists.matrix <- matrix(NA,ncol=length(data.types),nrow=length(patient.pairs));
	colnames(sample.paired.dists.matrix) <- data.types;
	rownames(sample.paired.dists.matrix) <- patient.pairs;
	for(data.type in data.types) {
		# filter to required patients 
		data.matrices[[data.type]] <- data.matrices[[data.type]][,sort(colnames(data.matrices[[data.type]])[colnames(data.matrices[[data.type]]) %in% patients])];

		# set up the matrix of distance pairs that are required
		dist.matrix <- matrix(NA,ncol=length(patients.for.correlations),nrow=length(patients.to.return));
		colnames(dist.matrix) <- patients.for.correlations;
		rownames(dist.matrix) <- patients.to.return;

		# determine the most efficient approach for calculating the distances
		dist.calc.operations <- list();
		pr <- length(patients.to.return);
		pc <- length(patients.for.correlations);
		opt.num.of.return.to.calc.at.once <- order(sapply(1:pr,function(k) {(pc+k)^2*(ceiling(pr/k)-1)+(pc+ifelse((pr%%k)== 0,k,pr%%k))^2}))[1];
		pr.tracker <- 0;
		while(pr.tracker < pr) {
			if((pr.tracker+opt.num.of.return.to.calc.at.once) < pr) {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):(pr.tracker+opt.num.of.return.to.calc.at.once)];
				pr.tracker <- pr.tracker+opt.num.of.return.to.calc.at.once;
				}
			else {
				dist.calc.operations[[as.character(pr.tracker)]] <- patients.to.return[(pr.tracker+1):pr];
				pr.tracker <- pr;
				}
			}

		# calculate distances and fill in matrix
		for(dist.op in 1:length(dist.calc.operations)) {
			if(class(dist.metrics[[data.type]]) == 'character') {
				if(dist.metrics[[data.type]] %in% c('pearson','spearman')) {
					dist.result <- as.dist(1 - cor(data.matrices[[data.type]][,unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)),drop=FALSE], use = 'pairwise', method = dist.metrics[[data.type]]));
					}
				else {
					dist.result <- BoutrosLab.dist.overload::dist(t(data.matrices[[data.type]][,unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)),drop=FALSE]),method=dist.metrics[[data.type]]);
					}
				}
			else if(class(dist.metrics[[data.type]]) == 'function') {
				dist.result <- as.dist((dist.metrics[[data.type]])(t(data.matrices[[data.type]][,unique(c(dist.calc.operations[dist.op][[1]],patients.for.correlations)),drop=FALSE])));
				}
			else {
				stop(paste0('invalid option for ',data.type,' distance metric: ',dist.metrics[[data.type]]));
				}
			dist.matrix[dist.calc.operations[dist.op][[1]],patients.for.correlations] <- as.matrix(dist.result)[dist.calc.operations[dist.op][[1]],patients.for.correlations];
			}


		sample.paired.dists[[data.type]] <- dist.matrix;
		patient.pairs <- as.character(sapply(1:(ncol(dist.matrix)),function(x) {paste0(colnames(dist.matrix)[x],':',rownames(dist.matrix))}));
		sample.paired.dists.matrix[patient.pairs,data.type] <- as.numeric(dist.matrix);
		}
	sample.paired.dists.matrix <- sample.paired.dists.matrix[apply(sample.paired.dists.matrix,1,function(x) {sum(!is.na(x))}) >=2,];

	split.rownames <- strsplit(rownames(sample.paired.dists.matrix),':');
	pair.patient1 <- sapply(1:length(split.rownames),function(i) {split.rownames[[i]][1]});
	pair.patient2 <- sapply(1:length(split.rownames),function(i) {split.rownames[[i]][2]});

	# calculate distance correlations between date types
	per.patient.data.type.corr <- matrix(NA, nrow=length(patients.to.return),ncol=length(data.types)*(length(data.types)-1)/2);
	rownames(per.patient.data.type.corr) <- patients.to.return;
	colnames(per.patient.data.type.corr) <- seq(1,ncol(per.patient.data.type.corr));
	data.type.pair.counter <- 0;
	for(i in 1:(length(data.types)-1)) {
		for(j in (i+1):length(data.types)) {
			data.type.pair.counter <- data.type.pair.counter +1;
			colnames(per.patient.data.type.corr)[data.type.pair.counter] <- paste0(data.types[i],':',data.types[j]);
			not.na.rows <- (!is.na(sample.paired.dists.matrix[,i])) & (!is.na(sample.paired.dists.matrix[,j]));
			for(patient in rownames(per.patient.data.type.corr)) {
				rows.to.use <- which(pair.patient2 == patient & not.na.rows);
				if(length(rows.to.use) > 1) {
					per.patient.data.type.corr[patient,data.type.pair.counter] <- cor(c(0,sample.paired.dists.matrix[rows.to.use,i]),c(0,sample.paired.dists.matrix[rows.to.use,j]),method=correlation.method);
					}
				}
			}
		}
	per.patient.data.type.corr <- per.patient.data.type.corr[which(apply(per.patient.data.type.corr,1,function(x) { sum(!is.na(x)) }) > 0),];

	return(per.patient.data.type.corr);
	}
