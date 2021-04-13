cluster.patients <- function(
	data.matrix,
	distance.metric,
	parent.output.dir,
	new.result.dir,
	subtype.table.file = NULL,
	max.num.subtypes = 12,
	clustering.reps = 1000,
	proportion.features = 0.8,
	proportion.patients = 0.8,
	verbose = FALSE,
	consensus.cluster.write.table = TRUE
	) {

	# move to the directory that we want the plots output to
	func.start.dir <- getwd();
	setwd(parent.output.dir);
	on.exit(setwd(func.start.dir));

	if(nrow(data.matrix) > 1) {

		# Defining function for when there is missing data here because
		# using variables from this function (ex. num.patients, distance.metric)
		# that are not passed to the function and don't have control of changing
		# function call within ConsensusCluster to pass more arguments to 
		# dianaWithMissingPatients. Creating the function here means the variables
		# be updated before being used.
		dianaWithMissingPatients <- function(dist,k) {
			num.patients <- nrow(as.matrix(dist));
			assignment <- rep(NA, num.patients);
			missing.patient <- apply(is.na(as.matrix(dist)),1,sum) == (num.patients - 1);
			missing.idx <- (num.patients - 1);
			while(any(is.na(as.matrix(dist)[!missing.patient,!missing.patient]))) {
				missing.patient <- apply(is.na(as.matrix(dist)), 1, sum) >= missing.idx;
				missing.idx <- missing.idx - 1;
				}
			if (sum(!missing.patient) >= 2) {
				clusters <- diana(as.dist(as.matrix(dist)[!missing.patient,!missing.patient]), metric = distance.metric);
				if (sum(!missing.patient) >= k) {
					assignment[!missing.patient] <- cutree(clusters, k);
					}
				}
			return(assignment);
			}

		# run ConsensusCluster Plus and have it output the analysis plots to a pdf
		cluster.result <- ConsensusClusterPlus(
			d = t(data.matrix),
			plot = 'pdf',
			maxK = max.num.subtypes,
			distance = distance.metric,
			verbose = verbose,
			writeTable = consensus.cluster.write.table,
			title = new.result.dir,
			seed = 17,
			finalLinkage = 'ward.D',
			innerLinkage = 'ward.D',
			clusterAlg = ifelse(any(is.na(data.matrix)),'dianaWithMissingPatients','hc'),
			pFeature = proportion.features,
			pItem = proportion.patients,
			reps = clustering.reps
			);

		# create a table of the subtypes determined for each number of clusters
		subtype.list <- list();
		for(i in 2:max.num.subtypes) {
			subtype.list[[paste('num_subtypes_', i ,sep='')]] <- cluster.result[[i]]$consensusClass;
			}
		subtype.table <- as.data.frame(subtype.list);

		}
	else if (nrow(data.matrix) == 1){
		cluster.result <- diana(t(data.matrix));

		# create a table of the subtypes determined for each number of clusters
		subtype.list <- list();
		for(i in 2:max.num.subtypes) {
			subtype.list[[paste('num_subtypes_', i, sep='')]] <- cutree(cluster.result, i);
			}
		subtype.table <- as.data.frame(subtype.list);
		rownames(subtype.table) <- colnames(data.matrix);
		}

	setwd(func.start.dir);

	if (!is.null(subtype.table.file)) {
		# save subtype table to file
		write.table(
			subtype.table,
			file = subtype.table.file,
			sep = '\t',
			row.names = TRUE,
			col.names = TRUE,
			quote = FALSE
			);
		}

	return(list(clustering_result = cluster.result, subtype_table = subtype.table));
	}
