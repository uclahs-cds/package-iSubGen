subtype.by.h.clustering <- function(
	aberration.matrix,
	distance.metric,
	parent.output.dir,
	new.result.dir,
	subtype.table.file=NULL,
	max.num.subtypes=12,
	clustering.reps=1000,
	pFeature=0.8,
	pItem=0.8,
	verbose=FALSE,
	consensus.cluster.write.table=TRUE
	) {

	# move to the directory that we want the plots output to
	curr.dir <- getwd();
	setwd(parent.output.dir);

	if(nrow(aberration.matrix) > 1) {
		dianaHook <- function(this_dist,k) {
			num.patients <- nrow(as.matrix(this_dist));
			assignment <- rep(NA,num.patients);
			missing.patient <- apply(is.na(as.matrix(this_dist)),1,sum) == (num.patients - 1);
			missing.idx <- (num.patients - 1);
			while(any(is.na(as.matrix(this_dist)[!missing.patient,!missing.patient]))) {
				missing.patient <- apply(is.na(as.matrix(this_dist)),1,sum) >= missing.idx;
				missing.idx <- missing.idx - 1;
				}
			if(sum(!missing.patient) >= 2) {
				tmp <- diana(as.dist(as.matrix(this_dist)[!missing.patient,!missing.patient]),metric=distance.metric);
				if(sum(!missing.patient) >= k) {
					assignment[!missing.patient] <- cutree(tmp,k);
					}
				}
			return(assignment);
			}

		# run ConsensusCluster Plus and have it output the analysis plots to a pdf
		results <- ConsensusClusterPlus(
			d=t(aberration.matrix),
			plot='pdf',
			maxK=max.num.subtypes,
			distance=distance.metric,
			verbose=verbose,
			writeTable=consensus.cluster.write.table,
			title=new.result.dir,
			seed=17,
			finalLinkage='ward.D',
			innerLinkage='ward.D',
			clusterAlg=ifelse(any(is.na(aberration.matrix)),'dianaHook','hc'),
			pFeature=pFeature,
			pItem=pItem,
			reps=clustering.reps
			);

		# create a table of the subtypes determined for each number of clusters
		subtype.list <- list();
		for(i in 2:max.num.subtypes) {
			subtype.list[[paste('num_subtypes_',i,sep='')]] <- results[[i]]$consensusClass;
			}
		subtype.table <- as.data.frame(subtype.list);

	} else if(nrow(aberration.matrix) == 1){
		results <- diana(t(aberration.matrix));

		# create a table of the subtypes determined for each number of clusters
		subtype.list <- list();
		for(i in 2:max.num.subtypes) {
			subtype.list[[paste('num_subtypes_',i,sep='')]] <- cutree(results,i);
			}
		subtype.table <- as.data.frame(subtype.list);
		rownames(subtype.table) <- colnames(aberration.matrix);
		}

	setwd(curr.dir);

	if(!is.null(subtype.table.file)) {
		# save subtype table to file
		write.table(
			subtype.table,
			file=subtype.table.file,
			sep='\t',
			row.names=TRUE,
			col.names=TRUE,
			quote=FALSE
			);
		}

	return(list(clustering_result=results,subtype_table=subtype.table));
	}
