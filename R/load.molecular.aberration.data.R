# returns the molecular aberration profiles/feature annotation
load.molecular.aberration.data <- function(file, patients=NULL, annotation.fields=NULL
	) {

	# read in matrix with annotation and mutation data
	mut.data.and.anno <- read.table(
		file,
		sep='\t',
		header=TRUE
		);

	### pull out the molecular annotation profiles ###
	mut.profiles <- NULL;
	if(!is.null(patients)) {
		# remove the mutation annotation and format as a matrix instead of a data frame
		if(any(!patients %in% colnames(mut.data.and.anno))) {
			warning(paste0('the following patients were not found in the given file:',paste(patients[which(!patients %in% colnames(mut.data.and.anno))],collapse=',')));
			patients <- patients[which(patients %in% colnames(mut.data.and.anno))];
			}
		mut.profiles <- matrix(as.numeric(as.matrix(mut.data.and.anno[,patients])),nrow=nrow(mut.data.and.anno));
		colnames(mut.profiles) <- patients;
		rownames(mut.profiles) <- rownames(mut.data.and.anno);
		}

	### pull out the feature annotation ###
	
	# find all the annotation columns that match the requested annotation fields
	colname.matches <- c();
	colname.repl <- c();
	for(i in annotation.fields) {
		match.idx <- grep(tolower(i), tolower(colnames(mut.data.and.anno)));

		colname.matches <- c(colname.matches,match.idx);
		if(length(match.idx) == 1) {
			colname.repl <- c(colname.repl,i);
			}
		if(length(match.idx) > 1) {
			colname.repl <- c(colname.repl,paste(rep(i,length(match.idx)),colnames(mut.data.and.anno)[match.idx],sep='.'));
			}
		}

	# select the requested annotation
	colname.matches <- unique(colname.matches);
	colname.repl <- unique(colname.repl);
	mut.anno <- NULL;
	if(length(colname.matches) > 0) {
		mut.anno <- mut.data.and.anno[,colname.matches];
		if(length(colname.matches) > 1) {
			colnames(mut.anno) <- colname.repl;
			}
		}
	else if(!is.null(annotation.fields)) {
		warning(paste('mutation.annotation.fields (',annotation.fields,') didn\'t match any of the column names. The options for ',file,' are: ',paste(colnames(mut.data.and.anno)[grep('\\d\\d\\d',colnames(mut.data.and.anno),invert=TRUE)],collapse=', '),sep=''));
		}


	### return the aberration profiles or feature annotation or both depending what was passed to the function ###
	if(!is.null(mut.profiles) & !is.null(mut.anno)) {
		return(list(aberration.profiles=mut.profiles,feature.annotation=mut.anno));
		}
	if(!is.null(mut.profiles)) {
		return(mut.profiles);
		}
	if(!is.null(mut.anno)) {
		return(mut.anno);
		}
	return(mut.data.and.anno);
	}

