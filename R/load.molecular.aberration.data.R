load.molecular.aberration.data <- function(
	file,
	patients = NULL,
	annotation.fields = NULL
	) {

	# read in matrix with patient annotation and aberration data
	aberration.data.and.anno <- read.table(
		file,
		sep = '\t',
		header = TRUE
		);

	### pull out the aberration data ###
	aberration.profiles <- NULL;
	if (!is.null(patients)) {
		# select aberration data and format as a matrix instead of a data frame
		if (any(!patients %in% colnames(aberration.data.and.anno))) {
			warning(paste0(
				'the following patients were not found in the given file:',
				paste(
					patients[!patients %in% colnames(aberration.data.and.anno)],
					collapse = ','
					)
				));
			patients <- patients[patients %in% colnames(aberration.data.and.anno)];
			}
		aberration.profiles <- matrix(
			data = as.numeric(as.matrix(aberration.data.and.anno[,patients])), 
			nrow = nrow(aberration.data.and.anno)
			);
		colnames(aberration.profiles) <- patients;
		rownames(aberration.profiles) <- rownames(aberration.data.and.anno);
		}

	### pull out the feature annotation ###
	
	# find all the annotation columns that match the requested annotation fields
	colname.matches <- c();
	colname.repl <- c();
	for(i in annotation.fields) {
		match.idx <- grep(tolower(i), tolower(colnames(aberration.data.and.anno)));

		colname.matches <- c(colname.matches, match.idx);
		if (length(match.idx) == 1) {
			colname.repl <- c(colname.repl, i);
			}
		if (length(match.idx) > 1) {
			colname.repl <- c(
				colname.repl, 
				paste(
					rep(i,length(match.idx)),
					colnames(aberration.data.and.anno)[match.idx],
					sep = '.'
					)
				);
			}
		}

	# select the requested annotation
	colname.matches <- unique(colname.matches);
	colname.repl <- unique(colname.repl);
	aberration.anno <- NULL;
	if (length(colname.matches) > 0) {
		aberration.anno <- aberration.data.and.anno[,colname.matches];
		if (length(colname.matches) > 1) {
			colnames(aberration.anno) <- colname.repl;
			}
		}
	else if (!is.null(annotation.fields)) {
		warning(paste(
			'annotation.fields (',annotation.fields,') didn\'t match any of the column names. The options for ',
			file,' are: ', paste(colnames(aberration.data.and.anno)[grep('\\d\\d\\d', colnames(aberration.data.and.anno),
			invert=TRUE)], collapse=', '),
			sep = ''
			));
		}

	### return the aberration data or feature annotation or both depending what was passed to the function ###
	if (!is.null(aberration.profiles) & !is.null(aberration.anno)) {
		return(list( aberration.profiles = aberration.profiles, feature.annotation = aberration.anno));
		}
	if (!is.null(aberration.profiles)) {
		return(aberration.profiles);
		}
	if (!is.null(aberration.anno)) {
		return(aberration.anno);
		}
	return(aberration.data.and.anno);
	}
