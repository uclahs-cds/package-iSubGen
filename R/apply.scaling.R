apply.scaling <- function(aberration.matrices, scaling.factors) {
	if(class(aberration.matrices) == 'matrix') {
		if(class(scaling.factors$center) == 'list') {
			scaling.factors$center <- scaling.factors$center[[1]];
		}
		if(class(scaling.factors$scale) == 'list') {
			scaling.factors$scale <- scaling.factors$scale[[1]];
		}	
		for(i in 1:nrow(aberration.matrices)) {
			center.adjustment <- scaling.factors$center[rownames(aberration.matrices)[i]];
			scale.adjustment <- 1;
			if(scaling.factors$scale[rownames(aberration.matrices)[i]] > 0) {
				scale.adjustment <- scaling.factors$scale[rownames(aberration.matrices)[i]];
			}
			aberration.matrices[i,] <- (aberration.matrices[i,] - center.adjustment)/scale.adjustment;
		}
		return(aberration.matrices);
	}
	for(aberration.type in names(aberration.matrices)) {
		for(i in 1:nrow(aberration.matrices[[aberration.type]])) {
			center.adjustment <- scaling.factors$center[[aberration.type]][rownames(aberration.matrices[[aberration.type]])[i]];
			scale.adjustment <- scaling.factors$scale[[aberration.type]][rownames(aberration.matrices[[aberration.type]])[i]];
			aberration.matrices[[aberration.type]][i,] <- (aberration.matrices[[aberration.type]][i,] - center.adjustment)/scale.adjustment;
		}
	}
	return(aberration.matrices);
}	
