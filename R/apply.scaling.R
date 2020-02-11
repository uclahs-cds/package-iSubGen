apply.scaling <- function(aberration.matrices, scaling.factors) {
	for(aberration.type in names(aberration.matrices)) {
		for(i in 1:nrow(aberration.matrices[[aberration.type]])) {
			center.adjustment <- scaling.factors$center[[aberration.type]][rownames(aberration.matrices[[aberration.type]])[i]];
			scale.adjustment <- scaling.factors$scale[[aberration.type]][rownames(aberration.matrices[[aberration.type]])[i]];
			aberration.matrices[[aberration.type]][i,] <- (aberration.matrices[[aberration.type]][i,] - center.adjustment)/scale.adjustment;
		}
	}
	return(aberration.matrices);
}	
