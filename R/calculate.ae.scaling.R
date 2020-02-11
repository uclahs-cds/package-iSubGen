calculate.ae.scaling <- function(aberration.matrices) {
	scale.factors.center <- list();
	scale.factors.scale <- list();
	for(aberration.type in names(aberration.matrices)) {
		scale.factors.center[[aberration.type]] <- apply(aberration.matrices[[aberration.type]],1,mean);
		scale.factors.scale[[aberration.type]] <- apply(aberration.matrices[[aberration.type]],1,sd);
        }
	return(list(center=scale.factors.center,scale=scale.factors.scale));
}
