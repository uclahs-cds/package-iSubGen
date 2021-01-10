calculate.scaling <- function(data.matrices) {
	if(class(data.matrices) == 'matrix') {
		return(list(
			center = apply(data.matrices,1,mean),
			scale = apply(data.matrices,1,sd)
			));
		}
	scale.factors.center <- list();
	scale.factors.scale <- list();
	for(data.type in names(data.matrices)) {
		scale.factors.center[[data.type]] <- apply(data.matrices[[data.type]],1,mean);
		scale.factors.scale[[data.type]] <- apply(data.matrices[[data.type]],1,sd);
		}
	return(list(center = scale.factors.center, scale = scale.factors.scale));
	}
