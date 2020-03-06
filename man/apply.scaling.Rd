\name{apply.scaling}
\alias{apply.scaling}
\title{Apply scaling factors}
\description{Apply scaling factors prior to autoencoder}
\usage{
apply.scaling(aberration.matrices, scaling.factors);
}
\arguments{
  \item{aberration.matrices}{a list of matrices. One for each aberration type to be scaled}
  \item{scaling.factors}{a list with a key \\"center\\" list and a key \\"scale\\" list}
}
\details{The keys for the aberration matrices and the center and scale list match. }
\value{A list of matrices of the same format as the aberration.matrices}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

molecular.data <- list();
scaling.factors <- list(
  center=list(),
  scale=list()
);
for(i in c('cna','csnvs','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9),paste0('EP0',10:30))
  );

  scaling.factors$center[[i]] <- apply(molecular.data[[i]],1,mean);
  scaling.factors$scale[[i]] <- apply(molecular.data[[i]],1,sd);
}

scaled.molecular.data <- apply.scaling(molecular.data,scaling.factors);

scaled.molecular.data2 <- apply.scaling(
	molecular.data[[1]],
	list(
		center=scaling.factors$center[[1]],
		scale=scaling.factors$scale[[1]]
	)
);
}
