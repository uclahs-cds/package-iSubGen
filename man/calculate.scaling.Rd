\name{calculate.ae.scaling}
\alias{calculate.ae.scaling}
\title{Calculate autoencoder scaling factors}
\description{Calculate autoencoder scaling factors}
\usage{
calculate.ae.scaling(aberration.matrices);
}
\arguments{
  \item{aberration.matrices}{a list of matrices. One for each aberration type to be scaled}
}
\details{The keys for the aberration matrices and the center and scale list match}
\value{a list with a key \\"center\\" list and a key \\"scale\\" list. The center and scale list keys match the aberration.matrices list keys}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

molecular.data <- list();
for(i in c('cna','csnvs','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9),paste0('EP0',10:30))
  );
}

scaling.factors <- calculate.ae.scaling(molecular.data);

scaling.factors2 <- calculate.ae.scaling(molecular.data[['methy']]);
}
