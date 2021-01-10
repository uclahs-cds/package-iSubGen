\name{apply.scaling}
\alias{apply.scaling}
\title{Apply scaling factors}
\description{Apply scaling factors prior to autoencoder}
\usage{
apply.scaling(data.matrices, scaling.factors);
}
\arguments{
  \item{data.matrices}{a list of matrices. The list has one matrix for each data type to be scaled}
  \item{scaling.factors}{a list with a key \\"center\\" list and a key \\"scale\\" list}
}
\details{The keys for the data matrices and the center and scale list match. }
\value{A list of matrices of the same format as the data.matrices}
\author{Natalie Fox}
\examples{

# Load molecular profiles for three data types and calculate scaling for each
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
scaling.factors <- list(
  center = list(),
  scale = list()
  );
for(i in c('cna','csnvs','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );

  scaling.factors$center[[i]] <- apply(molecular.data[[i]], 1, mean);
  scaling.factors$scale[[i]] <- apply(molecular.data[[i]], 1, sd);
  }

# Example 1: Transform the molecular profiles by the scaling factors
scaled.molecular.data <- apply.scaling(molecular.data, scaling.factors);

# Example 2: Transform one of the data types based on the scaling factors
scaled.molecular.data2 <- apply.scaling(
  molecular.data[[1]],
    list(
      center = scaling.factors$center[[1]],
      scale = scaling.factors$scale[[1]]
      )
    );
}
