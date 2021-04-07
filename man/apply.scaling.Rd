\name{apply.scaling}
\alias{apply.scaling}
\title{Apply scaling factors}
\description{Apply scaling factors prior to autoencoder}
\usage{
apply.scaling(data.matrices, scaling.factors);
}
\arguments{
  \item{data.matrices}{list, where each element is a matrix. The list has one matrix for each data type to be scaled}
  \item{scaling.factors}{list with two elements named: \\"center\\" and \\"scale\\", and each element is a named numerical vector or a list of named numerical vectors. If scaling.factors$center or scaling.factors$scale are a list then each element needs to correspond to a one of the data matrices. Finally, the named numerical vectors should match the row and rownames from the corresponding data matrix.}
}
\details{The names for the data matrices and the center and scale lists all must match.}
\value{A list of matrices of the same format as the data.matrices}
\author{Natalie Fox}
\examples{

# Load molecular profiles for three data types and calculate scaling for each
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
scaling.factors <- list();
for(i in c('cna','snv','methy')) {
  # Load molecular profiles from example files saved 
  # in the package as <data type>_profiles.txt
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );

  scaling.factors[[i]] <- list();

  scaling.factors[[i]]$center <- apply(molecular.data[[i]], 1, mean);
  scaling.factors[[i]]$scale <- apply(molecular.data[[i]], 1, sd);
  }

# Example 1: Transform the molecular profiles by the scaling factors
scaled.molecular.data <- apply.scaling(molecular.data, scaling.factors);

# Example 2: Transform one of the data types based on the scaling factors
scaled.molecular.data2 <- apply.scaling(
  molecular.data[[1]],
  scaling.factors[[1]]
  );
}
