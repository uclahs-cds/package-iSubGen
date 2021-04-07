\name{create.autoencoder.irf.matrix}
\alias{create.autoencoder.irf.matrix}
\title{Create matrix of independent reduced features}
\description{Create matrix of independent reduced features using autoencoders}
\usage{
create.autoencoder.irf.matrix(data.types, data.matrices,
autoencoders, filter.to.common.patients = FALSE,
patients.to.return = NULL)
}
\arguments{
  \item{data.types}{vector, where each element is a data type ID matching the names in data.matrices and dist.metrics}
  \item{data.matrices}{list, where each element is a matrix with features as rows and patients as columns}
  \item{autoencoders}{list, where each element is an autoencoder corresponding to each data type. Can be either an keras autoencoder object or the file where the autoencoder was saved.}
  \item{filter.to.common.patients}{logical, where TRUE indicates to filter out patients that don't have all data types.}
  \item{patients.to.return}{vector of patients to return correlations for. If NULL all patients/columns will be used.}
}
\value{matrix where rows are patients and columns are pairs of data types}
\author{Natalie Fox}
\examples{
\dontrun{

# Load three data types and create an autoencder for each
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
ae.result <- list();
for(i in c('cna','snv','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );
  ae.result[[i]] <- create.autoencoder(
    data.type = i,
    data.matrix = molecular.data[[i]],
    encoder.layers.node.nums = c(10,2)
    )$autoencoder;
  }

# Create a matrix of the bottleneck layers 
irf.matrix <- create.autoencoder.irf.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  autoencoders = ae.result
  );
}
}
