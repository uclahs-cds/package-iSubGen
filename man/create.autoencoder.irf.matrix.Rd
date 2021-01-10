\name{create.autoencoder.irf.matrix}
\alias{create.autoencoder.irf.matrix}
\title{Create matrix of compressed autoencoder features}
\description{Create matrix of compressed autoencoder features}
\usage{
create.autoencoder.irf.matrix(data.types, data.matrices,
autoencoders, number.of.layers.in.decode = 2, filter.to.common.patients = FALSE,
patients.to.return = NULL)
}
\arguments{
  \item{data.types}{a vector of the ids for the different data types that are the names of the lists for the data.matrices and dist.metrics}
  \item{data.matrices}{a list of the matrices with features (rows) by patients (columns)}
  \item{autoencoders}{a list of an autoencoder corresponding to each data type. Can be either an keras autoencoder object or the file where the autoencoder was saved.}
  \item{number.of.layers.in.decode}{number of layers to decode in the autoencoder}
  \item{filter.to.common.patients}{a boolean whether to remove patients that don't have all data types}
  \item{patients.to.return}{a vector patients to return correlations for. If NULL all patients/columns will be used.}
}
\value{matrix where rows are patients and columns are pairs of data types}
\author{Natalie Fox}
\examples{

# Load three data types and create an autoencder for each
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
ae.result <- list();
for(i in c('cna','csnvs','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );
  ae.result[[i]] <- create.autoencoder(
    data.type = i,
    data.matrix = molecular.data[[i]]
    )$autoencoder;
  }

# Create a matrix of the bottleneck layers 
compressed.feature.matrix <- create.autoencoder.irf.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  autoencoders = ae.result
  );
}
