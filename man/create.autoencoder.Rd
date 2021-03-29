\name{create.autoencoder}
\alias{create.autoencoder}
\title{Create an autoencoder for dimensionality reduction}
\description{Create an autoencoder for dimensionality reduction using keras and tensorflow packages}
\usage{
create.autoencoder(data.type, data.matrix, encoder.layers.node.nums = c(15,2),
autoencoder.activation = 'tanh', optimization.loss.function = 'mean_squared_error', 
model.file.output.dir = '.')
}
\arguments{
  \item{data.type}{data type ID. The ID will be used for naming the output file}
  \item{data.matrix}{matrix with data features as rows and patients as columns}
  \item{encoder.layers.node.nums}{vector with the number of nodes for each layer when the reducing the feature dimensions within the autoencoder. The autoencoder will be made symmetrically so the number of nodes in each layer will be used in reverse, not repeating the last layer to re encode the features in the autoencoder} 
  \item{autoencoder.activation}{activation function to use in the autoencoder}
  \item{optimization.loss.function}{loss function used for optimization while fitting the autoencoder}
  \item{model.file.output.dir}{file location for the autoencoder file}
}
\value{
  \item{autoencoder}{the autoencoder created by the keras package}
  \item{autoencoder.file}{the hdf5 file that the model was saved in and can be loaded from}
}
\author{Natalie Fox}
\examples{
\dontrun{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

ae.result <- create.autoencoder(
  data.type = 'cna',
  data.matrix = load.molecular.aberration.data(
    paste0(example.molecular.data.dir,'cna_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    ),
  encoder.layers.node.nums = c(15,5,2)
  );
}
}
