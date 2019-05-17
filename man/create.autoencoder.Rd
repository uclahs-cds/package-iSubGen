\name{create.autoencoder}
\alias{create.autoencoder}
\title{Create an autoencoder for dimensionality reduction}
\description{Create an autoencoder for dimensionality reduction using keras and tensorflow packages}
\usage{
create.autoencoder(data.type, data.matrix, decoding.layer.node.num=c(15,2), autoencoder.activation='tanh', optimization.loss.function='mean_squared_error', model.file.output.dir='.')
}
\arguments{
  \item{data.type}{the id of the different data type -  will be used for naming the output file}
  \item{data.matrix}{a matrix of omics features (rows) by patients (coloumns)}
  \item{decoding.layer.node.num}{a vector of the number of nodes for each layer in the reducing the feature dimensions within th autoencoder. The autoencoder will be made symmetrically so the number of nodes in each layer will be used in reverse, not repeating the last layer to re encode the features in the autoencoder} 
  \item{autoencoder.activation}{the type of activation function to use in the autoencoder}
  \item{optimization.loss.function}{the loss function used for optimization while fitting the autoencoder}
  \item{model.file.output.dir}{the file location for outputing the autoencoder file}
}
\value{
  \item{autoencoder}{the autoencoder created by the keras package}
  \item{autoencoder.file}{the hdf5 file that the model was saved in and can be loaded from}
}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

 ae.result <- create.autoencoder(
	data.type='cna',
	data.matrix=load.molecular.aberration.data(
		paste0(example.molecular.data.dir,'cna_profiles.txt'),
		patients=c(paste0('EP00',1:9),paste0('EP0',10:30))
	)
);

}
