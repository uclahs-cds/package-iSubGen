\name{combine.compressed.features}
\alias{combine.compressed.features}
\title{Combine compressed features}
\description{Combine a reduced features matrix (ex. from autoencoders) and pairwise similarity matrices into one compressed feature matrix.}
\usage{
combine.compressed.features(reduced.features, pairwise.similarity.features,
reduced.features.rescale.recenter = NA, pairwise.similarity.features.rescale.recenter = NA, 
reduced.features.rescale.denominator = NA, pairwise.similarity.features.rescale.denominator = NA,
reduced.features.weights = rep(1,ncol(reduced.features)), 
pairwise.similarity.features.weights = rep(1,ncol(pairwise.similarity.features)) )
}
\arguments{
  \item{reduced.features}{a matrix of a limited amount of multi-omics features with patients as rows and features as columns}
  \item{pairwise.similarity.features}{a matrix of a limited amount of multi-omics features with patients as rows and features as columns}
  \item{reduced.features.rescale.recenter}{either NA, "mean", a single number or a vector of numbers of length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.rescale.recenter}{either NA, "mean", a single number or a vector of numbers of length equal to the number of columns of pairwise.similarity.features}
  \item{reduced.features.rescale.denominator}{either NA, "sd", a single number or a vector of numbers of length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.rescale.denominator}{either NA, "sd", a single number or a vector of numbers of length equal to the number of columns of pairwise.similarity.features}
  \item{reduced.features.weights}{a single number or a vector of numbers of length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.weights}{a single number or a vector of numbers of length equal to the number of columns of pairwise.similarity.features}
}
\value{
  \item{compressed.feature.matrix}{a matrix of compressed features with patients as rows and features as columns}
  \item{reduced.features.rescale.recenter}{a numeric vector with length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.rescale.recenter}{a numeric vector with length equal to the number of columns of pairwise.similarity.features}
  \item{reduced.features.rescale.denominator}{a numeric vector with length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.rescale.denominator}{a numeric vector with length equal to the number of columns of pairwise.similarity.features}
  \item{reduced.features.weights}{a numeric vector with length equal to the number of columns of reduced.features}
  \item{pairwise.similarity.features.weights}{a numeric vector with length equal to the number of columns of pairwise.similarity.features}
}
\details{
The recenter values determine the how column centering is performed. If NA, no recentering is done. If the values equal "mean", then the mean of each column will be used. Otherwise, the numeric values specified will be used. The denominator values determine how column scaling is performed. If NA, no recentering is done. If the denomiator values equal "sd", then the standard deviation of each column will be used. Otherwise, the numeric values specified will be used. The values used are returned by the function along with the compressed feature matrix to be recorded for reproducibility purposes.  
}
\author{Natalie Fox}
\examples{

reduced.feature.matrix <- matrix(runif(25*4), ncol = 4);
rownames(reduced.feature.matrix) <- c(paste0('EP00',1:9), paste0('EP0',10:25));
pairwise.similarity.feature.matrix <- matrix(runif(25*6), ncol=6);
rownames(pairwise.similarity.feature.matrix) <- c(paste0('EP00',1:9), paste0('EP0',10:25));

compressed.feature.matrix <- combine.compressed.features(
  reduced.feature.matrix,
  pairwise.similarity.feature.matrix
  )$compressed.feature.matrix;

compressed.feature.matrix.rescaled.result <- combine.compressed.features(
  reduced.feature.matrix,
  pairwise.similarity.feature.matrix,
  reduced.features.rescale.recenter = 'mean',
  pairwise.similarity.features.rescale.recenter = 'mean',
  reduced.features.rescale.denominator = 'sd',
  pairwise.similarity.features.rescale.denominator = 'sd'
  );
compressed.feature.matrix.rescaled <- compressed.feature.matrix.rescaled.result$compressed.feature.matrix;

compressed.feature.matrix.reweighted.result <- combine.compressed.features(
  reduced.feature.matrix,
  pairwise.similarity.feature.matrix,
  reduced.features.weights = 1/4,
  pairwise.similarity.features.weights = 1/6
  );
compressed.feature.matrix.reweighted <- compressed.feature.matrix.reweighted.result$compressed.feature.matrix;
}
