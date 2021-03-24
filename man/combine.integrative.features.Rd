\name{combine.integrative.features}
\alias{combine.integrative.features}
\title{Combine iSubGen integrative features}
\description{Combine a independent reduced features matrix (ex. from autoencoders) and pairwise integrative similarity matrices into one integrative feature matrix.}
\usage{
combine.integrative.features(irf.matrix, cis.matrix,
irf.rescale.recenter = NA, cis.rescale.recenter = NA, 
irf.rescale.denominator = NA, cis.rescale.denominator = NA,
irf.weights = rep(1, ncol(irf.matrix)), 
cis.weights = rep(1, ncol(cis.matrix)))
}
\arguments{
  \item{irf.matrix}{matrix of independent reduced features with patients as rows and features as columns}
  \item{cis.matrix}{matrix of consensus integrative similarity or integrative similarity features with patients as rows and features as columns}
  \item{irf.rescale.recenter}{either NA, "mean", a single number or a vector of numbers of length equal to the number of columns of irf}
  \item{cis.rescale.recenter}{either NA, "mean", a single number or a vector of numbers of length equal to the number of columns of cis}
  \item{irf.rescale.denominator}{either NA, "sd", a single number or a vector of numbers of length equal to the number of columns of irf}
  \item{cis.rescale.denominator}{either NA, "sd", a single number or a vector of numbers of length equal to the number of columns of cis}
  \item{irf.weights}{single number or vector of numbers of length equal to the number of columns of irf}
  \item{cis.weights}{single number or vector of numbers of length equal to the number of columns of cis}
}
\value{
  \item{integrative.feature.matrix}{a matrix of compressed features with patients as rows and features as columns}
  \item{irf.rescale.recenter}{a numeric vector with length equal to the number of columns of irf}
  \item{cis.rescale.recenter}{a numeric vector with length equal to the number of columns of cis}
  \item{irf.rescale.denominator}{a numeric vector with length equal to the number of columns of irf}
  \item{cis.rescale.denominator}{a numeric vector with length equal to the number of columns of cis}
  \item{irf.weights}{a numeric vector with length equal to the number of columns of irf}
  \item{cis.weights}{a numeric vector with length equal to the number of columns of cis}
}
\details{
The recenter values determine the how column centering is performed. If NA, no recentering is done. If the values equal "mean", then the mean of each column will be used. Otherwise, the numeric values specified will be used. The denominator values determine how column scaling is performed. If NA, no recentering is done. If the denominator values equal "sd", then the standard deviation of each column will be used. Otherwise, the numeric values specified will be used. The values used are returned by the function along with the compressed feature matrix to be recorded for reproducibility purposes.  
}
\author{Natalie Fox}
\examples{

# Create matrices for combining
irf.matrix <- matrix(runif(25*4), ncol = 4);
rownames(irf.matrix) <- c(paste0('EP00',1:9), paste0('EP0',10:25));
cis.matrix <- matrix(runif(25*6), ncol=6);
rownames(cis.matrix) <- c(paste0('EP00',1:9), paste0('EP0',10:25));

# Example 1: Join the matrices without any weighting adjustments
isubgen.feature.matrix <- combine.integrative.features(
  irf.matrix,
  cis.matrix
  )$integrative.feature.matrix;

# Example 2: Combine matrices after scaling each column by subtracting the mean
# and dividing by the standard devation of the column
isubgen.feature.matrix.rescaled.result <- combine.integrative.features(
  irf.matrix,
  cis.matrix,
  irf.rescale.recenter = 'mean',
  cis.rescale.recenter = 'mean',
  irf.rescale.denominator = 'sd',
  cis.rescale.denominator = 'sd'
  );
isubgen.feature.matrix.2 <- isubgen.feature.matrix.rescaled.result$integrative.feature.matrix;

# Example 3: Combine matrices 
isubgen.feature.matrix.reweighted.result <- combine.integrative.features(
  irf.matrix,
  cis.matrix,
  irf.weights = 1/4,
  cis.weights = 1/6
  );
isubgen.feature.matrix.3 <- isubgen.feature.matrix.reweighted.result$integrative.feature.matrix;
}
