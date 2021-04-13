\name{cluster.patients}
\alias{cluster.patients}
\title{Clustering to find patient subtypes}
\description{A wrapper function for using consensus clustering to subtype patients}
\usage{
cluster.patients(data.matrix, distance.metric, parent.output.dir,
new.result.dir, subtype.table.file = NULL, max.num.subtypes = 12, 
clustering.reps = 1000, proportion.features = 0.8, proportion.patients = 0.8, 
verbose = FALSE, consensus.cluster.write.table = TRUE);
}
\arguments{
  \item{data.matrix}{matrix with patients as rows and features as columns}
  \item{distance.metric}{distance metric for comparing patient profiles. ex. euclidean}
  \item{parent.output.dir}{directory where the consensus clustering function will create a directory of results}
  \item{new.result.dir}{directory name for consensus clustering results}
  \item{subtype.table.file}{filename for subtype assignment table for different number of clusters}
  \item{max.num.subtypes}{maximum number of clusters to separate patients into}
  \item{clustering.reps}{number of subsamples for consensus clustering function}
  \item{proportion.features}{proportion of features to sample for each clustering iteration}
  \item{proportion.patients}{proportion of patients to sample for each clustering iteration}
  \item{verbose}{logical, where TRUE indicates to print messages to the screen to indicate progress}
  \item{consensus.cluster.write.table}{logical, where TRUE indicates for the ConsensusClusterPlus function to writeTable}
}
\value{
  \item{consensus_cluster_result}{consensus clustering function return value}
  \item{subtype_table}{the table written to subtype.table.file}
}
\author{Natalie Fox}
\examples{
\dontrun{

# For this example instead of clustering CIS and IRF matrices,
# create a data matrix to see how the function works without
# running through the whole iSubGen process.
# This example is created with to have 4 distinct clusters
set.seed(5);
ex.matrix <- matrix(
  c(
    sample(c(0,1), 30, replace = TRUE), rep(1,75), rep(0,25),
    sample(c(0,1), 30, replace = TRUE), rep(1,75), rep(0,25),
    sample(c(0,1), 30, replace = TRUE), rep(1,75), rep(0,25),
    sample(c(0,1), 30, replace = TRUE), rep(1,100),
    sample(c(0,1), 30, replace = TRUE), rep(1,100),
    sample(c(0,1), 30, replace = TRUE), rep(1,100),
    sample(c(0,1), 30, replace = TRUE), rep(0,100),
    sample(c(0,1), 30, replace = TRUE), rep(0,100),
    sample(c(0,1), 30, replace = TRUE), rep(0,100),
    sample(c(0,1), 30, replace = TRUE), rep(0,75), rep(1,25),
    sample(c(0,1), 30, replace = TRUE), rep(0,75), rep(1,25),
    sample(c(0,1), 30, replace = TRUE), rep(0,75), rep(1,25)
    ),
  nrow=130);
rownames(ex.matrix) <- paste0('gene',1:130);
colnames(ex.matrix) <- paste0('patient',LETTERS[1:12]); 

# Use Consensus clustering to subtype the patient profiles
subtyping.results <- cluster.patients(
  data.matrix = ex.matrix,
  distance.metric = 'euclidean',
  parent.output.dir = './',
  new.result.dir = 'example_subtyping',
  max.num.subtypes = 6,
  clustering.reps = 50,
  consensus.cluster.write.table = FALSE
  );	
}
}
