\name{subtype.by.h.clustering}
\alias{subtype.by.h.clustering}
\title{Subtype by heirarchial clustering}
\description{Runs consensus cluster using heirachial clustering to get patient groups based on the -omics data passed in.}
\usage{
subtype.by.h.clustering(aberration.matrix, distance.metric, parent.output.dir,
new.result.dir, subtype.table.file = NULL, renumbered.subtype.table.file = NULL,
max.num.subtypes = 12, clustering.reps = 1000, pFeature = 0.8, pItem = 0.8,
verbose = FALSE, consensus.cluster.write.table = TRUE);
}
\arguments{
  \item{aberration.matrix}{A matrix of patients (rows) by features (columns)}
  \item{distance.metric}{The distance metric for comparing patient profiles. ex. eucliean, jaccard}
  \item{parent.output.dir}{The directory where the consensus clustering function will create a directory of results}
  \item{new.result.dir}{The name of the directory for results that the consensus clustering function will create}
  \item{subtype.table.file}{Filename for subtype assignment table for different number of clusters}
  \item{renumbered.subtype.table.file}{Filename for subtype assignment table for classification that have been reordered to associate with aberration.matrix averages by subtype}
  \item{max.num.subtypes}{The maximum number of clusters to separate patients into}
  \item{clustering.reps}{Number of subsamples for consensus clustering function}
  \item{pFeature}{the proportion of features to sample for each clustering iteration}
  \item{pItem}{the proportion of patients to sample for each clustering iteration}
  \item{verbose}{a boolean specifying whether to print messages to the screen to indicate progress}
  \item{consensus.cluster.write.table}{a boolean specifying consensus cluster output}
}
\value{
  \item{consensus_cluster_result}{the consensus clustering function return value}
  \item{subtype_table}{the table written to subtype.table.file}
}
\author{Natalie Fox}
\examples{

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

subtypeing.results <- subtype.by.h.clustering(
  aberration.matrix = ex.matrix,
  distance.metric = 'jaccard',
  parent.output.dir = './',
  new.result.dir = 'example_subtyping',
  subtype.table.file = 'example_subtyping_assignment.txt',
  max.num.subtypes = 6,
  clustering.reps = 50
);
	
}
