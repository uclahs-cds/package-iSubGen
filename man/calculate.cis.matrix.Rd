\name{calculate.cis.matrix}
\alias{calculate.cis.matrix}
\title{Calculate consensus integrative correlation matrix}
\description{Calculate consensus pairwise correlations between patient distances}
\usage{
calculate.cis.matrix(data.types, data.matrices, dist.metrics,
correlation.method = "spearman", filter.to.common.patients = FALSE, 
patients.to.return = NULL, patients.for.correlations = NULL, 
patient.proportion = 0.8, feature.proportion = 1, num.iterations = 10, 
print.intermediary.similarity.matrices.to.file = TRUE, print.dir = '.',
patient.proportion.seeds = seq(1,num.iterations), 
feature.proportion.seeds = seq(1,num.iterations))
}
\arguments{
  \item{data.types}{vector of the IDs for the different data types that are the names of the lists for the data.matrices and dist.metrics}
  \item{data.matrices}{list of the matrices with features (rows) by patients (columns)}
  \item{dist.metrics}{list of the distance metrics for comparing patient profiles. ex. euclidean. Options are from philentropy::distance}
  \item{correlation.method}{specifies the type of correlation for similarity comparison. Options are pearson, spearman or kendall.}
  \item{filter.to.common.patients}{logical, where TRUE indicates to filter out patients that don't have all data types}
  \item{patients.to.return}{vector of patients to calculate CIS for. For example, this is the testing cohort patients when calculating CIS for the testing cohort using the training cohort patients. If NULL all patients/columns will be used.}
  \item{patients.for.correlations}{vector of patients to use to calculate the similarities. For example, this would be the training cohort patients when calculating CIS for the testing cohort.  If NULL all patients/columns will be used.}
  \item{patient.proportion}{proportion of patients.for.correlations to sample for each iteration (sampled without replacement).}
  \item{feature.proportion}{proportion of the features to sample for each iteration (sampled without replacement).}
  \item{num.iterations}{number of iterations to take the median from}
  \item{print.intermediary.similarity.matrices.to.file}{logical, where TRUE indicates that created intermediary integrative similarity matrix from each iteration should be printed to file}
  \item{print.dir}{directory for where to print the intermediary similarity matrices to file}
  \item{patient.proportion.seeds}{vector of scalars of the length num.iterations specifying the seeds used for random sampling for selecting the patient subsets at each iteration}
  \item{feature.proportion.seeds}{vector of scalars of the length num.iterations specifying the seeds used for random sampling for selecting the feature subsets at each iteration}
}
\value{CIS matrix where rows are patients and columns are pairs of data types}
\author{Natalie Fox}
\examples{

# Load molecular profiles for three data types from example files saved 
# in the package as <data type>_profiles.txt
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
for(i in c('cna','snv','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );
  }

# Example 1: calculate the consensus integrative similarity (CIS) matrix
corr.matrix <- calculate.cis.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    snv = 'euclidean',
    methy = 'euclidean'
    ),
  print.intermediary.similarity.matrices.to.file = FALSE
  );

# Example 2: calculate the CIS matrix for patients EP001 through EP009 in relation 
# to patients EP010 through EP030 meaning the profile of EP001 is correlated to 
# the profiles of EP010 through EP030 so when assessing new patients, they can be 
# compared to the training profiles
corr.matrix2 <- calculate.cis.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    snv = 'euclidean',
    methy = 'euclidean'
    ),
  patients.to.return = paste0('EP00',1:9),
  patients.for.correlations = paste0('EP0',10:30),
  print.intermediary.similarity.matrices.to.file = FALSE
  );

# Example 3: Adjusting the proportion of the features that will be used to correlate 
# the patient profiles
corr.matrix3 <- calculate.cis.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    snv = 'euclidean',
    methy = 'euclidean'
    ),
  patients.to.return = paste0('EP00',1:9),
  patients.for.correlations = paste0('EP0',10:30),
  feature.proportion = 0.6,
  print.intermediary.similarity.matrices.to.file = FALSE
  );

}
