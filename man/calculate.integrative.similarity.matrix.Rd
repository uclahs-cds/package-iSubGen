\name{calculate.integrative.similarity.matrix}
\alias{calculate.integrative.similarity.matrix}
\title{Calculate integrative similarity matrix}
\description{Calculate pairwise correlations between patient distances}
\usage{
calculate.integrative.similarity.matrix(data.types, data.matrices, dist.metrics,
correlation.method = "spearman", filter.to.common.patients = FALSE, 
patients.to.return = NULL, patients.for.correlations = NULL)
}
\arguments{
  \item{data.types}{vector, where each element is a data type ID matching the names in data.matrices and dist.metrics}
  \item{data.matrices}{list, where each element is a matrix with features as rows and patients as columns}
  \item{dist.metrics}{list, where each element is the distance metric to use for comparing patient profiles. ex. euclidean. Options are from philentropy::distance}
  \item{correlation.method}{specifies the type of correlation. Options are pearson, spearman or kendall.}
  \item{filter.to.common.patients}{logical, where TRUE indicates to filter out patients that don't have all data types}
  \item{patients.to.return}{vector, where each element a patient ID specifying the patients to calculate integrative similarity for. For example, this is the testing cohort patients when calculating integrative similarity for the testing cohort using the training cohort patients. If NULL all patients/columns will be used.}
  \item{patients.for.correlations}{vector, where each element a patient ID specifying the patients to use to calculate the similarities. For example, this would be the training cohort patients when calculating integrative similarity for the testing cohort. If NULL all patients/columns will be used.}
}
\value{matrix where rows are patients and columns are pairs of data types}
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

# Example 1: calculate integrative similarity between pairs of CNA, coding SNVs, methylation data
corr.matrix <- calculate.integrative.similarity.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    snv = 'euclidean',
    methy = 'euclidean'
    )
  );

# Example 2: calculate the integrative similarity for patients EP001 through EP009 
# in relation to patients EP010 through EP030 meaning the profile of EP001 is 
# correlated to the profiles of EP010 through EP030 so when assessing new patients,
# they can be compared to the training profiles
corr.matrix2 <- calculate.integrative.similarity.matrix(
  data.types = names(molecular.data),
  data.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    snv = 'euclidean',
    methy = 'euclidean'
    ),
  patients.to.return = paste0('EP00',1:9),
  patients.for.correlations = paste0('EP0',10:30)
  );

# Example 3: Calculate integrative similarity between CNA and methylation data
corr.matrix3 <- calculate.integrative.similarity.matrix(
  data.types=names(molecular.data)[c(1,3)],
  data.matrices=molecular.data[c(1,3)],
  dist.metrics=list(
    cna='euclidean',
    snv='euclidean',
    methy='euclidean'
    )[c(1,3)],
  patients.to.return=paste0('EP00',1:9),
  patients.for.correlations=paste0('EP0',10:30)
  );
}
