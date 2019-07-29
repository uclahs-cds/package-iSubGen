\name{calculate.consensus.integrative.correlation.matrix}
\alias{calculate.consensus.integrative.correlation.matrix}
\title{Calculate consensus integrative correlation matrix}
\description{Calculate consensus pairwise correlations between patient distances}
\usage{
calculate.consensus.integrative.correlation.matrix(aberration.types, aberration.matrices, dist.metrics,
correlation.method = "spearman", filter.to.common.patients = FALSE, 
patients.to.return = NULL, patients.for.correlations = NULL,patient.proportion = 0.8,
feature.proportion = 0.8,num.iterations=10,print.correlation.matrices.to.file=TRUE,print.dir='.')
}
\arguments{
  \item{aberration.types}{a vector of the ids for the different aberration types that are the names of the lists for the aberration.matrices and dist.metrics}
  \item{aberration.matrices}{a list of the matrices of omics features (rows) by patients (coloumns)}
  \item{dist.metrics}{a list of the distance metrics for comparing patient profiles. ex. eucliean, jaccard}
  \item{correlation.method}{Specifies the type of correlation. Options are pearson, spearman or kendall.}
  \item{filter.to.common.patients}{A boolean whether to remove patients that don't have all aberration types}
  \item{patients.to.return}{patients to return correlations for. If NULL all patients/columns will be used.}
  \item{patients.for.correlations}{patients to use to calculare the correlations. If NULL all patients/columns will be used.}
  \item{patient.proportion}{proportion of patients.for.correlations to sample for each bootstrap}
  \item{feature.proportion}{proportion of the features to sample for each bootstrap}
  \item{num.iterations}{number of iterations to take the median from}
  \item{print.correlation.matrices.to.file}{whether to print the sampled correlation matrices to file}
  \item{print.dir}{directory for where to print the correlation matrices to file}
}
\value{matrix where rows are patients and columns are pairs of aberration types}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

molecular.data <- list();
for(i in c('cna','csnvs','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9),paste0('EP0',10:30))
  );
}

corr.matrix <- calculate.integrative.correlation.matrix(
  aberration.types = names(molecular.data),
  aberration.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    csnvs = 'euclidean',
    methy = 'euclidean'
  )
);

corr.matrix2 <- calculate.integrative.correlation.matrix(
  aberration.types = names(molecular.data),
  aberration.matrices = molecular.data,
  dist.metrics = list(
    cna = 'euclidean',
    csnvs = 'euclidean',
    methy = 'euclidean'
  ),
  patients.to.return = paste0('EP00',1:9),
  patients.for.correlations = paste0('EP0',10:30)
);

corr.matrix3 <- calculate.integrative.correlation.matrix(
  aberration.types=names(molecular.data)[c(1,3)],
  aberration.matrices=molecular.data[c(1,3)],
  dist.metrics=list(
    cna='euclidean',
    csnvs='euclidean',
    methy='euclidean'
  )[c(1,3)],
  patients.to.return=paste0('EP00',1:9),
  patients.for.correlations=paste0('EP0',10:30)
);

}
