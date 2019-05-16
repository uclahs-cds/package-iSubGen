\name{calculate.integrated.correlation.matrix}
\alias{calculate.integrated.correlation.matrix}
\title{Calculate integrated correlation matrix}
\description{Calculate pairwise correlations between patient distances}
\usage{
calculate.integrated.correlation.matrix(data.types, data.matrices, dist.metrics, correlation.method = "spearman", filter.to.common.patients = FALSE,patients.to.return=NULL,patients.for.correlations=NULL)
}
\arguments{
  \item{data.types}{a vector of the ids for the different data types that are the names of the lists for the data.matrices and dist.metrics}
  \item{data.matrices}{a list of the matrices of omics features (rows) by patients (coloumns)}
  \item{dist.metrics}{a list of the distance metrics for comparing patient profiles. ex. eucliean, jaccard}
  \item{correlation.method}{Specifies the type of correlation. Options are pearson, spearman or kendall.}
  \item{filter.to.common.patients}{A boolean whether to remove patients that don't have all data types}
  \item{patients.to.return}{patients to return correlations for. If NULL all patients/columns will be used.}
  \item{patients.for.correlations}{patients to use to calculare the correlations. If NULL all patients/columns will be used.}
}
\value{matrix where rows are patients and columns are pairs of data types}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('BoutrosLab.subtype.assessment'),'/exdata/');

molecular.data <- list();
for(i in c('cna','csnvs','methy')) {
	molecular.data[[i]] <- load.molecular.aberration.data(paste0(example.molecular.data.dir,i,'_profiles.txt'),patients=c(paste0('EP00',1:9),paste0('EP0',10:30)));
}

corr.matrix <- calculate.integrated.correlation.matrix(
	data.types=names(molecular.data),
	data.matrices=molecular.data,
	dist.metrics=list(
		cna='euclidean',
		csnvs='euclidean',
		methy='euclidean'
	)
);

corr.matrix2 <- calculate.integrated.correlation.matrix(
	data.types=names(molecular.data),
	data.matrices=molecular.data,
	dist.metrics=list(
		cna='euclidean',
		csnvs='euclidean',
		methy='euclidean'
	),
	patients.to.return=paste0('EP00',1:9),
	patients.for.correlations=paste0('EP0',10:30)
);

corr.matrix3 <- calculate.integrated.correlation.matrix(
	data.types=names(molecular.data)[c(1,3)],
	data.matrices=molecular.data[c(1,3)],
	dist.metrics=list(
		cna='euclidean',
		csnvs='euclidean',
		methy='euclidean'
	)[c(1,3)],
	patients.to.return=paste0('EP00',1:9),
	patients.for.correlations=paste0('EP0',10:30)
);




}
