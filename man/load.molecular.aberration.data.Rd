\name{load.molecular.aberration.data}
\alias{load.molecular.aberration.data}
\title{Load molecular aberration data}
\description{Load the molecular aberration profiles/feature annotation}
\usage{
load.molecular.aberration.data(file, patients = NULL, annotation.fields = NULL);
}
\arguments{
  \item{file}{the file name of the matrix containing molecular and annotation data}
  \item{patients}{A vector of patients ids. Must match columnames from aberration file}
  \item{annotation.fields}{A vector referencing the column names for the feature annotation columns}
}
\details{The annotation.fields argument will look for any colnames which contain the values specificed in annotation.fields and then the column will be renamed to the value that matched from annotation.fields.}
\value{If the patients argument is specified then the patient molecular aberration profiles are returned. If the annotaiton.fields argument is specified then the feature annotation is returned. If both are specified then the two matrices are returned in a list. If neither is specified then the entire matrix with the mix of patients and annotation is returned.}
\author{Natalie Fox}
\examples{

example.aberration.data <- paste0(
  path.package('iSubGen'),
  '/exdata/cna_profiles.txt'
);

cna.profiles <- load.molecular.aberration.data(
  example.aberration.data,
  patients = c(paste0('EP00',1:9),paste0('EP0',10:30))
);

cna.annotation <- load.molecular.aberration.data(
  example.aberration.data,
  annotation.fields = c('gene','start','end')
);
}
