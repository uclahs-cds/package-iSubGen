\name{load.molecular.aberration.data}
\alias{load.molecular.aberration.data}
\title{Load molecular aberration data}
\description{Load the molecular aberration profiles/feature annotation}
\usage{
load.molecular.aberration.data(file, patients = NULL, annotation.fields = NULL);
}
\arguments{
  \item{file}{file name of the matrix containing molecular and annotation data. If it does not contain an _absolute_ path, the file name is _relative_ to the current working directory, 'getwd()' as in read.table.}
  \item{patients}{vector of patients IDs. Must match colnames from aberration file}
  \item{annotation.fields}{vector referencing the column names for the feature annotation columns}
}
\details{The annotation.fields argument will look for any colnames which contain the values specified in annotation.fields and then the column will be renamed to the value that matched from annotation.fields.}
\value{If the patients argument is specified then the patient molecular aberration profiles are returned. If the annotation.fields argument is specified then the feature annotation is returned. If both are specified then the two matrices are returned in a list. If neither is specified then the entire matrix with the mix of patients and annotation is returned.}
\author{Natalie Fox}
\examples{

example.aberration.data <- paste0(
  path.package('iSubGen'),
  '/exdata/cna_profiles.txt'
  );

# Load the CNA profiles for patients EP001 through EP030
cna.profiles <- load.molecular.aberration.data(
  example.aberration.data,
  patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
  );

# Load feature annotation for the CNA data
cna.annotation <- load.molecular.aberration.data(
  example.aberration.data,
  annotation.fields = c('gene','start','end')
  );
}
