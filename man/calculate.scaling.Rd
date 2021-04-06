\name{calculate.scaling}
\alias{calculate.scaling}
\title{Calculate scaling factors}
\description{Calculate scaling factors}
\usage{
calculate.scaling(data.matrices);
}
\arguments{
  \item{data.matrices}{list, where each element is a matrix. The list has one matrix for each data type to be scaled}
}
\details{The names for the data matrices and the center and scale lists all must match. }
\value{a list with two elements named: \\"center\\" and \\"scale\\", and each of these element is a named numerical vector or a list of named numerical vectors. If scaling.factors$center or scaling.factors$scale are a list then each element will correspond to a one of the data matrices. Finally, the named numerical vectors will match the row and rownames from the data matrices.}
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

# Example 1: Calculate scaling factors for all three data types
scaling.factors <- calculate.scaling(molecular.data);

# Example 2: Calculate scaling factors for only the methylation data
scaling.factors2 <- calculate.scaling(molecular.data[['methy']]);
}
