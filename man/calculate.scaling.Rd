\name{calculate.scaling}
\alias{calculate.scaling}
\title{Calculate scaling factors}
\description{Calculate scaling factors}
\usage{
calculate.scaling(data.matrices);
}
\arguments{
  \item{data.matrices}{a list of matrices. The list has one matrix for each data type to be scaled}
}
\details{The keys for the data.matrices and the center and scale list match}
\value{a list with a key \\"center\\" list and a key \\"scale\\" list. The center and scale list keys match the data.matrices list keys}
\author{Natalie Fox}
\examples{

# Load molecular profiles for three data types
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
for(i in c('cna','csnvs','methy')) {
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
