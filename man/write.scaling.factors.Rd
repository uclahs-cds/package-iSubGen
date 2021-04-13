\name{write.scaling.factors}
\alias{write.scaling.factors}
\title{Write scaling factors to file}
\description{Write scaling factors to file}
\usage{
write.scaling.factors(scaling.factors, scaling.factor.files.dir)
}
\arguments{
  \item{scaling.factors}{list with the scaling factors created by calculate.scaling}
  \item{scaling.factor.files.dir}{directory to output scaling factor files}
}
\details{Creates two files for each data type key. One file for the recentering values and one file for the rescaling values. Files have the names <data type>_gene_recenter.txt or <data type>_gene_rescale.txt}
\value{No return value, called for side effects}
\author{Natalie Fox}
\examples{

# load the aberration profiles for three data types 
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');
molecular.data <- list();
for(i in c('cna','snv','methy')) {
  molecular.data[[i]] <- load.molecular.aberration.data(
    paste0(example.molecular.data.dir,i,'_profiles.txt'),
    patients = c(paste0('EP00',1:9), paste0('EP0',10:30))
    );
  }

# calculate scaling factors for all three data types
scaling.factors <- calculate.scaling(molecular.data);

# save the scaling factors to file
write.scaling.factors(scaling.factors, example.molecular.data.dir);
}
