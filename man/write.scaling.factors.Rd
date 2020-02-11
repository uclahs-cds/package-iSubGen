\name{write.scaling.factors}
\alias{write.scaling.factors}
\title{Write scaling factors to file}
\description{Write scaling factors to file}
\usage{
write.scaling.factors(scaling.factors,scaling.factor.files.dir)
}
\arguments{
  \item{scaling.factors}{a list with the scaling factors created by calculate.ae.scaling}
  \item{scaling.factor.files.dir}{directory to output scaling factor files}
}
\details{Creates two files for each aberration type key. One file for the recentering values and one file for the rescaling values. Files have the names <aberration type>_gene_recenter.txt or <aberration type>_gene_rescale.txt}
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

scaling.factors <- calculate.ae.scaling(molecular.data);

write.scaling.factors(scaling.factors,example.molecular.data.dir);
}
