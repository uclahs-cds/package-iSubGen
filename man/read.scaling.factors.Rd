\name{read.scaling.factors}
\alias{read.scaling.factors}
\title{Read scaling factors from file}
\description{Read scaling factors from file}
\usage{
read.scaling.factors(scaling.factor.files.dir,aberration.types);
}
\arguments{
  \item{scaling.factor.files.dir}{the directory where the files were saved}
  \item{aberration.types}{a vactor of the aberration types with saved sacling factors}
}
\details{One scale and one center file is saved per aberration type}
\value{a list with a key \\"center\\" list and a key \\"scale\\" list. The center and scale list keys match the aberration.matrices list keys}
\author{Natalie Fox}
\examples{

example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

scaling.factors <- read.scaling.factors(example.molecular.data.dir,c('cna','csnvs','methy'));

}
