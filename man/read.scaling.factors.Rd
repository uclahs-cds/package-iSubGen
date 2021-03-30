\name{read.scaling.factors}
\alias{read.scaling.factors}
\title{Read scaling factors from file}
\description{Read scaling factors from file}
\usage{
read.scaling.factors(scaling.factor.files.dir,data.types);
}
\arguments{
  \item{scaling.factor.files.dir}{the directory where the files were saved}
  \item{data.types}{a vector of the data types with saved scaling factors}
}
\details{One scale and one center file is saved per data type}
\value{a list with a key \\"center\\" list and a key \\"scale\\" list. The center and scale list keys match the data.matrices list keys}
\author{Natalie Fox}
\examples{

# Get the path for the scaling provided in this R package
example.molecular.data.dir <- paste0(path.package('iSubGen'),'/exdata/');

# Example #1: reading scaling factors for a single data type
scaling.factors <- read.scaling.factors(example.molecular.data.dir, 'cna');

# Example #2: reading scaling factors for multiple data types
scaling.factors <- read.scaling.factors(example.molecular.data.dir, c('cna','snv','methy'));
}
