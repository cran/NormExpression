\name{getNormMatrix}
\alias{getNormMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getNormMatrix}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
getNormMatrix(data, norm.factors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{norm.factors}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, norm.factors)
{
    data * matrix(rep(norm.factors, dim(data)[1]), nrow = dim(data)[1],
        ncol = length(norm.factors), byrow = T)
  }
}
