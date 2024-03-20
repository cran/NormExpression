\name{estimateSizeFactorsForMatrix}
\alias{estimateSizeFactorsForMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{estimateSizeFactorsForMatrix}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
estimateSizeFactorsForMatrix(data, p = p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{p}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, p = p)
{
    loggeomeans <- rowMeans(.log(data), na.rm = TRUE)
    apply(data, 2, function(cnts) exp(quantile(.log(cnts) - loggeomeans,
        na.rm = TRUE, p = p)))
  }
}
