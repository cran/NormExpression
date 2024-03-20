\name{calcFactorRLE}
\alias{calcFactorRLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calcFactorRLE}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
calcFactorRLE(data, p = p)
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
    gm <- exp(rowMeans(.log(data), na.rm = TRUE))
    apply(data, 2, function(u) quantile((u/gm)[u != 0], na.rm = TRUE,
        p = p))
  }
}
