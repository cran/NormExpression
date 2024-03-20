\name{calcFactorUpperquartile}
\alias{calcFactorUpperquartile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calcFactorUpperquartile}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
calcFactorUpperquartile(data, lib.size, p = p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{lib.size}{
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
function (data, lib.size, p = p)
{
    y <- t(t(data)/lib.size)
    f <- apply(y, 2, function(x) quantile(x[x != 0], p = p))
  }
}
