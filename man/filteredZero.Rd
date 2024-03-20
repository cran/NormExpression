\name{filteredZero}
\alias{filteredZero}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{filteredZero}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
filteredZero(data, nonzeroRatio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{nonzeroRatio}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, nonzeroRatio)
{
    nozeroCount <- apply(data, 1, function(x) length(which(x !=
        0)))
    geneIndex <- which(nozeroCount >= ncol(data) * nonzeroRatio)
    return(geneIndex)
  }
}
