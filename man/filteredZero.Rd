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
%%     ~~Describe \code{data} here~~
}
  \item{nonzeroRatio}{
%%     ~~Describe \code{nonzeroRatio} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
