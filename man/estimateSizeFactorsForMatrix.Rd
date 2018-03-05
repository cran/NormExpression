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
%%     ~~Describe \code{data} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
