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
%%     ~~Describe \code{data} here~~
}
  \item{norm.factors}{
%%     ~~Describe \code{norm.factors} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
