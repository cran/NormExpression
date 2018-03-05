\name{identifyUbqRepeat}
\alias{identifyUbqRepeat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{identifyUbqRepeat}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
identifyUbqRepeat(data, pre_ratioC = NULL, lower_trimC = NULL, upper_trimC = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{pre_ratioC}{
%%     ~~Describe \code{pre_ratioC} here~~
}
  \item{lower_trimC}{
%%     ~~Describe \code{lower_trimC} here~~
}
  \item{upper_trimC}{
%%     ~~Describe \code{upper_trimC} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, pre_ratioC = NULL, lower_trimC = NULL, upper_trimC = NULL)
{
    qlower <- apply(data, 2, function(x) quantile(x[x != 0],
        p = lower_trimC))
    qupper <- apply(data, 2, function(x) quantile(x[x != 0],
        p = upper_trimC))
    ubq_genes <- NULL
    for (i in 1:nrow(data)) {
        genes_finded <- findGenes(data[i, ], qlower = qlower,
            qupper = qupper, pre_ratio = pre_ratioC)
        ubq_genes <- c(ubq_genes, genes_finded)
    }
    return(ubq_genes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
