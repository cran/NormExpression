\name{findGenes}
\alias{findGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{findGenes}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
findGenes(g, qlower = NULL, qupper = NULL, pre_ratio = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{qlower}{
%%     ~~Describe \code{qlower} here~~
}
  \item{qupper}{
%%     ~~Describe \code{qupper} here~~
}
  \item{pre_ratio}{
%%     ~~Describe \code{pre_ratio} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (g, qlower = NULL, qupper = NULL, pre_ratio = NULL)
{
    gene_name <- rownames(g)
    g <- unlist(g)
    seen <- which(g >= qlower & g <= qupper)
    counts <- length(seen)
    if (counts >= pre_ratio * length(g)) {
        gene_name
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
