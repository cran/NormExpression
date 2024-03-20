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
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{pre_ratioC}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{lower_trimC}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{upper_trimC}{
    Please refer to the file /inst/doc/readme.pdf.
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
