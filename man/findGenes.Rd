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
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{qlower}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{qupper}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{pre_ratio}{
    Please refer to the file /inst/doc/readme.pdf.
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
