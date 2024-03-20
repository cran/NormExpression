\name{CV2AUCVC}
\alias{CV2AUCVC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{CV2AUCVC}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
CV2AUCVC(data, cvResolution = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{cvResolution}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, cvResolution = 0.005)
{
    cv_cutoff <- NULL
    uniform_genes_counts <- NULL
    for (i in seq(0, 1, cvResolution)) {
        cv_cutoff <- c(cv_cutoff, i)
        gene_number <- length(which(data <= i))
        uniform_genes_counts <- c(uniform_genes_counts, gene_number)
    }
    getArea(cv_cutoff, uniform_genes_counts)
  }
}
