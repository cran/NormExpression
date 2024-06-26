\name{getCor}
\alias{getCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCor}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
getCor(data, method = c("spearman", "pearson", "kendall"), rounds = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{method}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{rounds}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, method = c("spearman", "pearson", "kendall"),
    rounds = 1e+06)
{
    sp_result <- NULL
    method <- match.arg(method)
    for (i in 1:rounds) {
        rg1 <- sample(1:nrow(data), size = 1)
        rg2 <- sample(1:nrow(data), size = 1)
        while (rg1 == rg2) {
            rg2 <- sample(1:nrow(data), size = 1)
        }
        gene1 <- unlist(data[rg1, ])
        gene2 <- unlist(data[rg2, ])
        sp_value <- cor(gene1, gene2, method = method)
        sp_result <- c(sp_result, sp_value)
    }
    return(sp_result)
  }
}
