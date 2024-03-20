\name{getCV}
\alias{getCV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCV}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
getCV(data, cvNorm = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{cvNorm}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, cvNorm = TRUE)
{
    if (!is.matrix(data))
        data <- as.matrix(data)
    if (cvNorm) {
        rawCV <- apply(data, 1, function(x) {
            sd(log2(x[x != 0]))/mean(log2(x[x != 0]))
        })
        (rawCV - min(rawCV))/(max(rawCV) - min(rawCV))
    }
    else {
        apply(data, 1, function(x) {
            sd(x)/mean(x)
        })
    }
  }
}
