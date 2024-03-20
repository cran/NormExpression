\name{getCorMedians}
\alias{getCorMedians}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getCorMedians}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
getCorMedians(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data)
{
    if (!is.data.frame(data))
        data <- data.frame(data)
    if (is.factor(data$Value))
        data$Value <- as.numeric(as.character(data$Value))
    sorted_result <- sort(tapply(data$Value, data$Methods, median),
        decreasing = FALSE)
    return(sorted_result)
  }
}
