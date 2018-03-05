\name{getAUCVCs}
\alias{getAUCVCs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{getAUCVCs}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
getAUCVCs(..., nonzeroRatio = NULL, cvNorm = TRUE, cvResolution = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{nonzeroRatio}{
%%     ~~Describe \code{nonzeroRatio} here~~
}
  \item{cvNorm}{
%%     ~~Describe \code{cvNorm} here~~
}
  \item{cvResolution}{
%%     ~~Describe \code{cvResolution} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., nonzeroRatio = NULL, cvNorm = TRUE, cvResolution = 0.005)
{
    matrices <- list(...)
    numMethod <- length(matrices)
    method_range <- seq(1, numMethod, 1)
    result <- NULL
    for (i in method_range) {
        AUCVC.result <- getAUCVC(matrices[[i]], nonzeroRatio = nonzeroRatio,
            cvNorm = cvNorm, cvResolution = cvResolution)
        result <- c(result, AUCVC.result)
        names(result)[i] <- names(matrices)[i]
    }
    sorted_AUCVCs <- sort(result, decreasing = TRUE)
    return(sorted_AUCVCs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
