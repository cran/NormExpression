\name{gridAUCVC4Matrices}
\alias{gridAUCVC4Matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gridAUCVC4Matrices}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gridAUCVC4Matrices(..., nonzeroRatios = NULL, cvNorm = TRUE, cvResolution = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{nonzeroRatios}{
%%     ~~Describe \code{nonzeroRatios} here~~
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
function (..., nonzeroRatios = NULL, cvNorm = TRUE, cvResolution = 0.005)
{
    if (is.null(nonzeroRatios)) {
        stop("Please provide nonzeroRatios!")
    }
    matrices <- list(...)
    numMethod <- length(matrices)
    grid_result <- NULL
    for (i in nonzeroRatios) {
        result.sorted <- getAUCVCs(..., nonzeroRatio = i, cvNorm = cvNorm,
            cvResolution = cvResolution)
        grid_record <- c(i, result.sorted)
        names(grid_record)[1] <- "NonzeroRatio"
        grid_result <- c(grid_result, names(grid_record), grid_record)
    }
    grid_result2 <- matrix(grid_result, ncol = numMethod + 1,
        byrow = TRUE)
    return(grid_result2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
