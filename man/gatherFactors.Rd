\name{gatherFactors}
\alias{gatherFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gatherFactors}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gatherFactors(data,
methods = c("HG7", "ERCC", "TN", "TC", "CR", "NR", "DESeq", "UQ", "TMM", "TU"),
HG7.size = NULL, ERCC.size = NULL, TN.size = NULL, TC.size = NULL,
CR.size = NULL, NR.size = NULL, pre_ratio = 0.5,
lower_trim = 0.05, upper_trim = 0.65, min_ubq = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{methods}{
%%     ~~Describe \code{methods} here~~
}
  \item{HG7.size}{
%%     ~~Describe \code{HG7.size} here~~
}
  \item{ERCC.size}{
%%     ~~Describe \code{ERCC.size} here~~
}
  \item{TN.size}{
%%     ~~Describe \code{TN.size} here~~
}
  \item{TC.size}{
%%     ~~Describe \code{TC.size} here~~
}
  \item{CR.size}{
%%     ~~Describe \code{CR.size} here~~
}
  \item{NR.size}{
%%     ~~Describe \code{NR.size} here~~
}
  \item{pre_ratio}{
%%     ~~Describe \code{pre_ratio} here~~
}
  \item{lower_trim}{
%%     ~~Describe \code{lower_trim} here~~
}
  \item{upper_trim}{
%%     ~~Describe \code{upper_trim} here~~
}
  \item{min_ubq}{
%%     ~~Describe \code{min_ubq} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, methods = c("HG7", "ERCC", "TN", "TC", "CR",
    "NR", "DESeq", "UQ", "TMM", "TU"), HG7.size = NULL, ERCC.size = NULL,
    TN.size = NULL, TC.size = NULL, CR.size = NULL, NR.size = NULL,
    pre_ratio = 0.5, lower_trim = 0.05, upper_trim = 0.65, min_ubq = 100)
{
    method1 <- as.list(methods)
    numMethod <- length(method1)
    method_range <- seq(1, numMethod, 1)
    for (i in method_range) {
        if (method1[[i]] == "HG7" || method1[[i]] == "ERCC" ||
            method1[[i]] == "TN" || method1[[i]] == "TC" || method1[[i]] ==
            "CR" || method1[[i]] == "NR") {
            size.name <- paste(method1[[i]], ".size", sep = "")
            out.name1 <- paste(method1[[i]], ".factors", sep = "")
            if (is.null(size.name)) {
                stop("Please provide", size.name, "!")
            }
            else {
                assign(out.name1, getFactors(data, method = "sizefactor",
                  lib.size = get(size.name)))
            }
        }
        if (method1[[i]] == "DESeq" || method1[[i]] == "RLE" ||
            method1[[i]] == "UQ" || method1[[i]] == "TMM") {
            out.name2 <- paste(method1[[i]], ".factors", sep = "")
            assign(out.name2, getFactors(data, method = method1[[i]]))
        }
        if (method1[[i]] == "TU") {
            TU.factors <- getFactors(data, method = "TU", pre_ratio = pre_ratio,
                lower_trim = lower_trim, upper_trim = upper_trim,
                min_ubq = min_ubq)
        }
    }
    factors.list <- NULL
    for (m in methods) {
        m.factors <- paste(m, ".factors", sep = "")
        factors.list <- c(factors.list, m.factors)
    }
    factors.result <- NULL
    for (i in method_range) {
        factors.result <- cbind(factors.result, get(factors.list[i]))
    }
    colnames(factors.result) <- methods
    return(factors.result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
