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
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{methods}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{HG7.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{ERCC.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TN.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TC.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{CR.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{NR.size}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{pre_ratio}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{lower_trim}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{upper_trim}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{min_ubq}{
    Please refer to the file /inst/doc/readme.pdf.
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
