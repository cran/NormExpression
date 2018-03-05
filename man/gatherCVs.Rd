\name{gatherCVs}
\alias{gatherCVs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gatherCVs}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gatherCVs(data,nonzeroRatio,HG7,ERCC,TN,TC,CR,NR,
DESeq,UQ,TMM,TU,GAPDH,cvNorm,cvResolution)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{nonzeroRatio}{
%%     ~~Describe \code{nonzeroRatio} here~~
}
  \item{HG7}{
%%     ~~Describe \code{HG7} here~~
}
  \item{ERCC}{
%%     ~~Describe \code{ERCC} here~~
}
  \item{TN}{
%%     ~~Describe \code{TN} here~~
}
  \item{TC}{
%%     ~~Describe \code{TC} here~~
}
  \item{CR}{
%%     ~~Describe \code{CR} here~~
}
  \item{NR}{
%%     ~~Describe \code{NR} here~~
}
  \item{DESeq}{
%%     ~~Describe \code{DESeq} here~~
}
  \item{UQ}{
%%     ~~Describe \code{UQ} here~~
}
  \item{TMM}{
%%     ~~Describe \code{TMM} here~~
}
  \item{TU}{
%%     ~~Describe \code{TU} here~~
}
  \item{GAPDH}{
%%     ~~Describe \code{GAPDH} here~~
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
function (data, nonzeroRatio = NULL, HG7 = NULL, ERCC = NULL,
    TN = NULL, TC = NULL, CR = NULL, NR = NULL, DESeq = NULL,
    UQ = NULL, TMM = NULL, TU = NULL, GAPDH = NULL, cvNorm = TRUE,
    cvResolution = 0.005)
{
    if (is.null(nonzeroRatio)) {
        stop("Please provide nonzeroRatio!")
    }
    methodsList <- list(HG7 = HG7, ERCC = ERCC, TN = TN, TC = TC,
        CR = CR, NR = NR, DESeq = DESeq, UQ = UQ, TMM = TMM,
        TU = TU, GAPDH = GAPDH)
    specifiedMethods <- methodsList[!unlist(lapply(methodsList,
        is.null))]
    numMethod <- length(specifiedMethods)
    method_range_tmp <- seq(1, numMethod, 1)
    cv_range_tmp <- seq(0, 1, cvResolution)
    method_range_times <- length(cv_range_tmp)
    cv_range_times <- length(method_range_tmp)
    method_range <- rep(method_range_tmp, each = round(method_range_times))
    cv_range <- rep(cv_range_tmp, times = round(cv_range_times))
    nozeroIndex <- filteredZero(data, nonzeroRatio = nonzeroRatio)
    for (j in method_range_tmp) {
        norm.matrix <- getNormMatrix(data, specifiedMethods[[j]])
        dataUse2CV <- norm.matrix[nozeroIndex, ]
        cv.result <- getCV(dataUse2CV, cvNorm = cvNorm)
        assign(paste(names(specifiedMethods)[j], ".cv", sep = ""),
            cv.result)
    }
    cv_uniform <- NULL
    cv_uniform_all <- mapply(function(i, j) {
        cv.result <- paste(names(specifiedMethods)[j], ".cv",
            sep = "")
        gene_number <- length(which(get(cv.result) <= i))
        cv_uniform_row <- c(i, gene_number, names(specifiedMethods)[j])
        rbind(cv_uniform, cv_uniform_row)
    }, cv_range, method_range)
    cv_uniform_all <- t(cv_uniform_all)
    colnames(cv_uniform_all) <- c("Cutoff", "Counts", "Methods")
    return(cv_uniform_all)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
