\name{gridAUCVC}
\alias{gridAUCVC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gridAUCVC}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gridAUCVC(data, dataType = c("bk", "sc"), HG7 = NULL, ERCC = NULL, TN = NULL,
TC = NULL, CR = NULL, NR = NULL, DESeq = NULL, UQ = NULL, TMM = NULL, TU = 0,
GAPDH = NULL, nonzeroRatios = c(0.7, 0.8, 0.9, 1), cvNorm = TRUE, cvResolution = 0.005)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{dataType}{
%%     ~~Describe \code{dataType} here~~
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
function (data, dataType = c("bk", "sc"), HG7 = NULL, ERCC = NULL,
    TN = NULL, TC = NULL, CR = NULL, NR = NULL, DESeq = NULL,
    UQ = NULL, TMM = NULL, TU = 0, GAPDH = NULL, nonzeroRatios = c(0.7,
        0.8, 0.9, 1), cvNorm = TRUE, cvResolution = 0.005)
{
    grid_result <- NULL
    if (length(TU) == 1 && TU == 1) {
        colnames_paraMatrix <- c("nonzeroRatio", "pre_ratio",
            "lower_trim", "upper_trim")
        write.table(t(as.matrix(colnames_paraMatrix)), file = "bestPara.txt",
            sep = "\t", row.names = FALSE, col.names = FALSE)
    }
    for (i in nonzeroRatios) {
        if (dataType == "sc") {
            if ((ncol(data) * i) <= 100) {
                cat("nonzeroRatio:", i, " is too small!\n")
                stop("We suggest that the minimal counts of
                nonzero samples should be greater than 100!")
            }
        }
        result <- nonzeroRatio2AUCVC(data = data, dataType = dataType,
            HG7 = HG7, ERCC = ERCC, TN = TN, TC = TC, CR = CR,
            NR = NR, DESeq = DESeq, UQ = UQ, TMM = TMM, TU = TU,
            GAPDH = GAPDH, nonzeroRatio = i, cvNorm = cvNorm,
            cvResolution = cvResolution)
        nonzeroM <- matrix(i, 1, 1, TRUE)
        colnames(nonzeroM) <- "NonzeroRatio"
        grid_record <- cbind(nonzeroM, result)
        grid_result <- rbind(grid_result, grid_record)
    }
    return(grid_result)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
