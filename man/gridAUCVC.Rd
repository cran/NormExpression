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
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{dataType}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{HG7}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{ERCC}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TN}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TC}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{CR}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{NR}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{DESeq}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{UQ}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TMM}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{TU}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{GAPDH}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{nonzeroRatios}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{cvNorm}{
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
