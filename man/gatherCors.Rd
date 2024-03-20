\name{gatherCors}
\alias{gatherCors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gatherCors}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gatherCors(data, cor_method = c("spearman", "pearson", "kendall"),
HG7 = NULL, ERCC = NULL, TN = NULL, TC = NULL, CR = NULL, NR = NULL,
DESeq = NULL, UQ = NULL, TMM = NULL, TU = NULL, GAPDH = NULL,
pre_ratio = 0.5, lower_trim = 0.05, upper_trim = 0.65, rounds = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{cor_method}{
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
  \item{pre_ratio}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{lower_trim}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{upper_trim}{
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
function (data, cor_method = c("spearman", "pearson", "kendall"),
    HG7 = NULL, ERCC = NULL, TN = NULL, TC = NULL, CR = NULL,
    NR = NULL, DESeq = NULL, UQ = NULL, TMM = NULL, TU = NULL,
    GAPDH = NULL, pre_ratio = 0.5, lower_trim = 0.05, upper_trim = 0.65,
    rounds = 1e+06)
{
    methodsList <- list(HG7 = HG7, ERCC = ERCC, TN = TN, TC = TC,
        CR = CR, NR = NR, DESeq = DESeq, UQ = UQ, TMM = TMM,
        TU = TU, GAPDH = GAPDH)
    specifiedMethods <- methodsList[!unlist(lapply(methodsList,
        is.null))]
    numMethod <- length(specifiedMethods)
    method_range <- seq(1, numMethod, 1)
    ubq_genes <- identifyUbq(data, pre_ratio = pre_ratio, lower_trim = lower_trim,
        upper_trim = upper_trim, min_ubq = 100)
    cor_value_method <- NULL
    for (j in method_range) {
        norm.matrix <- getNormMatrix(data, specifiedMethods[[j]])
        dataUse2Cor <- norm.matrix[ubq_genes, ]
        cor.result <- getCor(dataUse2Cor, method = cor_method,
            rounds = rounds)
        cor_vm <- cbind(cor.result, rep(names(specifiedMethods)[j],
            times = round(rounds)))
        cor_value_method <- rbind(cor_value_method, cor_vm)
    }
    colnames(cor_value_method) <- c("Value", "Methods")
    return(cor_value_method)
  }
}
