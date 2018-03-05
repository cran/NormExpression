\name{gatherCors4Matrices}
\alias{gatherCors4Matrices}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{gatherCors4Matrices}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
gatherCors4Matrices(..., raw_matrix, cor_method = c("spearman", "pearson", "kendall"),
pre_ratio = 0.5, lower_trim = 0.05, upper_trim = 0.65, rounds = 1e+06)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{raw_matrix}{
%%     ~~Describe \code{raw_matrix} here~~
}
  \item{cor_method}{
%%     ~~Describe \code{cor_method} here~~
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
  \item{rounds}{
%%     ~~Describe \code{rounds} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (..., raw_matrix, cor_method = c("spearman", "pearson",
    "kendall"), pre_ratio = 0.5, lower_trim = 0.05, upper_trim = 0.65,
    rounds = 1e+06)
{
    matrices <- list(...)
    numMethod <- length(matrices)
    method_range <- seq(1, numMethod, 1)
    ubq_genes <- identifyUbq(raw_matrix, pre_ratio = pre_ratio,
        lower_trim = lower_trim, upper_trim = upper_trim, min_ubq = 100)
    cor_value_method <- NULL
    for (j in method_range) {
        dataUse2Cor <- matrices[[j]][ubq_genes, ]
        cor.result <- getCor(dataUse2Cor, method = cor_method,
            rounds = rounds)
        cor_vm <- cbind(cor.result, rep(names(matrices)[j], times = round(rounds)))
        cor_value_method <- rbind(cor_value_method, cor_vm)
    }
    colnames(cor_value_method) <- c("Value", "Methods")
    return(cor_value_method)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
