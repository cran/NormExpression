\name{calcFactorWeighted}
\alias{calcFactorWeighted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{calcFactorWeighted}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
calcFactorWeighted(obs, ref, libsize.obs, libsize.ref, logratioTrim,
sumTrim, doWeighting, Acutoff)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obs}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{ref}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{libsize.obs}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{libsize.ref}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{logratioTrim}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{sumTrim}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{doWeighting}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{Acutoff}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
    sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10)
{
    if (all(obs == ref))
        return(1)
    obs <- as.numeric(obs)
    ref <- as.numeric(ref)
    if (is.null(libsize.obs))
        nO <- sum(obs)
    else nO <- libsize.obs
    if (is.null(libsize.ref))
        nR <- sum(ref)
    else nR <- libsize.ref
    logR <- log2((obs/nO)/(ref/nR))
    absE <- (log2(obs/nO) + log2(ref/nR))/2
    v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
    fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
    logR <- logR[fin]
    absE <- absE[fin]
    v <- v[fin]
    n <- length(logR)
    loL <- floor(n * logratioTrim) + 1
    hiL <- n + 1 - loL
    loS <- floor(n * sumTrim) + 1
    hiS <- n + 1 - loS
    keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >=
        loS & rank(absE) <= hiS)
    if (doWeighting) {
        2^(sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],
            na.rm = TRUE))
    }
    else {
        2^(mean(logR[keep], na.rm = TRUE))
    }
  }
}
