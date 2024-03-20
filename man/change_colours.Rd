\name{change_colours}
\alias{change_colours}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{change_colours}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
change_colours(p, palette, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{palette}{
    Please refer to the file /inst/doc/readme.pdf.
}
  \item{type}{
    Please refer to the file /inst/doc/readme.pdf.
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (p, palette, type)
{
    n <- nlevels(p$data[[deparse(p$mapping$group)]])
    tryCatch(as.character(palette), error = function(e) stop("be vector",call. = FALSE))
    if (n > length(palette))
        stop("Not enough colours in palette.")
    if (missing(type))
        type <- grep("colour|fill", names(p$layers[[1]]$mapping),
            value = TRUE)[1]
    pal <- function(n) palette[seq_len(n)]
    p + discrete_scale(type, "foo", pal)
  }
}
