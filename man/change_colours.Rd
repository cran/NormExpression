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
%%     ~~Describe \code{p} here~~
}
  \item{palette}{
%%     ~~Describe \code{palette} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
