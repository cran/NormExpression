\name{plotCors}
\alias{plotCors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plotCors}
\description{
Please refer to the file /inst/doc/readme.pdf.
}
\usage{
plotCors(data, methods = c("None", "HG7", "ERCC", "TN", "TC", "CR", "NR", "DESeq",
"UQ", "TMM", "TU"), legend.position = c(0.15, 0.56))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
%%     ~~Describe \code{data} here~~
}
  \item{methods}{
%%     ~~Describe \code{methods} here~~
}
  \item{legend.position}{
%%     ~~Describe \code{legend.position} here~~
}
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data, methods = c("None", "HG7", "ERCC", "TN", "TC",
    "CR", "NR", "DESeq", "UQ", "TMM", "TU"), legend.position = c(0.15,
    0.56))
{
    if (!is.data.frame(data))
        data <- data.frame(data)
    if (is.factor(data$Value))
        data$Value <- as.numeric(as.character(data$Value))
    data$Methods <- factor(data$Methods, levels = methods, labels = methods)
    change_colours(ggplot(data = data, aes(x = Value, y = ..count../sum(..count..))) +
        geom_freqpoly(aes(group = Methods, color = Methods),
            size = 3, bins = 50) + xlab("Spearman correlation") +
        ylab("Fraction of gene pairs") + theme_bw() + theme(panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 48), axis.title.y = element_text(size = 48),
        axis.text.x = element_text(size = 38), axis.text.y = element_text(size = 38),
        legend.text = element_text(size = 39), legend.title = element_text(size = 43),
        legend.position = legend.position, legend.background = element_blank(),
        legend.key = element_blank(), legend.key.height = unit(1.8,
            "cm"), plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")) +
        scale_x_continuous(expand = c(0.01, 0.01), breaks = round(seq(-1,
            1, 0.25), 2)) + scale_y_continuous(expand = c(0.01,
        0)) + guides(color = guide_legend(title = NULL)), c("olivedrab",
        "blue", "red", "violet", "orange", "yellow", "magenta",
        "peru", "black", "maroon", "lightblue", "darkslateblue",
        "seashell4", "tan2", "darkgreen", "springgreen"))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
