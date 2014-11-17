\name{ClusterCols}
\alias{ClusterCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Colors for ClusterPlot function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ClusterCols(x, Data, nrclusters = 7, cols = my_palette2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, Data, nrclusters = 7, cols = my_palette2) 
{
    if (length(cols) < nrclusters) {
        stop("Not for every cluster a color is specified")
    }
    Clustdata = cutree(Data, nrclusters)
    Clustdata = Clustdata[Data$order]
    ordercolors = Clustdata
    order = seq(1, nrclusters)
    for (k in 1:length(unique(Clustdata))) {
        select = which(Clustdata == unique(Clustdata)[k])
        ordercolors[select] = order[k]
    }
    colfunc = function(x, cols) {
        indextemp = which(attr(Data$diss, "Labels") == x)
        index1 = which(Data$order == indextemp)
        index2 = ordercolors[index1]
        color = cols[index2]
        return(color)
    }
    if (is.leaf(x)) {
        label <- attr(x, "label")
        attr(x, "nodePar") <- list(pch = NA, lab.col = colfunc(label, 
            cols), lab.cex = 0.9, font = 2)
        attr(x, "edgePar") <- list(lwd = 2, col = colfunc(label, 
            cols))
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
