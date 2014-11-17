\name{ClusterPlot}
\alias{ClusterPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the dendrogram but colored by a result of choice
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ClusterPlot(Data1, Data2, nrclusters = 7, cols = my_palette2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data1}{
%%     ~~Describe \code{Data1} here~~
}
  \item{Data2}{
%%     ~~Describe \code{Data2} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (Data1, Data2, nrclusters = 7, cols = my_palette2, ...) 
{
    x = Data1$Clust
    Data = Data2$Clust
    d_temp <- dendrapply(as.dendrogram(x, hang = 0.02), ClusterCols, 
        Data, nrclusters, cols)
    plot(d_temp, nodePar = list(pch = NA), edgePar = list(lwd = 2), 
        ylab = "Height", font.axis = 2, font.lab = 2, font = 2, 
        ...)
    axis(side = 2, lwd = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
