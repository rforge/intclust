\name{LabelPlot}
\alias{LabelPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Color the leaves of the dendrogram
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LabelPlot(Data, Sel1, Sel2 = NULL, col1, col2 = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{Sel1}{
%%     ~~Describe \code{Sel1} here~~
}
  \item{Sel2}{
%%     ~~Describe \code{Sel2} here~~
}
  \item{col1}{
%%     ~~Describe \code{col1} here~~
}
  \item{col2}{
%%     ~~Describe \code{col2} here~~
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
function (Data, Sel1, Sel2 = NULL, col1, col2 = NULL, ...) 
{
    x = Data$Clust
    d_temp <- dendrapply(as.dendrogram(x, hang = 0.02), labelCols, 
        Sel1, Sel2, col1, col2)
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
