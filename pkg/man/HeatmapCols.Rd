\name{HeatmapCols}
\alias{HeatmapCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heatmap of two dendrograms 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HeatmapCols(Data1, Data2, names = rownames(fingerprintMat), nrclusters = 7, cols = my_palette)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data1}{
%%     ~~Describe \code{Data1} here~~
}
  \item{Data2}{
%%     ~~Describe \code{Data2} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
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
function (Data1, Data2, names = rownames(fingerprintMat), nrclusters = 7, 
    cols = my_palette) 
{
    data1 = Data1$Clust
    data2 = Data2$Clust
    DistM = distanceheatmaps(data1, data2, names, nrclusters)
    heatmap.2(DistM, Rowv = as.dendrogram(data1), Colv = as.dendrogram(data2), 
        trace = "none", col = cols, key = FALSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
