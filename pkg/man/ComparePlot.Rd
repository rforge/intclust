\name{ComparePlot}
\alias{ComparePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare the result of multiple aggregated data outcomes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComparePlot(List, nrclusters = 7, cols = my_palette2, fusionsLog = FALSE, WeightClust = FALSE, names = NULL, reverse = FALSE, margins = c(8.1, 3.1, 3.1, 4.1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{fusionsLog}{
%%     ~~Describe \code{fusionsLog} here~~
}
  \item{WeightClust}{
%%     ~~Describe \code{WeightClust} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{reverse}{
%%     ~~Describe \code{reverse} here~~
}
  \item{margins}{
%%     ~~Describe \code{margins} here~~
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
function (List, nrclusters = 7, cols = my_palette2, fusionsLog = FALSE, 
    WeightClust = FALSE, names = NULL, reverse = FALSE, margins = c(8.1, 
        3.1, 3.1, 4.1), ...) 
{
    if (length(List) == 1 & (class(List[[1]]) == "Weighted" | 
        class(List[[1]]) == "CEC" | class(List[[1]]) == "WeightedSim") & 
        reverse == TRUE) {
        List[[1]]$Results = rev(List[[1]]$Results)
    }
    MatrixColors = MatrixFunction(List, nrclusters, fusionsLog, 
        WeightClust, names)
    Names = ColorsNames(MatrixColors, cols)
    nobs = dim(MatrixColors)[2]
    nmethods = dim(MatrixColors)[1]
    if (is.null(names)) {
        for (j in 1:nmethods) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    similar = round(SimilarityMeasure(MatrixColors), 2)
    par(mar = margins)
    color2D.matplot(MatrixColors, cellcolors = Names, show.values = FALSE, 
        axes = FALSE, xlab = "", ylab = "", ...)
    axis(1, at = seq(0.5, (nobs - 0.5)), labels = colnames(MatrixColors), 
        las = 2, cex.axis = 0.7)
    axis(2, at = seq(0.5, (nmethods - 0.5)), labels = rev(names), 
        cex.axis = 0.65, las = 2)
    axis(4, at = seq(0.5, (nmethods - 0.5)), labels = rev(similar), 
        cex.axis = 0.65, las = 2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
