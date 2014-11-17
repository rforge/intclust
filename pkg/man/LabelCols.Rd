\name{LabelCols}
\alias{LabelCols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for colors in LabelPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LabelCols(x, Sel1, Sel2 = NULL, col1, col2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
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
function (x, Sel1, Sel2 = NULL, col1, col2 = NULL) 
{
    colfunc = function(x, Sel1, Sel2, col1, col2) {
        if (x \%in\% Sel1) {
            return(col1)
        }
        else if (x \%in\% Sel2) {
            return(col2)
        }
        else {
            return("black")
        }
    }
    if (is.leaf(x)) {
        label <- attr(x, "label")
        attr(x, "nodePar") <- list(pch = NA, lab.col = colfunc(label, 
            Sel1, Sel2, col1, col2), lab.cex = 0.9, font = 2)
        attr(x, "edgePar") <- list(lwd = 2, col = colfunc(label, 
            Sel1, Sel2, col1, col2))
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
