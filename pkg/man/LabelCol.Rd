\name{labelCol}
\alias{labelCol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for colors in LabelPlot
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
labelCol(x, Sel)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{Sel}{
%%     ~~Describe \code{Sel} here~~
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
function (x, Sel) 
{
    if (is.leaf(x)) {
        label <- attr(x, "label")
        attr(x, "nodePar") <- list(pch = NA, lab.col = ifelse(label \%in\% 
            Sel, "red", "black"), lab.cex = 0.9, font = 2)
        attr(x, "edgePar") <- list(lwd = 2, col = ifelse(label \%in\% 
            Sel, "red", "black"))
    }
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
