\name{ColorsNames}
\alias{ColorsNames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that annotates colors to their names
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ColorsNames(MatrixColors, cols = my_palette2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MatrixColors}{
%%     ~~Describe \code{MatrixColors} here~~
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
function (MatrixColors, cols = my_palette2) 
{
    Names = c()
    for (i in 1:dim(MatrixColors)[2]) {
        for (j in 1:dim(MatrixColors)[1]) {
            temp = MatrixColors[j, i]
            Color = cols[temp]
            Names = c(Names, Color)
        }
    }
    return(Names)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
