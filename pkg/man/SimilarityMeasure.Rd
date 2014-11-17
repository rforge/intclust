\name{SimilarityMeasure}
\alias{SimilarityMeasure}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A measure of similarity over the functions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SimilarityMeasure(List, nrclusters = 7, fusionsLog = TRUE, WeightClust = TRUE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
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
function (List, nrclusters = 7, fusionsLog = TRUE, WeightClust = TRUE, 
    names = NULL) 
{
    if (class(List) != "list") {
        MatrixColors = List
    }
    else {
        MatrixColors = MatrixFunction(List, nrclusters, fusionsLog, 
            WeightClust, names)
    }
    Similarity = c()
    for (i in 1:dim(MatrixColors)[1]) {
        Shared = 0
        for (j in 1:dim(MatrixColors)[2]) if (MatrixColors[i, 
            j] == MatrixColors[1, j]) {
            Shared = Shared + 1
        }
        Similarity = c(Similarity, Shared/ncol(MatrixColors))
    }
    return(Similarity)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
