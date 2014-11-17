\name{BoxPlotDist}
\alias{BoxPlotDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Box Plot of Distances
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
BoxPlotDist(Data1, Data2, lab1, lab2, limits1, limits2, type = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data1}{
%%     ~~Describe \code{Data1} here~~
}
  \item{Data2}{
%%     ~~Describe \code{Data2} here~~
}
  \item{lab1}{
%%     ~~Describe \code{lab1} here~~
}
  \item{lab2}{
%%     ~~Describe \code{lab2} here~~
}
  \item{limits1}{
%%     ~~Describe \code{limits1} here~~
}
  \item{limits2}{
%%     ~~Describe \code{limits2} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
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
function (Data1, Data2, lab1, lab2, limits1, limits2, type = 1) 
{
    Dist1 <- Data1$DistM
    Dist2 <- Data2$DistM
    Dist1lower <- Dist1[lower.tri(Dist1)]
    Dist2lower <- Dist2[lower.tri(Dist2)]
    Categorize <- function(Distlower, limits) {
        Cat = c(rep(0, length(Distlower)))
        for (j in 1:(length(limits) + 1)) {
            if (j == 1) {
                Cat[Distlower <= limits[j]] = j
            }
            else if (j <= length(limits)) {
                Cat[Distlower > limits[j - 1] & Distlower <= 
                  limits[j]] = j
            }
            else {
                Cat[Distlower > limits[j - 1]] = j
            }
        }
        Cat <- factor(Cat)
        return(Cat)
    }
    Dist1cat <- Categorize(Dist1lower, limits1)
    Dist2cat <- Categorize(Dist2lower, limits2)
    dataBox1 <- data.frame(D1 = Dist1lower, C2 = Dist2cat)
    p1 <- ggplot(dataBox1, aes(factor(C2), D1))
    p1 <- p1 + geom_boxplot(outlier.shape = NA) + geom_point(color = "blue", 
        size = 2, shape = 19, position = "jitter", cex = 1.5) + 
        xlab(lab2) + ylab(lab1)
    dataBox2 <- data.frame(D2 = Dist2lower, C1 = Dist1cat)
    p2 <- ggplot(dataBox2, aes(factor(C1), D2))
    p2 <- p2 + geom_boxplot(outlier.shape = NA) + geom_point(color = "blue", 
        size = 2, shape = 19, position = "jitter", cex = 1.5) + 
        xlab(lab1) + ylab(lab2)
    if (type == 3) {
        dev.new(width = 14, height = 7)
        grid.arrange(p1, p2, ncol = 2, nrow = 1)
    }
    else if (type == 1) {
        dev.new()
        print(p1)
    }
    else if (type == 2) {
        dev.new()
        print(p2)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
