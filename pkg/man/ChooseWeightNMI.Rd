\name{ChooseWeightNMI}
\alias{ChooseWeightNMI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the weight for weighted similarity clustering with NMI
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChooseWeightNMI(List, type = c("data", "clusters"), w = 
seq(0, 1, by = 0.01), nrclusters = NULL, distmeasure = 
c("euclidean", "tanimoto"), clust = "agnes", linkage = "ward", 
gap = FALSE, maxK = 50, names = c("B", "FP"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{w}{
%%     ~~Describe \code{w} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{clust}{
%%     ~~Describe \code{clust} here~~
}
  \item{linkage}{
%%     ~~Describe \code{linkage} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
  \item{maxK}{
%%     ~~Describe \code{maxK} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
