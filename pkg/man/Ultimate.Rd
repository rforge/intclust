\name{Ultimate}
\alias{Ultimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that performs any aggregated data function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ultimate(List, distmeasure, NN = 20, alpha = 0.5, T = 20, t = 10, 
r = NULL, nrclusters = 7, nrclusterssep = c(7, 7), nrclustersseq 
= NULL, weight = NULL, Clustweight = 0.5, clust = "agnes", linkage
 = "ward", gap = FALSE, maxK = 50, IntClust = c("ADC", "ADECa", "ADECb",
  "ADECc", "WonM", "CECa", "CECb", "CECc", "WeightedClust", 
  "WeightedSim", "SNFa", "SNFb", "SNFc"), fusionsLog = TRUE, WeightClust
   = TRUE, PlotCompare = FALSE, cols = Colors2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{NN}{
%%     ~~Describe \code{NN} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{nrclusterssep}{
%%     ~~Describe \code{nrclusterssep} here~~
}
  \item{nrclustersseq}{
%%     ~~Describe \code{nrclustersseq} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{Clustweight}{
%%     ~~Describe \code{Clustweight} here~~
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
  \item{IntClust}{
%%     ~~Describe \code{IntClust} here~~
}
  \item{fusionsLog}{
%%     ~~Describe \code{fusionsLog} here~~
}
  \item{WeightClust}{
%%     ~~Describe \code{WeightClust} here~~
}
  \item{PlotCompare}{
%%     ~~Describe \code{PlotCompare} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
