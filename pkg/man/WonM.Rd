\name{WonM}
\alias{WonM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighting on Membership
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WonM(List, distmeasure = c("tanimoto", "tanimoto"), nrclusters = seq(5, 25, 1), clust = "agnes", linkage = "ward")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{clust}{
%%     ~~Describe \code{clust} here~~
}
  \item{linkage}{
%%     ~~Describe \code{linkage} here~~
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
function (List, distmeasure = c("tanimoto", "tanimoto"), nrclusters = seq(5, 
    25, 1), clust = "agnes", linkage = "ward") 
{
    Distance = function(data, distmeasure) {
        data <- data + 0
        tanimoto = function(m) {
            S = matrix(0, nrow = dim(m)[1], ncol = dim(m)[1])
            for (i in 1:dim(m)[1]) {
                for (j in 1:i) {
                  N.A = sum(m[i, ])
                  N.B = sum(m[j, ])
                  N.C = sum(m[i, (m[i, ] == m[j, ])])
                  if (N.A == 0 & N.B == 0) {
                    coef = 1
                  }
                  else {
                    coef = N.C/(N.A + N.B - N.C)
                  }
                  S[i, j] = coef
                  S[j, i] = coef
                }
            }
            D = 1 - S
            return(D)
        }
        if (distmeasure == "jaccard") {
            dist = dist.binary(data, method = 1)
            dist = as.matrix(dist)
        }
        else if (distmeasure == "tanimoto") {
            dist = tanimoto(data)
            dist = as.matrix(dist)
            rownames(dist) <- rownames(data)
        }
        else if (distmeasure == "euclidean") {
            dist = daisy(data, metric = "euclidean")
            dist = as.matrix(dist)
        }
        else {
            stop("Incorrect choice of distmeasure. Must be one of: tanimoto, jaccard or euclidean.")
        }
        return(dist)
    }
    DistM = vector("list", length(List))
    for (i in 1:length(List)) {
        DistM[[i]] = Distance(List[[i]], distmeasure[i])
        rownames(DistM[[i]]) = rownames(List[[i]])
        colnames(DistM[[i]]) = rownames(List[[i]])
    }
    HClustering = list()
    if (clust == "agnes" & linkage == "ward") {
        for (i in 1:length(List)) {
            HClustering[[i]] = agnes(DistM[[i]], diss = TRUE, 
                method = linkage)
        }
    }
    Consensus = list()
    for (i in 1:length(List)) {
        Consensus[[i]] = matrix(0, dim(List[[i]])[1], dim(List[[i]])[1])
        rownames(Consensus[[i]]) = rownames(List[[i]])
        colnames(Consensus[[i]]) = rownames(List[[i]])
    }
    for (k in nrclusters) {
        for (j in 1:length(List)) {
            Temp = cutree(HClustering[[j]], k)
            CM = matrix(0, dim(List[[j]])[1], dim(List[[j]])[1])
            for (l in 1:length(Temp)) {
                label = Temp[l]
                sameclust = which(Temp == label)
                CM[l, sameclust] = 1
            }
            Consensus[[j]] = Consensus[[j]] + CM
        }
    }
    OverallConsensus = matrix(0, dim(List[[1]])[1], dim(List[[1]])[1])
    for (i in 1:length(Consensus)) {
        OverallConsensus = OverallConsensus + Consensus[[i]]
    }
    OverallClustering = agnes(OverallConsensus, diss = FALSE, 
        method = linkage)
    out = list(DistanceM = DistM, ClustSep = HClustering, Consensus = OverallConsensus, 
        Clust = OverallClustering)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
