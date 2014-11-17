\name{WeightedSimClust}
\alias{WeightedSimClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted similarity clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WeightedSimClust(List, type = c("data", "clusters"), w = seq(0, 1, 0.01), clust = "agnes", linkage = "ward", distmeasure = c("euclidean", "tanimoto"), gap = FALSE, maxK = 50, nrclusters = NULL, names = c("B", "FP"), AllClusters = FALSE)
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
  \item{clust}{
%%     ~~Describe \code{clust} here~~
}
  \item{linkage}{
%%     ~~Describe \code{linkage} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
  \item{maxK}{
%%     ~~Describe \code{maxK} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{AllClusters}{
%%     ~~Describe \code{AllClusters} here~~
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
function (List, type = c("data", "clusters"), w = seq(0, 1, 0.01), 
    clust = "agnes", linkage = "ward", distmeasure = c("euclidean", 
        "tanimoto"), gap = FALSE, maxK = 50, nrclusters = NULL, 
    names = c("B", "FP"), AllClusters = FALSE) 
{
    if (length(w) == 1) {
        if (type == "data") {
            for (a in 1:length(distmeasure)) {
                if (distmeasure[[a]] == "euclidean") {
                  stand <- function(c) {
                    minc = min(c)
                    maxc = max(c)
                    c1 = (c - minc)/(maxc - minc)
                    return(c1)
                  }
                  List[[a]] = apply(List[[a]], 2, stand)
                }
            }
            Clust1 <- Cluster(List[[1]], distmeasure = distmeasure[1], 
                clust, linkage, gap, maxK)
            Clust2 <- Cluster(List[[2]], distmeasure = distmeasure[2], 
                clust, linkage, gap, maxK)
            Dist1 = Clust1$DistM
            Dist2 = Clust2$DistM
        }
        else {
            Dist1 = List[[1]]$DistM
            Dist2 = List[[2]]$DistM
            if ((all(range(Dist1) == c(0, 1)) & (range(Dist2)[1] != 
                0 & range(Dist2)[2] != 1)) | ((range(Dist1)[1] != 
                0 & range(Dist1)[2] != 1) & all(range(Dist2) == 
                c(0, 1)))) {
                print("warning: the distance matrices do not have the same range. It might be that that one data set is bainary and the other is continuous. Standardzie the variables.")
            }
        }
        Weight = w
    }
    else {
        temp = ChooseWeight(List, type, w, nrclusters, distmeasure, 
            clust, linkage, gap, maxK, names)
        Dist1 = temp[[1]]$DistM
        Dist2 = temp[[2]]$DistM
        Weight = temp[[4]]
    }
    DistW = Weight * Dist1 + (1 - Weight) * Dist2
    if (clust != "agnes" | linkage != "ward") {
        print("Only agglomerative hierarchical clustering with ward linkage implemented.")
    }
    WeightedSimCluster = agnes(DistW, diss = TRUE, method = linkage)
    out = list(Dist1, Dist2, Weight, DistW, Clust = WeightedSimCluster)
    names(out) = c("Dist1", "Dist2", "Weight", "DistW", "Clust")
    AllCluster = list()
    if (AllClusters == TRUE) {
        for (a in 1:length(w)) {
            W1 = w[a] * Dist1 + (1 - w[a]) * Dist2
            AllCluster[[a]] = agnes(W1, diss = TRUE, method = linkage)
            names(AllCluster)[a] = paste("Weight_", w[a], sep = "")
        }
        out = list(Dist1, Dist2, Weight, DistW, Clust = WeightedSimCluster, 
            AllCluster)
        names(out) = c("Dist1", "Dist2", "Weight", "DistW", "Clust", 
            "Results")
    }
    class(out) = "WeightedSim"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
