\name{WeightedClust}
\alias{WeightedClust}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Weighted Clustering
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
WeightedClust(List, distmeasure = c("tanimoto", "tanimoto"), weight = seq(1, 0, -0.1), Clustweight = 0.5, clust = "agnes", linkage = "ward")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
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
function (List, distmeasure = c("tanimoto", "tanimoto"), weight = seq(1, 
    0, -0.1), Clustweight = 0.5, clust = "agnes", linkage = "ward") 
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
    DistM = vector("list", length(List))
    for (i in 1:length(List)) {
        DistM[[i]] = Distance(List[[i]], distmeasure[i])
        rownames(DistM[[i]]) = rownames(List[[i]])
        colnames(DistM[[i]]) = rownames(List[[i]])
    }
    if (length(List) != 2) {
        stop("Method implemented for 2 data matrices only.")
    }
    WeightedDist = list()
    for (i in 1:length(weight)) {
        alpha = weight[i]
        WeightedDist[[i]] = alpha * DistM[[1]] + (1 - alpha) * 
            DistM[[2]]
    }
    WeightedClust = list()
    for (i in 1:length(WeightedDist)) {
        WeightedClust[[i]] = agnes(WeightedDist[[i]], diss = TRUE, 
            method = linkage)
        names(WeightedClust)[i] = paste("Weight", weight[i], 
            sep = " ")
        if (weight[i] == Clustweight) {
            Clust = WeightedClust[[i]]
        }
    }
    out = list(Dist = DistM, WeightedDist = WeightedDist, Results = WeightedClust, 
        Clust = Clust)
    class(out) = "Weighted"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
