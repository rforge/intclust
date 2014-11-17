\name{ChooseWeight}
\alias{ChooseWeight}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine the weight for weighted similarity clustering with jaccord coefficient
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChooseWeight(List, type = c("data", "clusters"), w = seq(0, 1, by = 0.01), nrclusters = NULL, distmeasure = c("tanimoto", "tanimoto"), clust = "agnes", linkage = "ward", gap = FALSE, maxK = 50, names = c("B", "FP"))
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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (List, type = c("data", "clusters"), w = seq(0, 1, by = 0.01), 
    nrclusters = NULL, distmeasure = c("tanimoto", "tanimoto"), 
    clust = "agnes", linkage = "ward", gap = FALSE, maxK = 50, 
    names = c("B", "FP")) 
{
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
        if (is.null(nrclusters)) {
            if (gap == FALSE) {
                stop("Specify a number of clusters of put gap to TRUE")
            }
            else {
                nrclusters = ceiling(mean(Clust1$k$Tibs2001SEmax, 
                  Clust2$k$Tibs2001SEmax))
            }
        }
        out <- list(Clust1, Clust2)
    }
    else {
        Dist1 = List[[1]]$DistM
        Dist2 = List[[2]]$DistM
        if ((all(range(Dist1) == c(0, 1)) & (range(Dist2)[1] != 
            0 & range(Dist2)[2] != 1)) | ((range(Dist1)[1] != 
            0 & range(Dist1)[2] != 1) & all(range(Dist2) == c(0, 
            1)))) {
            print("warning: the distance matrices do not have the same range. It might be that that one data set is bainary and the other is continuous. Standardzie the variables.")
        }
        out <- list(Clust1 = List[[1]], Clust2 = List[[2]])
    }
    labels <- c(paste("w", names[1], sep = ""), paste("J(sim", 
        names[1], ",simW)", sep = ""), paste("J(sim", names[2], 
        ",simW)", sep = ""), "Ratio")
    ResultsWeight <- data.frame(col1 = numeric(), col2 = numeric(), 
        col3 = numeric(), col4 = numeric())
    colnames(ResultsWeight) = labels
    for (i in 1:length(w)) {
        DistW = w[i] * Dist1 + (1 - w[i]) * Dist2
        hclust1 <- cutree(agnes(Dist1, diss = FALSE, method = linkage), 
            k = nrclusters)
        hclust2 <- cutree(agnes(Dist2, diss = FALSE, method = linkage), 
            k = nrclusters)
        hclust3 <- cutree(agnes(DistW, diss = FALSE, method = linkage), 
            k = nrclusters)
        Counts = function(clusterlabs1, clusterlabs2) {
            index = c(1:length(clusterlabs1))
            allpairs = combn(index, 2, simplify = FALSE)
            n11 = n10 = n01 = n00 = 0
            for (j in 1:length(allpairs)) {
                pair = allpairs[[j]]
                if (clusterlabs1[pair[1]] == clusterlabs1[pair[2]]) {
                  if (clusterlabs2[pair[1]] == clusterlabs2[pair[2]]) {
                    n11 = n11 + 1
                  }
                  else {
                    n10 = n10 + 1
                  }
                }
                else {
                  if (clusterlabs2[pair[1]] == clusterlabs2[pair[2]]) {
                    n01 = n01 + 1
                  }
                  else {
                    n00 = n00 + 1
                  }
                }
            }
            Jac = n11/(n11 + n01 + n10)
            return(Jac)
        }
        Jac1 = Counts(hclust1, hclust3)
        Jac2 = Counts(hclust2, hclust3)
        Ratio = Jac1/Jac2
        ResultsWeight[i, ] = c(w[i], Jac1, Jac2, Ratio)
    }
    Weight = ResultsWeight[which.min(abs(ResultsWeight[, "Ratio"] - 
        1)), 1]
    out[[3]] = ResultsWeight
    out[[4]] = Weight
    names(out) = c("Clust1", "Clust2", "Result", "Weight")
    q <- qplot(ResultsWeight[, 1], ResultsWeight[, 4], xlab = bquote(w[.(names[1])]), 
        ylab = bquote(J[(list(sim[.(names[1])], sim[W]))]/J[(list(sim[.(names[2])], 
            sim[W]))]))
    q <- q + geom_abline(intercept = 1, slope = 0, colour = "red", 
        lwd = 1.2) + geom_point(aes(Weight, ResultsWeight[ResultsWeight[, 
        1] == Weight, 4]), shape = 19, size = 2.5, colour = "red") + 
        scale_y_continuous(breaks = seq(0, max(ResultsWeight[, 
            4]), 1))
    q <- q + scale_x_continuous(breaks = c(seq(0, max(ResultsWeight[, 
        1]), 0.25), Weight))
    q
    print(q)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
