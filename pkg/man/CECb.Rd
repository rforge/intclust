\name{CECb}
\alias{CECb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Complementary Ensemble Clustering - version b
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CECb(List, distmeasure = c("tanimoto", "tanimoto"), nrclusters = seq(5, 25, 1), weight = NULL, clust = "agnes", linkage = "ward", Clustweight = 0.5)
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
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{clust}{
%%     ~~Describe \code{clust} here~~
}
  \item{linkage}{
%%     ~~Describe \code{linkage} here~~
}
  \item{Clustweight}{
%%     ~~Describe \code{Clustweight} here~~
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
    25, 1), weight = NULL, clust = "agnes", linkage = "ward", 
    Clustweight = 0.5) 
{
    if (class(List) != "list") {
        stop("Data must be of type lists")
    }
    if (length(List) != 2) {
        stop("Method is just implemented for 2 data modalities.")
    }
    if (is.null(nrclusters)) {
        stop("Give a number of cluters to cut the dendrogram into for each data modality.")
    }
    if (clust != "agnes" | linkage != "ward") {
        print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting\n\t\t\t\t\t\tcoassociation matrix.")
        clust = "agnes"
        linkage = "ward"
    }
    A1 = List[[1]]
    A2 = List[[2]]
    A = list(A1, A2)
    Incidence = list()
    for (i in 1:length(List)) {
        Incidence[[i]] = matrix(0, dim(List[[i]])[1], dim(List[[i]])[1])
        rownames(Incidence[[i]]) = rownames(List[[i]])
        colnames(Incidence[[i]]) = rownames(List[[i]])
    }
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
    DistM = list()
    for (i in 1:length(A)) {
        DistM[[i]] = Distance(A[[i]], distmeasure[i])
    }
    HClust_A = list()
    for (i in 1:length(DistM)) {
        HClust_A[[i]] = agnes(DistM[[i]], diss = TRUE, method = linkage)
    }
    for (k in 1:length(nrclusters)) {
        MembersofClust = list()
        for (i in 1:length(HClust_A)) {
            Temp = cutree(HClust_A[[i]], nrclusters[k])
            MembersofClust[[i]] = matrix(0, dim(List[[i]])[1], 
                dim(List[[i]])[1])
            for (l in 1:length(Temp)) {
                label = Temp[l]
                sameclust = which(Temp == label)
                MembersofClust[[i]][l, sameclust] = 1
            }
            Incidence[[i]] = Incidence[[i]] + MembersofClust[[i]]
        }
    }
    if (is.null(weight)) {
        weight = seq(1, 0, -0.1)
    }
    IncidenceComb = list()
    for (i in 1:length(weight)) {
        IncidenceComb[[i]] = weight[i] * Incidence[[1]] + (1 - 
            weight[i]) * Incidence[[2]]
    }
    CEC = list()
    for (i in 1:length(IncidenceComb)) {
        CEC[[i]] = agnes(IncidenceComb[[i]], diss = FALSE, method = linkage)
        names(CEC)[i] = paste("Weight", weight[i], sep = "")
        if (weight[i] == Clustweight) {
            Clust = CEC[[i]]
        }
    }
    out = list(Incidence = Incidence, IncidenceComb = IncidenceComb, 
        Results = CEC, Clust = Clust)
    class(out) = "CEC"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
