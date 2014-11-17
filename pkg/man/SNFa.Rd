\name{SNFa}
\alias{SNFa}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Similarity Network Fusion - version a
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SNFa(List, distmeasure = c("tanimoto", "tanimoto"), NN = 20, alpha = 0.5, T = 20, clust = "agnes", linkage = "ward")
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
function (List, distmeasure = c("tanimoto", "tanimoto"), NN = 20, 
    alpha = 0.5, T = 20, clust = "agnes", linkage = "ward") 
{
    if (class(List) != "list") {
        stop("Data must be of type list")
    }
    if (alpha < 0.3 | alpha > 1) {
        print("Warning: alpha is recommended to be between 0.3 and 1 for the SNF method. Default is 0.5.")
    }
    if (clust != "agnes" | linkage != "ward") {
        print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting\n\t\t\t\t\t\tfused matrix.")
        clust = "agnes"
        linkage = "ward"
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
    DistM = vector("list", length(List))
    for (i in 1:length(List)) {
        DistM[[i]] = Distance(List[[i]], distmeasure[i])
    }
    AffM = vector("list", length(List))
    for (i in 1:length(List)) {
        AffM[[i]] = affinityMatrix(DistM[[i]], NN, alpha)
    }
    SNF_FusedM = SNF(AffM, NN, T)
    rownames(SNF_FusedM) = rownames(List[[1]])
    colnames(SNF_FusedM) = rownames(List[[1]])
    if (clust == "agnes" & linkage == "ward") {
        HClust = agnes(SNF_FusedM, diss = FALSE, method = linkage)
    }
    out = list(SNF_FusedM = SNF_FusedM, Clust = HClust)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
