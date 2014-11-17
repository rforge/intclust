\name{SNFc}
\alias{SNFc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Similarity Network Fusion - version c
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SNFc(List, distmeasure = c("tanimoto", "tanimoto"), NN = 20, alpha = 0.5, T = 20, clust = "agnes", linkage = "ward")
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
    AffinityMatrix.2 = function(Diff, K = 20, sigma = 0.5) {
        N = nrow(Diff)
        diag(Diff) = 0
        sortedColumns = as.matrix(t(apply(Diff, 2, sort)))
        finiteMean <- function(x) {
            mean(x[is.finite(x)])
        }
        means = apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
            .Machine$double.eps
        sum1 <- function(x, y) (x + y)
        Sig = outer(means, means, sum1)/3 + Diff/3 + .Machine$double.eps
        Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
        densities = dnorm(Diff, 0, sigma * Sig, log = FALSE)
        densities = densities * (0.5 * Sig * sqrt(2 * pi))
        W = (densities + t(densities))/2
        W = densities
        return(W)
    }
    AffM = vector("list", length(List))
    for (i in 1:length(List)) {
        AffM[[i]] = AffinityMatrix.2(DistM[[i]], NN, alpha)
    }
    snf.2 = function(Wall, K = 20, t = 20) {
        LW = length(Wall)
        normalize <- function(X) {
            NMatrix = matrix(0, dim(X)[1], dim(X)[2])
            for (i in 1:dim(X)[1]) {
                row = X[i, ]
                row[i] = 0
                D = sum(row)
                for (j in 1:dim(X)[2]) {
                  N = X[i, j]
                  NMatrix[i, j] = N/(2 * D)
                  if (i == j) {
                    NMatrix[i, j] = 1/2
                  }
                }
            }
            return(NMatrix)
        }
        PMatrix <- vector("list", LW)
        SMatrix <- vector("list", LW)
        nextW <- vector("list", LW)
        for (i in 1:LW) {
            PMatrix[[i]] = normalize(Wall[[i]])
        }
        for (i in 1:LW) {
            zero <- function(x) {
                s = sort(x, index.return = TRUE)
                x[s$ix[1:(length(x) - K)]] = 0
                return(x)
            }
            SMatrix[[i]] = matrix(0, dim(PMatrix[[i]])[1], dim(PMatrix[[i]])[2])
            for (k in 1:nrow(SMatrix[[i]])) {
                SMatrix[[i]][k, ] = zero(PMatrix[[i]][k, ])
            }
        }
        for (i in 1:t) {
            for (j in 1:LW) {
                sumWJ = matrix(0, dim(PMatrix[[j]])[1], dim(PMatrix[[j]])[2])
                for (k in 1:LW) {
                  if (k != j) {
                    sumWJ = sumWJ + PMatrix[[k]]
                  }
                }
                nextW[[j]] = SMatrix[[j]] \%*\% (sumWJ/(LW - 1)) \%*\% 
                  t(SMatrix[[j]])
            }
            for (j in 1:LW) {
                PMatrix[[j]] = nextW[[j]]
                PMatrix[[j]] = normalize(PMatrix[[j]])
            }
        }
        W = matrix(0, nrow(PMatrix[[1]]), ncol(PMatrix[[1]]))
        for (i in 1:LW) {
            W = W + PMatrix[[i]]
        }
        W = W/LW
        return(W)
    }
    SNF_FusedM = snf.2(AffM, NN, T)
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
