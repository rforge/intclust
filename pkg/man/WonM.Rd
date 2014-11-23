\name{WonM}
\alias{WonM}
\title{
Weighting on Membership
}
\description{
Weighting on membership is performed with the \code{WonM} function. The first step is to compute the
appropriate distance matrices for each data source and to use these for hierarchical clustering. This
is executed with the agnes function and the ward link. The user may specify a range of values for
the number of clusters to cut the resulting dendrograms in. For each value of number of clusters,
an incidence matrix is computed and these are added for each data source separately. Eventually,
the sums of the incidence matrices are joined together as well, resulting in one consensus matrix.
Hierarchical clustering is performed on the consensus matrix to obtain the final clustering result.
}
\usage{
WonM(List, distmeasure = c("tanimoto", "tanimoto"), nrclusters = 
seq(5, 25, 1), clust = "agnes", linkage = "ward")
}
\arguments{
  \item{List}{
A list of data matrices of the same type. It is assumed the rows are corresponding
with the objects.
}
  \item{distmeasure}{
A vector of the distance measures to be used on each data matrix.
}
  \item{nrclusters}{
A sequence of the number of clusters to cut the dendrogram in.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
}
\value{
The returned value is list with four elements:
\item{DistM}{A list with the distance matrix for each
data structure}
\item{ClustSep }{The hierarchical clustering result on each
data set}
\item{Consensus }{The computed consensus matrix over
all data sources}
\item{Clust}{The resulting clustering}
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_WonM=WonM(L,distmeasure=c("tanimoto","tanimoto"),nrclusters=seq(5,25),
clust="agnes",linkage="ward")
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}
\keyword{ Weighting on Membership }
\keyword{ Integrative Clustering }
