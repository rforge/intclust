\name{ClusterCols}
\alias{ClusterCols}
\title{
Internal Function of \code{Clusterplot}
}
\description{
The ClusterCols function is capable of coloring the leaves of a dendrogram of one method with
the colors (and thus cluster) the compounds had in the clustering of another method. 
}
\usage{
ClusterCols(x, Data, nrclusters = NULL, cols = Colors2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The dendrogram of method 1 to be colored
}
  \item{Data}{
The output of an agnes function, i.e. the resulting clustering of method 2 on
which the colors should be based.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The colors for the clusters.
}
}
\details{
The function helps to see how compounds have changed clusters between these two methods. 
Given a dendrogram to color of method 1 and the resulting clustering of method 2 
(the output of an agnes function), each leaf (compound) will be given the color of the cluster 
it belongs to in method 2. If groups of compounds with the same color appear in the dendrogram, 
this implies that the clustering of these is similar to the original clustering. If method 1
and method 2 are the same, the dendrogram is just colored by cluster.
}
\value{
A dendrogram of method 1 with colored leaves of which the colors are determined by method 2.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

d <- dendrapply(as.dendrogram(MCF7_F$Clust,hang=0.02),ClusterCols,MCF7_F$Clust,7,Colors2)
}
\seealso{
\code{\link{ClusterPlot}}
}
\keyword{ Dendrogram }


