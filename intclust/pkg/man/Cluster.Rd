\name{Cluster}
\alias{Cluster}
\title{
Perforn clustering on a single data source
}
\description{
The function \code{Cluster} was written to perform clustering on
a single source of information, i.e one data matrix. The option is
available to compute the gap statistic to determine the optimal
number of clusters. 
}
\usage{
Cluster(Data, distmeasure = "tanimoto", clust = "agnes", linkage =
 "ward", gap = TRUE, maxK = 50)
}
\arguments{
  \item{Data}{
 A matrix containing the data. It is assumed the rows are corresponding
 zith the objects.
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character).
Should be one of "tanimoto", "euclidean", "jaccard".
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
  \item{gap}{
Logical. Indicator if gap statistics should be computed. 
Setting to $FALSE$ will greatly reduce the computation time.
}
  \item{maxK}{
The maximum number of clusters to be considered during the gap. 
}
}
\details{
The gap statistic is determined by the criteria described by the cluster package:
firstSEmax, globalSEmax, firstmax,globalmax, Tibs2001SEmax. The number of 
iterations is set to a default of 500. The implemented distances to be used for
the dissimilarity matrix are jaccard, tanimoto and euclidean. The jaccard distances
were computed with the \code{dist.binary(\ldots,method=1)} function in the ade4
package and the euclidean ones with the \code{daisy} function in again the cluster
package. The Tanimoto distances were implemented manually. 
}
\value{
The returned value is a list with two elements:
\item{DistM}{The distance matrix of the data matrix}
\item{Clust}{The resulting clustering}
If the gap option was indicated to be true, another 3 elements are joined to the list.
Clust\_gap contains the output from the function to compute the gap statistics and
gapdata is a subset of this output. Both can be used to make plots to visualize the 
gap statistic. The final component is k which is a matrix containing the optimal number
of clusters determined by each criterion mentioned earlier.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

}
\note{
For now, the
only option is to carry out agglomerative hierarchical clustering 
with the ward link as it was implemented in the \code{agnes}
function in the cluster package.
}
\keyword{ Clustering }
