\name{WeightedClust}
\alias{WeightedClust}
\title{
Weighted clustering
}
\description{
Weighted clustering is performed with the function \code{WeightedClust}. Given a
list of the data matrices, a dissimilarity matrix is computed of each with
the provided distance measures. These matrices are then combined resulting
in a weighted dissimilarity matrix. Hierarchical clustering is performed
on this weighted combination with the agnes function and the ward link
}
\usage{
WeightedClust(List, distmeasure = c("tanimoto", "tanimoto"), weight = 
seq(1, 0, -0.1), Clustweight = 0.5, clust = "agnes", linkage = "ward")
}
\arguments{
  \item{List}{
A list of data matrices of the same type.It is assumed the rows are corresponding
with the objects.
}
  \item{distmeasure}{
A vector of the distance measures to be used on each data matrix.
}
  \item{weight}{
The weight is a sequence from 0 to 1 and a result is produced for each by default.
A specific weight can be provided by the user.
}
  \item{Clustweight}{
A weight for which the result will be put aside of the other results. This was
done for comparative reason and easy access.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
}
\value{
The returned value is a list of four elements:
\item{DistM}{A list with the distance matrix for each
data structure}
\item{WeightedDist }{A list with the weighted distance matrices}
\item{Results}{The hierarchical clustering result for each element in 
IncidenceComb}
\item{Clust}{The result for the weight specified in Clustweight}
The value has class 'Weighted'
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_Weighted=WeightedClust(L, distmeasure=c("tanimoto","tanimoto"),
weight=seq(1,0,-0.1),Clustweight=0.5,clust="agnes",linkage="ward")

}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}

\keyword{ Weighted Clustering }
