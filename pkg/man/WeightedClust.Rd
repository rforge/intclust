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
seq(1, 0, -0.1), WeightClust = 0.5, clust = "agnes", linkage = "ward")
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
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weight is sequence from 0 to 1  with steps of 0.1 and a result is produced 
 for each. It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
}
  \item{WeightClust}{
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
\details{
The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
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
weight=seq(1,0,-0.1),WeightClust=0.5,clust="agnes",linkage="ward")

}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}

\keyword{ Weighted Clustering }
