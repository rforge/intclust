\name{CECb}
\alias{CECb}
\title{
Complementary Ensemble Clustering - version b
}
\description{
Function \code{CECb} performs complementary ensemble clustering in which
in every iteration the total number of features are used in the clustering
procedure. However, the function is capable of cutting the resulting
dendrogram several times, each time into a different number of cluster.
}
\usage{
CECb(List, distmeasure = c("tanimoto", "tanimoto"), nrclusters = 
seq(5, 25, 1), weight = NULL, clust = "agnes", linkage = "ward", 
WeightClust = 0.5)
}
\arguments{
  \item{List}{
A list of data matrices. It is assumed the  rows are corresponding
 with the objects.
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard".
}
  \item{nrclusters}{
A sequence of numbers of clusters to cut the dendrogram in.
}
  \item{weight}{
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weight is sequence from 0 to 1  with steps of 0.1 and a result is produced 
 for each. It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
  \item{WeightClust}{
A weight for which the result will be put aside of the other results.
This was done for comparative reason and easy access.
}
}
\details{
Ensemble clustering is performed on each data matrix. This comes down to 
repeatedly applying hierarchical clustering. All features will be used 
in every iteration. Variation is inserted by not splitting the dendrogram
a single time into one specific number of clusters but multiple times and
for a range of numbers of clusters. Afterwards the two incidence matrices
are combined in a weighted sum and hierarchical clustering is performed
once more. More information can be found in Fodeh et al. (2013).

The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
The returned value is a list with the following four elements.
\item{Incidence}{The summed incidence matrices for each data matrix}
\item{IncidenceComb}{The co-association matrix after a weighted sum of the
elements of Incidence for each weight}
\item{Results}{The hierarchical clustering result for each element in 
IncidenceComb}
\item{Clust}{The result for the weight specified in Clustweight}
The value has class 'CEC'
}
\references{
FODEH, J. S., BRANDT, C., LUONG, B. T., HADDAD, A., SCHULTZ, M., MURPHY, T.,
KRAUTHAMMER, M. (2013). Complementary Ensemble Clustering of Biomedical Data.
J Biomed Inform. 46(3) pp.436-443.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)
MCF7_CECb=CECb(L,distmeasure=c("tanimoto","tanimoto"),nrclusters=seq(5,25,1),
clust="agnes",linkage="ward")
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}
\seealso{
\code{\link{CECa}},\code{\link{CECc}}
}
\keyword{ Complementary Ensemble Clustering }
\keyword{ Integrative Clustering }
