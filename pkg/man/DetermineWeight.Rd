\name{DetermineWeight}
\alias{DetermineWeight}
\title{
Determines a weight for weighted similarity clustering.
}
\description{
The function \code{DetermineWeight} determines an optimal weight to perform weighted
similarity clustering with. For each given weight, is each separate clustering
compared to the clustering on a weighted dissimilarity matrix and a Jaccard 
coefficient is calculated. The ratio of the Jaccard coefficients closets to one,
indicates an optimal weight.
}
\usage{
DetermineWeight(List, type = c("data", "clusters"), weight=seq(0, 1, by = 0.01), 
nrclusters = NULL, distmeasure = c("tanimoto", "tanimoto"), 
clust = "agnes", linkage = "ward", gap = FALSE, maxK = 50, 
names = c("B", "FP"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of two data matrices or a list of the outcomes of the \code{Cluster}
function.
}
  \item{type}{
Indicates whether data matrices or clustering results were given in List.
Should be one of "data" or "clusters".
}
  \item{weight}{
 Optional. A list of different weight combinations for the data sets in List.
 If NULL, the weight is sequence from 0 to 1  with steps of 0.1 and a result is produced 
 for each. It is further possible to fix weights for some data matrices and to
 let it vary randomly for the remaining data sets. An example is provided in the details.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in. This is necessary for the
computation of the Jaccard coefficient.
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard".
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
  \item{gap}{
Logical. Whether or not to calculate the gap statistic in the clustering
on each data matrix separately. Only if type="data".
}
  \item{maxK}{
The maximal number of clusters to consider in calculating the gap
statistic. Only if type="data".
}
  \item{names}{
The labels to give to the elements in List.
}
}
\details{
If the type of List is data, an hierarchical clustering is performed on each
data matrix separately. After obtaining clustering results for the two data
matrices, the distance matrices are extracted. If these are not calculated
with the same distance measure, they are normalized to be in the same range.
For each weight, a weighted linear combination of the distance matrices is
taken and hierarchical clustering is performed once again. The resulting
clustering is compared to each of the separate clustering results and a 
Jaccard coefficient is computed.The ratio of the Jaccard coefficients
closets to one, indicates an optimal weight. A plot of all the ratios is
produced with an extra indication for the optimal weight. 

The weight combinations should be provided as elements in a list. For three data
matrices an example could be: weights=list(c(0.5,0.2,0.3),c(0.1,0.5,0.4)). To provide
a fixed weight for some data sets and let it vary randomly for others, the element "x" 
indicates a free parameter. An example is weights=list(c(0.7,"x","x")). The weight 0.7 
is now fixed for the first data matrix while the remaining 0.3 weight will be divided over
the other two data sets. This implies that every combination of the sequence from 0 to 0.3
with steps of 0.1 will be reported and clustering will be performed for each.
}
\value{
The returned value is a list with four elements:
\item{Clust1}{The result of \code{Cluster} for the first element of List}
\item{Clust2}{The result of \code{Cluster} for the second element of List}
\item{Result}{A data frame with the Jaccard coefficients and their ratios 
for each weight}
\item{Weight}{The optimal weight}
}
\references{
PERUALILA-TAN, N., SHKEDY, Z., TALLOEN, W., GOEHLMANN, H. W. H.,
QSTAR Consortium, VAN MOERBEKE, M., KASIM, A., (in press). Weighted-Similarity Based
Clustering of Chemical Structure and Bioactivity Data in Early Drug Discovery. Journal of
Bioinformatics and Computational Biology. 
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

L=list(MCF7_F,MCF7_T)

MC7_Weight=DetermineWeight(L,type="clusters",weight=seq(0,1,by=0.01),nrclusters=7,
distmeasure=c("tanimoto","tanimoto"),clust="agnes",linkage="ward",gap=FALSE,
maxK=50,names=c("FP","TP"))

}
\seealso{
\code{\link{WeightedSimClust}}
}
\keyword{ Weighted Similarity Clustering }
\keyword{ Integrated Data Clustering }
