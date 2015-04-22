\name{WeightedSimClust}
\alias{WeightedSimClust}
\title{
Weighted similarity clustering
}
\description{
The \code{WeightedSimClust} function performs weighted 
similarity clustering. The input can be data matrices of which
the distance mattices are computed or clustering results where from
the distance matrices are extracted. An optimal weight is chosen with the \code{
ChooseWeight} function or can be specified by the user. With the found 
weight the distance matrices are linearly combined and hierarchical clustering
is performed. 
}
\usage{
WeightedSimClust(List, type = c("data", "clusters"), weight = seq(0, 1, 0.01),
 clust = "agnes", linkage = "ward", distmeasure = c("euclidean",
  "tanimoto"), gap = FALSE, maxK = 50, nrclusters = NULL, 
  names = c("B", "FP"), AllClusters = FALSE)
}

\arguments{
 \item{List}{
A list of data matrices of the same type. It is assumed the rows are corresponding
with the objects.
}
  \item{type}{
Indicates whether data matrices or clustering results were given in List.
Should be one of "data" or "clusters".
}
  \item{weight}{
One specific weight to perform clustering on or a list with different weight
combinations. If different weight combinations are provided, the function
\code{Chooseweight} is called and an optimal combination is chosen.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard".
}
  \item{gap}{
Logical. Whether or not to calculate the gap statistic in the clustering
on each data matrix separately. Only if type="data".
}
  \item{maxK}{
The maximal number of clusters to consider in calculating the gap
statistic. Only if type="data".
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
 \item{names}{
The labels to give to the elements in List.
}
  \item{AllClusters}{
Logical. Whether clustering should be performed for every weight.
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
The returned value is a list with four elements:
\item{Dist1}{The distance matrix of the first data object}
\item{Dist2}{The distance matrix of the second data object}
\item{Weight}{The optimal weight}
\item{DistW }{The weighted distance matrices for the optimal 
weight}
\item{Clust}{The resulting clustering}
If AllClusters was specified to be TRUE, a sixth element appears
containing the clustering results for all weights. The value has
class 'WeightedSimClust'
}
\references{
RAVINDRANATH, A. C.,PERUALILA-TAN, N., KASIM, A.,DRAKAKIS, G., LIGGI, S., 
BREWERTON, S. C.,MASON, D., BODKIN, M. J., EVANS, D. A., BHAGWAT, A.
TALLOEN, W., GOHLMANN, H. W. H., QSTAR Consortium, SHKEDY, Z., BENDER, A. 
(2015). Connecting gene expression data from connectivity map and in silico
target predictions for small molecule mechanism-of-action analysis.
Mol. BioSyst. Available at: <http://pubs.rsc.org/En/content/
articlelanding/2015/mb/c4mb00328d#!divAbstract>
}
\author{
Marijke Van Moerbeke
}
\examples{

data(fingerprintMat)
data(targetMat)

L=list(fingerprintMat,targetMat)

MCF7_WeightSim=WeightedSimClust(L,type="data", weight=seq(0,1,0.01),clust="agnes",
linkage="ward",distmeasure=c("tanimoto","tanimoto"), gap=FALSE,maxK=50,
nrclusters=7,names=c("FP","B"),AllClusters=FALSE)

}
\seealso{
\code{\link{DetermineWeight}}
}
\keyword{ Weighted Similarity Clustering }
\keyword{ Integrated Data Clustering }
