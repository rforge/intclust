\name{ChooseWeight}
\alias{ChooseWeight}
\title{
Determines a weight for weighted similarity clustering.
}
\description{
The function \code{ChooseWeight} determines an optimal weight to perform weighted
similarity clustering with. For each given weight, is each separate clustering
compared to the clustering on a weighted dissimilarity matrix and a Jaccard 
coefficient is calculated. The ratio of the Jaccard coefficients closets to one,
indicates an optimal weight.
}
\usage{
ChooseWeight(List, type = c("data", "clusters"), w=seq(0, 1, by = 0.01), 
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
  \item{w}{
A sequence of weights.
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

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

L=list(MCF7_F,MCF7_T)

MC7_Weight=ChooseWeight(L,type="clusters",w=seq(0,1,by=0.01),nrclusters=7,
distmeasure=c("tanimoto","tanimoto"),clust="agnes",linkage="ward",gap=FALSE,
maxK=50,names=c("FP","B"))

}
\seealso{
\code{\link{WeightedSimClust}}
}
\keyword{ Weighted Similarity Clustering }
\keyword{ Integrated Data Clustering }
