\name{CECa}
\alias{CECa}
\title{
Complementary Ensemble Clustering - version a
}
\description{
Function \code{CECa} performs complementary ensemble clustering in
which in every iteration the number of random samples taken is 
randomly set between m/2 and m-1 with m the total number of features.
The number of features to sample can also be prespecified by the user.
}
\usage{
CECa(List, distmeasure = c("tanimoto", "tanimoto"), t = 10, r = NULL, 
nrclusters = NULL, weight = NULL, clust = "agnes", linkage = "ward", 
Clustweight = 0.5)
}
\arguments{
  \item{List}{
A list of data matrices (max two matrices). It is assumed the 
rows are corresponding with the objects.
}
  \item{distmeasure}{
A character vector with the distance measure for each data matrix.
Should be one of "tanimoto", "euclidean", "jaccard".
}
  \item{t}{
The number of iterations.
}
  \item{r}{
Optional. The number of features to take for the random sample.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{weight}{
 Optional. A specific weight to give to the first co-association matrix.
 If NULL, the weight is sequence from 0 to 1 and a result is produced 
 for each.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
  \item{Clustweight}{
A weight for which the result will be put aside of the other results.
This was done for comparative reason and easy access.
}
}
\details{
Ensemble clustering is performed on each data matrix. This comes down to 
repeatedly applying hierarchical clustering. A random sample of features is taken 
in each application. Afterwards the two incidence matrices are combined in 
a weighted sum and hierarchical clustering is performed once more.
More information can be found in Fodeh et al. (2013).
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

MCF7_CECa=CECa(L,distmeasure=c("tanimoto","tanimoto"),t=25,r=NULL,,nrclusters=c(7,7),
clust="agnes",linkage="ward")
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}
\seealso{
\code{\link{CECb}},\code{\link{CECc}}
}
\keyword{ Complementary Ensemble Clustering }
\keyword{ Integrative Clustering }
