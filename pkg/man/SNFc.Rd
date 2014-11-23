\name{SNFc}
\alias{SNFc}
\title{
Similarity Network Fusion - version c
}
\description{
Function \code{SNFc} performs SNF but first a normalization over all objects
is performed before taking the k neighbours of each object as a subset in
 obtaining the kernel matrix.
}
\usage{
SNFc(List, distmeasure = c("tanimoto", "tanimoto"), NN = 20, alpha = 0.5,
 T = 20, clust = "agnes", linkage = "ward")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of data matrices of the same type.It is assumed the rows are corresponding
with the objects.
}
  \item{distmeasure}{
A vector of the distance measures to be used on each data matrix.
}
  \item{NN}{
The number of neighbours to be used in the procedure. 
}
  \item{alpha}{
The parameter epsilon. The value is recommended to be between 0.3 and 0.8.
}
  \item{T}{
The number of iterations.
}
  \item{clust}{
Choice of clustering function (character). Defaults to "agnes".
}
  \item{linkage}{
Choice of inter group dissimilarity (character). Defaults to "ward".
}
}
\value{
The returned value is a list with two elements:
 \item{FusedM }{The fused similarity matrix}
 \item{Clus}{The resulting clustering}
}
\references{
WANG, B., MEZLINI, M. A., DEMIR, F., FIUME, M., TU, Z., BRUDNO, M.,
HAIBE-KAINS, B., GOLDENBERG, A. (2014). Similarity Network Fusion
for aggregating data types on a genomic scale. Nature. 11(3) pp. 
333-337.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
L=list(fingerprintMat,targetMat)

MCF7_SNFc=SNFc(list(fingerprintMat,targetMat),distmeasure=c("tanimoto",
"tanimoto"),NN=10,alpha=0.5,T=20,clust="agnes",linkage="ward")
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}

\seealso{
\code{\link{SNFa}},\code{\link{SNFb}}
}
\keyword{ Similarity Network Fusion}
\keyword{ Integrative Clustering }