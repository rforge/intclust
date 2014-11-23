\name{SNFa}
\alias{SNFa}
\title{
Similarity Network Fusion - version a
}
\description{
The function \code{SNFa} performs similarity network fusion as
implemented by the package SNFtool. The overall method is 
described in the paper by Wang et al (2014).
}
\usage{
SNFa(List, distmeasure = c("tanimoto", "tanimoto"), NN = 20, alpha = 0.5,
 T = 20, clust = "agnes", linkage = "ward")
}
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

MCF7_SNFa=SNFa(list(fingerprintMat,targetMat),distmeasure=c("tanimoto",
"tanimoto"),NN=10,alpha=0.5,T=20,clust="agnes",linkage="ward")
}
\note{
For now, only hierarchical clustering with the \code{agnes} 
function and the ward link is implemented.
}

\seealso{
\code{\link{SNFb}},\code{\link{SNFc}}
}
\keyword{ Similarity Network Fusion}
\keyword{ Integrative Clustering }
