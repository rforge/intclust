\name{SharedComps}
\alias{SharedComps}
\title{
Intersection of compounds over multiple methods over compounds
}
\description{
The \code{SharedComps} function is an easy function to select the compounds over a number of methods.
To be used on cluster elements of the \code{DiffGenes} elements.
}
\usage{
SharedComps(List)
}
\arguments{
  \item{List}{
A list of the outputs of clusters for different methods.
}
}
\value{
A vector containing the shared compounds of all listed elements.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)
data(GeneInfo)
data(ListGO)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

L=list(MCF7_F,MCF7_T)
names=c('FP','TP')

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,
"limma",0.05,top=10)

Comps=SharedComps(list(MCF7_DiffGenes_FandT10$`Method 1`$"Cluster 1",
MCF7_DiffGenes_FandT10$`Method 2`$"Cluster 1"))

}
