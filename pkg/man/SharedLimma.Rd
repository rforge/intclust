\name{SharedLimma}
\alias{SharedLimma}
\title{
Intersection of genes over multiple methods
}
\description{
Finds the shared genes over outputs of the \code{DiffGenes} function.
}
\usage{
SharedLimma(DataLimma, names = NULL)
}
\arguments{
  \item{DataLimma}{
The output of a \code{DiffGenes} function
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The result of the \code{Shared} function is a list with two elements. The first element Table is a table
indicating how many genes were found to be differentially expressed and how many of these are shared. 
The table also contains the number of compounds shared between the clusters of the different methods.
The second element Which is another list with a component per cluster. Each component consists of two
vectors: sharedcomps indicating which objects were shared across the method and sharedgenes represents
 the shared genes.
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

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,"limma",0.05,top=10)

MCF7_SharedLimma=SharedLimma(DataLimma=MCF7_DiffGenes_FandT10)
str(MCF7_SharedLimma)
}
\seealso{
\code{\link{Shared}},\code{\link{SharedMLP}},\code{\link{SharedSelection}}
\code{\link{SharedSelectionLimma}},\code{\link{SharedSelectionMLP}}
}

