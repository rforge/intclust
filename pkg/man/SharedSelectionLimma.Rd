\name{SharedSelectionLimma}
\alias{SharedSelectionLimma}
\title{
Intersection of genes over multiple methods for a selection of compounds
}
\description{
Finds the shared genes over outputs of the \code{DiffGenesSelection} function.
}
\usage{
SharedSelectionLimma(DataLimma = NULL, names = NULL)
}
\arguments{
  \item{DataLimma}{
The output of a \code{DiffGenesSelection} function
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The result of the \code{SharedSelection} function is a list with two elements.
The first element Table is a table indicating how many genes were found to be
differentially expressed and how many of these are shared. The table also
contains the number of compounds shared between the clusters of the different
methods. The second element Which is another list with a component per cluster.
Each component consists of two vectors: sharedcomps indicating which objects
were shared across the method and sharedgenes represents the shared genes.
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

Comps=FindCluster(L,nrclusters=7,select=c(1,1))
Comps

MCF7_SelFLimma=DiffGenesSelection(L,Comps,names=names,GeneExpr=geneMat,nrclusters=7,
method="limma",sign=0.05,top=NULL,fusionsLog=TRUE,WeightClust=TRUE)

MCF7_SelSharedLimma=SharedSelectionLimma(MCF7_SelFLimma)
}
\seealso{
\code{\link{SharedLimma}},\code{\link{SharedMLP}},\code{\link{Shared}}
\code{\link{SharedSelection}},\code{\link{SharedSelectionMLP}}
}
