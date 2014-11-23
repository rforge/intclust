\name{SharedSelectionMLP}
\alias{SharedSelectionMLP}
\title{
Intersection of genes and pathways over multiple methods for a selection of compounds
}
\description{
Finds the shared genes and pathways of outputs of the \code{Geneset.intersectSelection}
function.
}
\usage{
SharedSelectionMLP(DataMLP = NULL, names = NULL)
}
\arguments{
  \item{DataMLP}{
The output of \code{Geneset.intersectSelection} function.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The function \code{SharedSelectionMLP} will provide this information. Given the output of
the \code{Geneset.intersectSelection} function, it investigates how many genes and pathways
are expressed by each cluster per method, how many of these are shared over the methods and which
ones are shared including their respective p-values of each method and a mean p-value. This is very
handy to look into the shared genes and pathways of clusters that share many objects but also of 
those that only share only a few. Further, the result also includes the number of compounds per
cluster per method and how many of these are shared over the methods.
}
\author{
Marijke Van Moerbeke
}
\examples{
\dontrun{
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


MCF7_SelFPaths=PathwaysSelectionIter(L,Comps,names=names,GeneExpr=geneMat,nrclusters=7,
method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,
topG=NULL,GENESET=ListGO,sign=0.05,niter=2,fusionsLog=TRUE,WeightClust=TRUE)

MCF7_SelFPaths_Inter=Geneset.intersectSelection(MCF7_SelFPaths,sign=0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_SelShared=SharedSelectionMLP(MCF7_SelFPaths_Inter)
}
}
\seealso{
\code{\link{SharedLimma}},\code{\link{SharedMLP}},\code{\link{Shared}}
\code{\link{SharedSelection}},\code{\link{SharedSelectionLimma}}
}

