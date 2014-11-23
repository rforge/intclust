\name{SharedSelection}
\alias{SharedSelection}
\title{
Intersection of genes and pathways over multiple methods for a selection of compounds
}
\description{
The function \code{SharedSelection} performs the same procedure as Shared but only
for a specific selection of compounds and only for the cluster that contains the
maximum number of compounds of the selection. It works on the output of
\code{DiffGenesSelection} and \code{Geneset.intersectSelection}.
}
\usage{
SharedSelection(DataLimma = NULL, DataMLP = NULL, names = NULL)
}
\arguments{
  \item{DataLimma}{
The output of a \code{DiffGenesSelection} function.
}
  \item{DataMLP}{
The output of \code{Geneset.intersectSelection} function.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The result of the \code{SharedSelection} function is a list with 2 elements. 
The first element Table is a table indicating how many genes and/or pathways were
found to be differentially expressed and how many of these are shared. The table
also contains the number of compounds shared between the cluster of the different methods.
The second element Which is another list with one component. It consists of three vectors:
sharedcomps indicating which objects were shared across the methods, sharedgenes
represents the shared genes and the last one sharedpaths shows the shared pathways. The elements
pvalsgenes and pvalspaths contain the mean p-values of the shared genes and pathways.
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

MCF7_SelFPaths=PathwaysSelectionIter(L,Comps,names=names,GeneExpr=geneMat,nrclusters=7,
method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,
topG=NULL,GENESET=ListGO,sign=0.05,niter=2,fusionsLog=TRUE,WeightClust=TRUE)

MCF7_SelFPaths_Inter=Geneset.intersectSelection(MCF7_SelFPaths,sign=0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_SelShared=SharedSelection(MCF7_SelFLimma,MCF7_SelFPaths_Inter)
}
\seealso{
\code{\link{SharedLimma}},\code{\link{SharedMLP}},\code{\link{Shared}}
\code{\link{SharedSelectionLimma}},\code{\link{SharedSelectionMLP}}
}

