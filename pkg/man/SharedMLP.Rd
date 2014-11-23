\name{SharedMLP}
\alias{SharedMLP}
\title{
Intersection of genes and pathways over multiple methods
}
\description{
Finds the shared genes and pathways of outputs of the \code{Geneset.intersect} function.
}
\usage{
SharedMLP(DataMLP,names=NULL)
}
\arguments{
  \item{DataMLP}{
The output of \code{Geneset.intersect} function.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The function \code{SharedMLP} will provide this information. Given the outputs of the
\code{Geneset.intersect} function, it investigates how many genes and pathways
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

MCF7_Paths_FandT=PathwaysIter(L,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),
ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,
niter=2,fusionsLog=TRUE,WeightClust=TRUE,names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_SharedMLP=SharedMLP(DataMLP=MCF7_Paths_intersection)
str(MCF7_SharedMLP)
}

}

\seealso{
\code{\link{SharedLimma}},\code{\link{Shared}},\code{\link{SharedSelection}}
\code{\link{SharedSelectionLimma}},\code{\link{SharedSelectionMLP}}
}
