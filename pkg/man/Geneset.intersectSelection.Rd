\name{Geneset.intersectSelection}
\alias{Geneset.intersectSelection}
\title{
Intersection over resulting gene sets of \code{PathwaysIterSelection} function
}
\description{
The function \code{Geneset.intersectSelection} performs the same procedure as 
\code{Geneset.intersect} but only for a specific selection of compounds and 
only for the cluster that contains the maximum number of compounds of the 
selection. It works on the output of \code{PathwaysIterSelection}.
}
\usage{
Geneset.intersectSelection(list.output, sign, names = NULL, 
seperatetables = FALSE, separatepvals = FALSE)
}
\arguments{
  \item{list.output}{
The output of the \code{PathwaysIterSelection} function
}
  \item{sign}{
The significance level to be handled for cutting of the pathways.
}
  \item{names}{
Optional. Names of the methods.
}
  \item{seperatetables}{
Logical. If TRUE, a separate element is created per cluster. containing
the pathways for each iteration.
}
  \item{separatepvals}{
Logical. If TRUE, the p-values of the each iteration of each pathway
in the intersection is given. If FALSE, only the mean p-value is provided.
}
}
\value{
The output is a list with an element per method. For each method, it is
portrayed for the cluster with the maximum number of compounds which 
pathways belong to the intersection over all iterations and their 
corresponding mean p-values.
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

Comps=FindCluster(L,nrclusters=7,select=c(1,4))
Comps

MCF7_Paths=PathwaysSelectionIter(L,Comps,names=names,GeneExpr=geneMat,
nrclusters=7,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = 
"GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,niter=2,fusionsLog=TRUE,
WeightClust=TRUE)

MCF7_Paths_Inter=Geneset.intersectSelection(MCF7_Paths,sign=0.05,
names=names,seperatetables=FALSE,separatepvals=FALSE)
}
}
