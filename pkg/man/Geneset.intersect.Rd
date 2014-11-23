\name{Geneset.intersect}
\alias{Geneset.intersect}
\title{
Intersection over resulting gene sets of \code{PathwaysIter} function
}
\description{
The function \code{Geneset.intersect} puts per method the results of the 
\code{PathwaysIter} function together for each cluster and takes the intersection
over the iterations per cluster per method. This is to see if over the 
different resamplings of the data, similar pathways were discovered.
}
\usage{
Geneset.intersect(list.output, sign, names = NULL, seperatetables = FALSE, 
separatepvals = FALSE)
}
\arguments{
  \item{list.output}{
The output of the \code{PathwaysIter} function.
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
The output is a list with an element per method. For each method, it
is portrayed per cluster which pathways belong to the intersection over
all iterations and their corresponding mean p-values.
}
\author{
Marijke Van Moerbeke
}
\seealso{
\code{\link{PathwaysIter}}
}\examples{
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

MCF7_Paths_FandT=PathwaysIter(list(MCF7_F,MCF7_T),GeneExpr=geneMat,nrclusters=7,
method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,
topG=NULL,GENESET=ListGO,sign=0.05,niter=2,fusionsLog=TRUE,WeightClust=TRUE,
names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=c("FP",
"TP"),seperatetables=FALSE,separatepvals=FALSE)

str(MCF7_Paths_intersection)
}
}
