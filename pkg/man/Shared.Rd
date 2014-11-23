\name{Shared}
\alias{Shared}
\title{
Intersection of genes and pathways over multiple methods
}
\description{
It is interesting to investigate exactly which and how many differently expressed genes and pathways
are shared by the clusters over the different methods. The function \code{Shared} will provide this
information. Given the outputs of the \code{DiffGenes} function and/or the \code{PathwaysIter} function,
it investigates how many genes and/or pathways are expressed by each cluster per method, how many
of these are shared over the methods and which ones are shared including their respective p-values
of each method and a mean p-value. This is very handy to look into the shared genes and pathways
of clusters that share many objects but also of those that only share only a few. Further, the result
also includes the number of compounds per cluster per method and how many of these are shared
over the methods.
}
\usage{
Shared(DataLimma = NULL, DataMLP = NULL, names = NULL)
}
\arguments{
  \item{DataLimma}{
The output of a \code{DiffGenes} function.
}
  \item{DataMLP}{
The output of \code{Geneset.intersect} function.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The result of the \code{Shared} function is a list with two elements. The first element Table is a table
indicating how many genes and/or pathways were found to be differentially expressed and how
many of these are shared. The table also contains the number of compounds shared between the
clusters of the different methods. The second element Which is another list with a component per
cluster. Each component consists of three vectors: sharedcomps indicating which objects were shared
across the methods, sharedgenes represents the shared genes and the last one sharedpaths shows the
shared pathways. The elements pvalsgenes and pvalspaths contain the mean p-values of the shared
genes and pathways.
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

MCF7_Paths_FandT=PathwaysIter(L,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),
ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,
niter=2,fusionsLog=TRUE,WeightClust=TRUE,names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,"limma",0.05,top=10)

MCF7_Shared10=Shared(DataLimma=MCF7_DiffGenes_FandT10,DataMLP=MCF7_Paths_intersection)
str(MCF7_Shared10)
}
\seealso{
\code{\link{SharedLimma}},\code{\link{SharedMLP}},\code{\link{SharedSelection}}
\code{\link{SharedSelectionLimma}},\code{\link{SharedSelectionMLP}}
}

