\name{DiffGenesSelection}
\alias{DiffGenesSelection}
\title{
Differential gene expression for a selection of compounds
}
\description{
The function \code{DiffGenesSelection} performs the same procedure as 
\code{DiffGenes} but only for a specific selection of compounds and 
only for the cluster that contains the maximum number of
compounds of the selection over the multiple methods.
}
\usage{
DiffGenesSelection(List, Selection, GeneExpr = geneMat, nrclusters = 7,
 method = "limma", sign = 0.05, top = NULL, fusionsLog = TRUE, 
 WeightClust = TRUE, names = NULL)
}
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of 
the list will be used as the reference.
}
  \item{Selection}{
The selection of objects to follow.
}
  \item{GeneExpr}{
The gene expression matrix of the objects.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{method}{
The method to applied to look for DE genes. For now, only the limma method
is available
}
  \item{sign}{
The significance level to be handled.
}
  \item{top}{
Overrules sign. The number of top genes to be shown.
}
  \item{fusionsLog}{
To be handed to \code{MatrixFunction}.
}
  \item{WeightClust}{
To be handed to \code{MatrixFunction}.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The returned value is a list with an element for each
method. Per method, there is only the cluster that
contains the maximum number of compounds and this contains
another list of the following elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)

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

}
\seealso{
\code{\link{DiffGenes}}
}
\keyword{Differential Gene Expression }
\keyword{ limma }
