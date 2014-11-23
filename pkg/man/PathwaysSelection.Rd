\name{PathwaysSelection}
\alias{PathwaysSelection}
\title{
Pathway analysis for a selection of compounds
}
\description{
The function \code{PathwaysIterSelection} performs the same procedure as \code{Pathways}
but only for a specific selection of compounds and only for the cluster that contains
the maximum number of compounds of the selection.
}
\usage{
PathwaysSelection(List,Selection, GeneExpr = geneMat, nrclusters = 7, method = c("limma",
 "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
 topG = NULL, GENESET = ListGO, sign = 0.05, fusionsLog = TRUE,
  WeightClust = TRUE, names = NULL)
}
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of the
list will be used as the reference.
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
is available.
}
  \item{ENTREZID}{
Vector containing the ENTREZID's of the genes. If not provided, the rownames
of the matrix will be considered.
}
  \item{geneSetSource}{
The source for the getGeneSets function ("GOBP", "GOMF","GOCC", "KEGG"
or "REACTOME").
}
  \item{top}{
Overrules sign. The number of pathways to display for each cluster. 
If not specified, only the significant genes are shown.
}
  \item{topG}{
 Overrules sign. The number of top genes to be returned in the result.
 If not specified, only the significant genes are shown.
}
  \item{GENESET}{
Optional. Can provide own candidate gene sets.
}
  \item{sign}{
The significance level to be handled.
}
  \item{fusionsLog}{
To be handed to MatrixFunction.
}
  \item{WeightClust}{
To be handed to MatrixFunction.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
The output of this function is a list with an element for each method and per method, one for the
cluster with the maximum number of compounds. This element is again a list with the following
four elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Characteristics}{The found (top) characteristics of the feauture data}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
\item{Pathways}{A list with the element ranked.genesets.table which is a data frame containing
the genesets, their p-values and their descriptions. The second element is nr.genesets and
contains the used and total number of genesets.}
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

MCF7_SelFPaths=PathwaysSelection(L,Comps,GeneExpr=geneMat,nrclusters=7,
method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource =
"GOBP",top=NULL, topG=NULL,GENESET=ListGO,sign=0.05,fusionsLog=TRUE,
WeightClust=TRUE,names=names)
}
}
\seealso{
\code{\link{PathwaysSelectionIter}},\code{\link{PathwaysIter}}
}
\keyword{ Pathway Analysis }
