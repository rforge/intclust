\name{PathwaysIter}
\alias{PathwaysIter}
\title{
Iterations of the pathway analysis
}
\description{
The MLP method to perform pathway analysis is based on resampling of the data. Therefore it is
recommended to perform the pathway analysis multiple times to observe how much the results are
influenced by a different resample. The function \code{PathwaysIter} performs the pathway analysis
as described in \code{Pathways} a specified number of times. The input can be one data set or
a list as in \code{Pathway.2} and \code{Pathways}.
}
\usage{
PathwaysIter(List, Selection=NULL, GeneExpr = geneMat, nrclusters = NULL, 
method = c("limma", "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP",
top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, niter = 10, 
fusionsLog = TRUE, WeightClust = TRUE, names = NULL)
}
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of the
list will be used as the reference.
}
 \item{Selection}{
If pathway analysis should be conducted for a specific selection of compounds, this selection can be provided here.
Selection can be of the type "character" (names of the compounds) or "numeric" (the number of specific cluster).
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
  \item{niter}{
The number of times to perform pathway analysis.
}
  \item{fusionsLog}{
To be handed to ReorderToReference.
}
  \item{WeightClust}{
To be handed to ReorderToReference.
}
  \item{names}{
Optional. Names of the methods.
}
}
\value{
This element is again a
list with the following four elements:
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

MCF7_PathsFandT=PathwaysIter(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
topG = NULL, GENESET = ListGO, sign = 0.05,niter=2,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names)
}
}
\keyword{ Pathway Analysis}

