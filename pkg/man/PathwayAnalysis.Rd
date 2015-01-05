\name{PathwayAnalysis}
\alias{PathwayAnalysis}
\title{
Pathway Analysis
}
\description{
The \code{PathwayAnalysis} function combines the functions \code{PathwaysIter} and \code{Geneset.intersect} such that only one function should be called.
}
\usage{
PathwayAnalysis(List, Selection=NULL, GeneExpr = geneMat, nrclusters = NULL,
method = c("limma", "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP",
top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, niter = 10, 
fusionsLog = TRUE, WeightClust = TRUE, names = NULL,seperatetables=FALSE,
separatepvals=FALSE)
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

MCF7_PathsFandT=PathwaysAnalysis(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
topG = NULL, GENESET = ListGO, sign = 0.05,niter=2,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names,seperatetables=FALSE,separatepvals=FALSE)
}
}
\keyword{ Pathway Analysis}

