\name{Pathways.2}
\alias{Pathways.2}
\title{
Pathways analysis for a list of length one
}
\description{
The \code{Pathways.2} function finds pathway analysis for a list of one
element. It is also an internal function of \code{Pathways} if only one
method is specified.
}
\usage{
Pathways.2(Data, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
topG = NULL, GENESET = ListGO, sign = 0.05,names=NULL)
}
\arguments{
  \item{Data}{
One clustering result of one of the integrated data clustering functions
of the package.The output of \code{ChooseFeatures} is also
accepted.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet  of the objects.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{method}{
The method to applied to look for DE genes. For now, only the limma method
is available
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
  \item{names}{
Optional. Name of the method.
}
}
\value{
The returned value is a list with an element per cluster. This element is again a
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
data(geneMat)
data(GeneInfo)
data(ListGO)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

MCF7_PathsF=Pathways.2(MCF7_F, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
topG = NULL, GENESET = ListGO, sign = 0.05)
}
}
\seealso{
\code{\link{Pathways}}
}
\keyword{ Pathway Analysis }

