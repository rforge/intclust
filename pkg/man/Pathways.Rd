\name{Pathways}
\alias{Pathways}
\title{
Pathway analysis for multiple clustering results
}
\description{
A pathway analysis per the cluster per method is conducted.
}
\usage{
Pathways(List, GeneExpr = geneMat, nrclusters = 7, method = c("limma", "MLP"),
 ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, topG = NULL, 
 GENESET = ListGO, sign = 0.05, fusionsLog = TRUE, WeightClust = TRUE, 
 names = NULL)
}
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of the
list will be used as the reference. The output of \code{ChooseFeatures} is also
accepted.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
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
\details{
After finding differently expressed genes, it can be investigated whether pathways are related to
those genes. This can be done with the help of the function \code{Pathways} which makes use of the 
\code{MLP} function of the MLP package. Given the output of a method, the cutree function is performed
which results into a specific number of clusters. For each cluster, the limma method is performed
comparing this cluster to the other clusters. This to obtain the necessary p-values of the genes.
These are used as the input for the \code{MLP} function to find interesting pathways. By default 
the candidate gene sets are determined by the \code{AnnotateEntrezIDtoGO} function. The default
source will be GOBP, but this can be altered.  Further, it is also possible to provide own candidate
gene sets in the form of a list of pathway categories in which each component contains a vector of
Entrez Gene identifiers related to that particular pathway. The default values for the minimum and
maximum number of genes in a gene set for it to be considered were used. For MLP this is 
respectively 5 and 100. If a list of outputs of several methods is provided as data input, the
cluster numbers are rearranged according to a reference method. The first method is taken as the 
reference and MatrixFunction is applied to get the correct ordering. When the clusters haven been
re-appointed, the pathway analysis as described above is performed for each cluster of each method.
}
\value{
The returned value is a list with an element per cluster per method. This element is again a
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

MCF7_PathsFandT=Pathways(L, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
"MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, 
topG = NULL, GENESET = ListGO, sign = 0.05,fusionsLog = TRUE, WeightClust = TRUE, 
 names =names)
 }
}
\seealso{
\code{\link{PathwaysIter}}
}
\keyword{ Pathway Analysis }

