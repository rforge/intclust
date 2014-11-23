\name{PathwaysCluster}
\alias{PathwaysCluster}
\title{
Pathway analysis for the output of \code{ChooseFeatures}
}
\description{
The function \code{PathwaysCluster} performs the same procedure as \code{Pathways}, 
\code{PathwaysIter} and \code{Geneset.intersect} but is specifically designed to 
perform on the output of the \code{ChooseFeatures} function.
}
\usage{
PathwaysCluster(Object,GeneExpr,topG,topP=NULL,method=c("limma","MLP"),
ENTREZID=NULL,geneSetSource = "GOBP",GENESET=ListGO,sign=0.05,niter=10)
}
\arguments{
  \item{Object}{
The output of the \code{ChooseFeatures} function.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
}
  \item{topG}{
Overrules sign. The number of top genes to be returned in the result.
If not specified, only the significant genes are shown.
}
  \item{topP}{
Overrules sign. The number of pathways to display for each cluster. 
If not specified, only the significant genes are shown.
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
The source for the getGeneSets function("GOBP", "GOMF","GOCC", "KEGG"
or "REACTOME").
}
  \item{GENESET}{
Optional. Can provide own candidate gene sets.
}
  \item{sign}{
The significance level to be handled
}
  \item{niter}{
The number of times to perform pathway analysis.
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
MCF7_Cluster=ChooseFeatures(Interactive=TRUE,LeadCpds=NULL,ClusterResult=MCF7_F,
ClusterColors=MCF7_F,BinData=list(fingerprintMat),Datanames=c("FP"),
geneMat,topChar = 20, topG = 20,nrclusters=7,N=2)

T=PathwaysCluster(MCF7_Cluster,geneMat,NULL,NULL,method=c("limma","MLP"),
ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",GENESET=ListGO,sign=0.05,
niter=2)

}
}
\keyword{ Pathway Analysis}

