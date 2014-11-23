\name{ProfilePlotAll}
\alias{ProfilePlotAll}
\title{
Plotting multiple gene profiles
}
\description{
In \code{ProfilePlotAll}, the gene profiles of the significant genes for a specific
cluster are shown on 1 plot. Therefore, each gene is standardized by subtracting its
the mean.
}
\usage{
ProfilePlotAll(Genes = Genes, Comps = Comps, GeneExpr = geneMat, 
Raw = FALSE, Order = NULL, Color = NULL, nrclusters = 7, 
Clusters = NULL, cols = NULL, AddLegend = TRUE, 
margins = c(8.1, 4.1, 1.1, 6.5), extra = 5, ...)
}
\arguments{
  \item{Genes}{
The genes to be plotted.
}
  \item{Comps}{
The objects to be plotted or to be separated from the other objects.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
}
  \item{Raw}{
Logical. Should raw p-values be plotted?
}
  \item{Order}{
Optional. If the compounds are to set in a specific order of a specific method.
}
  \item{Color}{
The clustering result that determines the color of the labels of the
objects in the plot.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{Clusters}{
Optional. A list of clusters to which Obs belongs for each method considered.
The observations that are not part of Comps will be given a different color per
method.
}
  \item{cols}{
Optional. The color to use for the objects in Clusters for each method.
}
  \item{AddLegend}{
Optional. Whether a legend of the colors should be added to the plot.
}
  \item{margins}{
Optional. Margins to be used for the plot.
}
  \item{extra}{
The space between the plot and the legend.
}
  \item{\dots}{
Optional parameter to be handed to the plot function.
}
}
\value{
A plot which contains multiple gene profiles. A distinction is made between the values for the
objects in Comps and the others.
}
\author{
Marijke Van Moerbeke
}
\seealso{
 \code{\link{ProfilePlot}}
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

MCF7_Paths_FandT=PathwaysIter(L,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),
ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,
niter=2,fusionsLog=TRUE,WeightClust=TRUE,names=c("FP","TP"))

MCF7_Paths_intersection=Geneset.intersect(MCF7_Paths_FandT,0.05,names=names,
seperatetables=FALSE,separatepvals=FALSE)

MCF7_DiffGenes_FandT10=DiffGenes(list(MCF7_F,MCF7_T),geneMat,nrclusters=7,"limma",0.05,top=10)

MCF7_Shared10=Shared(DataLimma=MCF7_DiffGenes_FandT10,DataMLP=MCF7_Paths_intersection)

Comps=SharedComps(list(MCF7_DiffGenes_FandT10$`Method 1`$"Cluster 1",
MCF7_DiffGenes_FandT10$`Method 2`$"Cluster 1"))

MCF7_SharedGenes=FindGenes(DataLimma=MCF7_DiffGenes_FandT10,names=c("FP","TP"))

Genes=names(MCF7_SharedGenes[[1]])[-c(2,4,5)]

ListC=list(MCF7_DiffGenes_FandT10[[1]][[1]]$Compounds$LeadCpds,
MCF7_DiffGenes_FandT10[[2]][[1]]$Compounds$LeadCpds)

colscl=ColorPalette(colors=c("red","green","purple","brown","blue","orange"),ncols=9)

ProfilePlotAll(Genes=Genes,Comps=Comps,GeneExpr=geneMat,Raw=FALSE,Order=MCF7_F,
Color=NULL,nrcluster=7,Clusters=ListC,cols=colscl,AddLegend=TRUE,
usedgenes=Genes,margins=c(8.1,4.1,1.1,6.5),cex=0.75)
}
}
\keyword{ Gene Profile }

