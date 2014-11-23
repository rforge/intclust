\name{ProfilePlot}
\alias{ProfilePlot}
\title{
Plotting a gene profile
}
\description{
The function \code{ProfilePlot} shows the gene expression of a specific gene over the objects of the
data set and has the ability to separate the values a set of objects from the others. It can be used
for example to investigate the gene expression of a gene that was found to be shared between the
clusters over the methods. The values of the shared objects of these clusters are then plotted in red
at the left side of the plot while the other values gained a blue color and are put to the right. 
Further, it is optional to give the objects that belong to these clusters but are not common to
the cluster over the methods, a specific color for each method. This way, it is known how much
the gene expression of these objects resembles those of the shared components.
}
\usage{
ProfilePlot(Gene, Comps = Comps, GeneExpr = geneMat, Clusters = NULL,
 cols = NULL, AddLegend = TRUE, names = NULL,
  margins = c(8.1, 4.1, 1.1, 6.5), extra, ...)
}
\arguments{
  \item{Gene}{
The gene to be plotted.
}
  \item{Comps}{
The objects to be plotted or to be separated from the other objects.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
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
  \item{names}{
Optional but necessary when Addlegend=TRUE. Indicates the labels in the legend.
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
A plot of the gene profile. A distinction is made by the use of different colors
if not all objects are listed in the parameter Comps.
}
\author{
Marijke Van Moerbeke
}

\seealso{
\code{\link{ProfilePlotAll}}
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
str(MCF7_Shared10)

Comps=MCF7_Shared10$Which$`Cluster 2`$sharedcomps
Genes=MCF7_Shared10$Which$`Cluster 2`$sharedgenes


ProfilePlot(Gene=Genes[1],Comps,GeneExpr=geneMat,Clusters=NULL,
cols=NULL,margins=c(8.1,4.1,1.1,4.1),extra=2,names=names)
}
}
\keyword{ Gene profile }
