\name{CompareSvsM}
\alias{CompareSvsM}
\title{
Comparison of clustering results for the single and multiple source clustering.
}
\description{
A visual comparison of all methods is handy to see which compounds will always
cluster together independent of the applied methods. The function \code{CompareSvsM} 
plots the \code{ComparePlot} of the single source clustering results on the left and 
that of the multiple source clustering results on the right such that a visual
comparison is possible. 
}
\usage{
CompareSvsM(ListS,ListM, nrclusters = NULL, cols = Colors2, fusionsLogS=FALSE,
fusionsLogM=FALSE,WeightClustS=FALSE,WeightClustM=FALSE, namesS = NULL,
namesM=NULL, margins = c(8.1, 3.1, 3.1, 4.1), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ListS}{
A list of the outputs from the single source clusterings to be compared. The first element of the
list will be used as the reference.
}
  \item{ListM}{
A list of the outputs from the multiple source clusterings to be compared. The first element of the
list will be used as the reference.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The hex codes of the colors to be used.
}
  \item{fusionsLogS}{
The fusionslog parameter for the elements in ListS. To be handed to \code{ReorderToReference}.
}
  \item{fusionsLogM}{
The fusionsLog parameter for the elements in ListM. To be handed to \code{ReorderToReference}.
}
  \item{WeightClustS}{
The WeightClust parameter for the elements in ListS. To be handed to \code{ReorderToReference}.
}
  \item{WeightClustM}{
The WeightClust parameter for the elements in ListM. To be handed to \code{ReorderToReference}.
}
  \item{namesS}{
Optional. Names of the single source clusterings to be used as labels for the columns.
}
  \item{namesM}{
Optional. Names of the multiple source clusterings to be used as labels for the columns.
}
  \item{margins}{
Optional. Margins to be used for the plot.
}
  \item{\dots}{
Other options which can be given to the \code{color2D.matplot} function.
}
}
\details{
This function relies on \code{ComparePlot} to plot both the results of the single source
clusterings as the multiple source clusterings. 
}
\value{
The returned value is a plot with on the left the comparison over the objects in ListS and on the right
a comparison over the objects in ListM.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

L=list(fingerprintMat,targetMat)

MCF7_W=WeightedClust(L, distmeasure=c("tanimoto","tanimoto"),
weight=seq(1,0,-0.1),WeightClust=0.5,clust="agnes",linkage="ward")

ListM=list(MCF7_W)
namesM=seq(1.0,0.0,-0.1)

ListS=list(MCF7_F,MCF7_T)
namesS=c("FP","TP")

CompareSvsM(ListS,ListM,nrclusters=7,cols=Colors2,fusionsLogS=FALSE,
fusionsLogM=FALSE,WeightClustS=FALSE,WeightClustM=FALSE,namesS,
namesM,reverse=FALSE)
}
\seealso{
\code{\link{ComparePlot}}
}
