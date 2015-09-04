\name{ClusterPlot}
\alias{ClusterPlot}
\title{
Plot a dendrogram with leaves colored by a result of choice
}
\description{
The above described function \code{ClusterCols} is used in the function
Clusterplot which actually plots the dendrogram made by ClusterCols.
Further, given the outputs of any other functions, it is capable of 
selectimg the elements needed for ClusterCols.
}
\usage{
ClusterPlot(Data1, Data2, nrclusters = NULL, cols = NULL,plottype="new",
location=NULL, ...)
}
\arguments{
  \item{Data1}{
The resulting list method 1 which contains the dendrogram to
be colored.
}
  \item{Data2}{
The resulting list method 2 , i.e. the resulting clustering
on which the colors should be based.
}
  \item{nrclusters}{
Optional. The number of clusters to cut the dendrogram in. If not specified
the dendrogram will be drawn without colors to discern the different clusters.
}
  \item{cols}{
The colors for the clusters.
}
  \item{plottype}{
Should be one of "pdf","new" or "sweave". If "pdf", a location should be provided in
"location" and the figure is saved there. If "new" a new graphic device is opened and if "sweave",
the figure is made compatible to appear in a sweave or knitr document.
}
\item{location}{
If plottype is "pdf", a location should be provided in "location" and the figure is saved there. 
}
  \item{\dots}{
Other options which can be given to the plot function.
}
}
\details{
This function relies ion the \code{ClusterCols} function.
}
\value{
A plot of the dendrogram of method 1 with colored leaves of
which the colors are determined by method 2.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)
MCF7_T = Cluster(targetMat,type="data",distmeasure="tanimoto",normalize=FALSE,
method=NULL,clust="agnes",linkage="ward",gap=FALSE,maxK=55,StopRange=FALSE)

ClusterPlot(MCF7_T ,MCF7_F,nrcluster=7,cols=Colors2,plottype="new",location=NULL,
main="Clustering on Fingerprints: Dendrogram",ylim=c(-0.1,1.8))
}
\keyword{ Dendrogram }
