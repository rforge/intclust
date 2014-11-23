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
ClusterPlot(Data1, Data2, nrclusters = 7, cols = Colors2, ...)
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
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The colors for the clusters.
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
\seealso{
\code{\link{ClusterCols}}
}
\examples{
data(fingerprintMat)
data(targetMat)
data(Colors2)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

ClusterPlot(MCF7_T ,MCF7_F,nrcluster=7,cols=Colors2,main="Clustering on Fingerprints: Dendrogram",
ylim=c(-0.1,1.8))
}
\keyword{ Dendrogram }
