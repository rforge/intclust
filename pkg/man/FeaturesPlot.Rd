\name{FeaturesPlot}
\alias{FeaturesPlot}
\title{
Plot of a selection of features
}
\description{
The function \code{FeaturesPlot} plots the binary data matrix for a 
selection of features. It is possible to separate between compounds
of interest and the other compounds. This is a visualization to see which
characterist are (not) expressed in a specific cluster.
}
\usage{
FeaturesPlot(LeadCpds, OrderedCpds, Features, Data, Color, nrclusters = 7, 
cols = Colors2, name = c("FP"))
}
\arguments{
  \item{LeadCpds}{
A character vector containing the compounds one wants to separate 
from the others.
}
  \item{OrderedCpds}{
A character vector containing all compounds and in which order thet should
be plotted after putting the LeadCpds aside. 
}
  \item{Features}{
A character vector containing the selection of features to be plotted.
}
  \item{Data}{
The data matrix the features are derived from.
}
  \item{Color}{
The clustering result that determines the color of the labels of the
objects in the plot.
}
  \item{nrclusters}{
The number of clusters to consider.
}
  \item{cols}{
The colors for the labels of the objects.
}
  \item{name}{
The overall name to give the features.
}
}
\value{
A plot indicating the values of the features of the LeadCpds in green 
and those of the others in blue.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
\dontrun{
MCF7_Interactive=ChooseFeatures(Interactive=TRUE,LeadCpds=NULL,ClusterResult=MCF7_T,
ClusterColors=MCF7_F,BinData=list(fingerprintMat),Datanames=c("FP"),geneMat,
topChar = 20, topG = 20,nrclusters=7,N=1)

Lead=MCF7_Interactive$"Choice 1"$Compounds$LeadCpds
Order=MCF7_Interactive$"Choice 1"$Compounds$OrderedCpds
Feat=MCF7_Interactive$"Choice 1"$Characteristic$FP

FeaturesPlot(LeadCpds=Lead,OrderedCpds=Order,Features=Feat,Data=fingerprintMat,Color=MCF7_F,
nrclusters=7,cols=Colors1,name=c("FP"))
}
}
\keyword{Features Plot}
