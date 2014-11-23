\name{LabelCols}
\alias{LabelCols}
\title{
Internal function of Labelplot
}
\description{
Sometimes only one or two particular clusters are of interest if for example these
contain an object of great importance. Then it would be nice not to color all
clusters but just the compounds of these clusters and see where these are now.
This can be done with the function \code{LabelCols}.
}
\usage{
LabelCols(x, Sel1, Sel2 = NULL, col1, col2 = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
The dendrogram to be colored.
}
  \item{Sel1}{
The selection of objects to be colored.
}
  \item{Sel2}{
An optional second selection to be colored.
}
  \item{col1}{
The color for the first selection.
}
  \item{col2}{
The color for the optional second selection.~~
}
}
\value{
A dendrogram of which the leaves of the selection(s) are colored.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

ClustF_6=cutree(MCF7_F$Clust,6)

SelF=rownames(fingerprintMat)[ClustF_6==6]
SelF

d_F <- dendrapply(as.dendrogram(MCF7_F$Clust,hang=0.02),LabelCols,
Sel1=SelF,Sel2=NULL,col1="darkorchid3")

d_F
}
\keyword{ Dendrogran }

