\name{ComparePlot}
\alias{ComparePlot}
\title{
Comparison of clustering results over multiple results
}
\description{
A visual comparison of all methods is handy to see which compounds will always
cluster together independent of the applied methods. To this aid the function 
\code{ComparePlot} has been written.
}
\usage{
ComparePlot(List, nrclusters = 7, cols = Colors2, fusionsLog = FALSE, 
WeightClust = FALSE, names = NULL, reverse = FALSE, 
margins = c(8.1, 3.1, 3.1, 4.1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
A list of the outputs from the methods to be compared. The first element of the
list will be used as the reference.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{cols}{
The hex codes of the colors to be used.
}
  \item{fusionsLog}{
To be handed to \code{MatrixFunction}.
}
  \item{WeightClust}{
To be handed to \code{MatrixFunction}n.
}
  \item{names}{
Optional. Names of the methods to be used as labels for the columns.
}
  \item{reverse}{
Logical. Should the last element of the List be used as reference? 
Mostly to be used for CEC or weighted clustering.
}
  \item{margins}{
Optional. Margins to be used for the plot.
}
  \item{\dots}{
Other options which can be given to the \code{color2D.matplot} function.
}
}
\details{
This functions makes use of the functions \code{MatrixFunction} and \code{Colorsnames}. 
Given a list with the outputs of several methods, the first step is to call upon 
\code{ MatrixFunction} and to produce a matrix of which the columns are ordered according
to the ordering of the objects of the first method in the list. Each cell represent the
number of the cluster the object belongs to for a specific method indicated by the rows.
The clusters are arranged in such a way that these correspond to that one cluster of the
referenced method that they have the most in common with. The function \code{color2D.matplot}
produces a plot of this matrix but needs a vector indicating the names of the colors to be used.
This is where \code{ColorsNames} comes in. A vector of the color names of the output of the
\code{MatrixFunction} is created and handed to \code{color2D.matplot}. It is optional to adjust
the margins of the plot and to give a vector with the names of the methods which will be used
as labels for the rows in the plot. The labels for the columns are the names of the object in
the order of clustering of the referenced method. Further, the similarity measures of the methods
compared to the reference will be computed and shown on the right side of the plot.
}
\value{
A plot which translates the matrix output of the function \code{MatrixFunction} in which the columns
represent the objects in the ordering the referenced method and the rows the outputs of the 
given methods. Each cluster is given a distinct color. This way it can be easily observed which
objects will cluster together. The labels on the right side of the plot are the similarity
measures computed by \code{SimilarityMeasure}.
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

L=list(MCF7_F,MCF7_T)
names=c("FP","TP")

ComparePlot(L,nrclusters=7,cols=Colors2,fusionsLog=TRUE,WeightClust=TRUE,names=names,
margins=c(9.1,4.1,4.1,4.1))
}
\seealso{
\code{\link{MatrixFunction}},\code{\link{ColorsNames}}
}
