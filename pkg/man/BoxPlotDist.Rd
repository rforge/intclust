\name{BoxPlotDist}
\alias{BoxPlotDist}
\title{
Box Plots of one distance matrix categorized against another distance matrix.
}
\description{
Given two distance matrices, the function categorizes one distance matrix and
produces a box plot from the other distance matrix against the created 
categories. The option is available to choose one of the plots or to have both
plots. The function also works on outputs from ADEC and CEC functions which do 
not have distance matrices but incidence matrices.
}
\usage{
BoxPlotDist(Data1, Data2, lab1, lab2, limits1=NULL, limits2=NULL, type = 1)
}
\arguments{
  \item{Data1}{
The outcome of one of the integrated data cluster functions of the package.
}
  \item{Data2}{
The outcome of one of the integrated data cluster functions of the package.
}
  \item{lab1}{
The label to plot for Data1.
}
  \item{lab2}{
The label to plot for Data2.
}
  \item{limits1}{
The limits for the categories of Data1.
}
  \item{limits2}{
The limits for the categories of Data2.
}
  \item{type}{
The type of plot:
1 - Plot the values of Data1 versus the categories of Data2.
2 - Plot the values of Data2 versus the categories of Data1.
3 - Plot both types 1 and 2.
}
}
\value{
One/multiple box plots.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
data(targetMat)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)
MCF7_T = Cluster(targetMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)


BoxPlotDist(MCF7_F,MCF7_T,"FP","TP",limits1=c(0.3,0.7),limits2=c(0.3,0.7),type=1)
}
\keyword{ Box Plot}