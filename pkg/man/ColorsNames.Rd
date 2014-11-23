\name{ColorsNames}
\alias{ColorsNames}
\title{
Function that annotates colors to their names
}
\description{
The \code{ColorsNames} function is used on the output of the \code{MatrixFunction}
and matches the cluster numbers indicated by the cell with the names of the colors. 
This is necessary to produce the plot of the \code{ComparePlot} function and is
therefore an internal function of this function but can also be applied separately.
}
\usage{
ColorsNames(MatrixColors, cols = Colors2)
}
\arguments{
  \item{MatrixColors}{
A matrix which is the output of the MatrixFunction function.
}
  \item{cols}{
The hex codes of the colors to be used.
}
}
\value{
A vector containing the hex code of the color that corresponds to each cell of 
the matrix. This function is called upon by the \code{ComparePlot} function.
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

MatrixColors=MatrixFunction(L,nrclusters=7,fusionsLog=TRUE,WeightClust=TRUE,
names=names)

Names=ColorsNames(MatrixColors,cols=Colors2)
}
\seealso{
 \code{\link{MatrixFunction}}
}
\keyword{ hex color codes }

