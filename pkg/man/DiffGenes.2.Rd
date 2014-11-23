\name{DiffGenes.2}
\alias{DiffGenes.2}
\title{
Differential gene expression for a list of length one
}
\description{
The \code{DiffGenes.2} function finds differential gene expression
for a list of one element. It is also an internal function of \code{DiffGenes}
if only one method is specified.
}
\usage{
DiffGenes.2(Data, GeneExpr = geneMat, nrclusters = 7, method = "limma", 
sign = 0.05, top = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
One clustering result of one of the integrated data clustering functions
of the package.
}
  \item{GeneExpr}{
The gene expression matrix or ExpressionSet of the objects.
}
  \item{nrclusters}{
The number of clusters to cut the dendrogram in.
}
  \item{method}{
The method to applied to look for DE genes. For now, only the limma method
is available.
}
  \item{sign}{
The significance level to be handled.
}
  \item{top}{
Overrules sign. The number of top genes to be shown.
}
}
\details{
For each cluster, it is investigated which genes are differential expressed by
the limma method. The cluster is tested versus all other compounds combined.
}
\value{
The returned value is a list with an element per cluster.
Each element is a list with the following elements:
\item{Compounds}{A list with the elements LeadCpds (the compounds of interest) and 
OrderedCpds (all compounds in the order of the clustering result)}
\item{Genes}{A list with the elements TopDE (a table with information on the top
genes) and AllDE (a table with information on all genes)}
}
\references{
SMYTH, G. K. (2004). Linear models and empirical Bayes methods for assessing
differential expression in microarray experiments. Statistical Applications
in Genetics and Molecular Biology. 3(1).
}
\author{
Marijke Van Moerbeke
}
\seealso{
\code{\link{DiffGenes}}
}
\examples{
data(fingerprintMat)
data(geneMat)

MCF7_F = Cluster(fingerprintMat,distmeasure="tanimoto",clust="agnes",linkage="ward",
gap=FALSE,maxK=55)

MCF7_F_DE = DiffGenes.2(Data=MCF7_F,GeneExpr=geneMat,nrclusters=7,method="limma",
sign=0.05,top=10)

}
\keyword{Differential Gene Expression }
\keyword{ limma }
