\name{Distance}
\alias{Distance}
\title{
Distance function
}
\description{
The \code{Distance} function was written calculates the distances between the data objects. The included
distance measures are euclidean for continuous data and the tanimoto coefficient or jaccard index for binary data.
}
\usage{
Distance(Data,distmeasure=c("tanimoto","jaccard","euclidean","hamming","cont tanimoto"),
normalize=FALSE,method=NULL)
}
\arguments{
  \item{Data}{
A data matrix. It is assumed the  rows are corresponding with the objects.
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character).
Should be one of "tanimoto", "euclidean", "jaccard","hamming","cont tanimoto".
}
  \item{normalize}{
Logical. Indicates whether to normalize the distance matrices or not.
This is recommended if different distance types are used. More details
on standardization in \code{Normalization}.
}
  \item{method}{
A method of normalization. Should be one of "Quantile","Fisher-Yates",
"standardize","Range" or any of the first letters of these names.
}
}
\details{
The euclidean distance distance is included for continuous matrices while
for binary matrices, one has the choice of either the jaccard index, the
tanimoto coeffcient or the hamming distance. The hamming distance is obtained
by applying the \code{hamming.distance} function of the \pkg{e1071} package.
It will compute the hamming distance between the rows of the data matrix. The
hamming distance counts the number of times where two rows differ in their
zero and one values. The Jaccard index is calcaluted as determined by the
formula of the \code{dist.binary} function in the \pkg{a4} package and the
tanimoto coefficient as described by \cite{Li2011}. For both, first the
similarity is calculated as\deqn{s=frac{n11}{n11+n10+n01}}
with n11 the number of features the 2 compounds have in common, n10
the number of features of the first compound and n01 the number of features
of the second compound. These similarities are converted to distances by:
\deqn{J=\sqrt{1-s}}
for the jaccard index and by:
\deqn{T=1-s}
for the tanimoto coefficient. The lower the similarity values s are, the more features are
shared between the two objects and the more alike they are. Since clustering is
based on dissimilarity, the conversion to distances is performed.
If normalize=TRUE and the distance meausure is euclidean, the data matrix is
normalized beforehand. Further, a version of the tanimoto coefficient is also available for continuous data.

}
\value{
The returned value is a distance matrix.
}
\references{
LI, Y., TU, K., ZHENG, S., WANG, J., LI, Y., LI, X. (2011). Association of
Feature Gene Expression with Structural Fingerprints of Chemical Compounds.
Journal of Bioinformatics and Computational biology. 9(4). pp. 503-519.
MAECHLER, M., ROUSSEEUW, P., STRUYF, A., HUBERT, M. (2014).
cluster: Cluster Analysis Basics and Extensions. R package version
1.15.3.
TALLOEN, W., VERBEKE, T. (2011). a4: Automated Affymetrix Array
Analysis Umbrella Package. R package version 1.14.0

}
\author{
Marijke Van Moerbeke
}
\examples{
data(fingerprintMat)
Dist_F=Distance(fingerprintMat,distmeasure="tanimoto",normalize=FALSE,method=NULL)
}
\keyword{ Distances }