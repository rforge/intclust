\name{Distance}
\alias{Distance}
\title{
Distance function
}
\description{
The \code{Distance} function was written calculates the distances between the data objects. The included
distance measures are euclidean for continuous data and Tanimoto or Jaccard for binary data.
}
\usage{
Distance(Data, distmeasure = "tanimoto")
}
\arguments{
  \item{Data}{
The data matrix. It is assumed the 
rows are corresponding with the objects.
}
  \item{distmeasure}{
Choice of metric for the dissimilarity matrix (character).
Should be one of "tanimoto", "euclidean", "jaccard".
}
}
\details{
The tanimoto distance is calcuated as \deqn{1-\frac{N_C}{N_A+N_B-N_C}} as presented in
Li et al (2011) with N_A the number of features for object A, N_B the number of features in
object B and N_C the number of features that A and B have in common.
The Jaccard index is calcaluted as determined by the formula of the \code{dist.binary} 
function in the a4 package. First the similarity is calculated as \deqn{s=frac{n11}{n11+n10+n01}}. Then
the distance is determined as \deqn{d=\sqrt{1-s}}.
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
Dist_F=Distance(fingerprintMat,distmeasure="tanimoto")
}
\keyword{ Distances }