\name{BioassayPlot}
\alias{BioassayPlot}
\title{
Plot of bio-assays
}
\description{
The function \code{BioassayPlot} plots the bioassay values for the compounds.
It is possible to separate between compounds of interest and the other
compounds. This is a visualization to see how bio-assays differ from cluster 
to cluster. 
}
\usage{
BioassayPlot(LeadCpds,OrderedCpds,Data,nrclusters,Color=NULL,cols=Colors2,
ylab="bio-assays")
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
  \item{Data}{
The bio-assay data matrix.
}
  \item{nrclusters}{
The number of clusters to consider.
}
  \item{Color}{
The clustering result that determines the color of the labels of the
objects in the plot. If NULL, the labels are black.
}
  \item{cols}{
The colors for the labels of the objects.
}
  \item{ylab}{
The lable of the y-axis.
}
}
\value{
A plot in which the values of the bio-assays for the LeadCpds are separeted
from the others.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(Colors2)
Comps=c("Cpd1", "Cpd2", "Cpd3", "Cpd4", "Cpd5")

Data=matrix(sample(15, size = 50*5, replace = TRUE), nrow = 50, ncol = 5)
rownames(Data)=rownames(Data, do.NULL = FALSE, prefix = "row")
for(i in 1:50){
	rownames(Data)[i]=paste("Cpd",i,sep="")
}

BioassayPlot(LeadCpds=Comps,OrderedCpds=rownames(Data),Color=NULL,Data=Data,
nrclusters=7,cols=Colors2,ylab="bio-assays")



}
\keyword{Features Plot}
