\name{AnnotateEntrezIDtoGO}
\alias{AnnotateEntrezIDtoGO}
\title{
Annotate ENTREZID from genes to GO terms
}
\description{
The function \code{AnnotateEntrezIDtoGO} relies on the biomaRt package.
The function annotates ENTREZID's from genes to GO locations of which 
the database can be specified by the users. The available
databases and data sets can be seen by the function \code{listMarts()
} and \code{listDatasets("name of selected database")} respectively.
}
\usage{
AnnotateEntrezIDtoGO(entrez,database,attributes=c("entrezgene","go_id",
"description"),filters="entrezgene",species="Human")
}
\arguments{
  \item{entrez}{
The ENTREZ ID's from the genes.
}
  \item{database}{
A character vector with the first element the selected database and the
 second the data set.
}
  \item{attributes}{
A character vector of the attributes from the genes one wishes to select
 from the data set. 
Defaults to c("entrezgene","go_id","description").
}
  \item{filters}{
A character string one wishes to filter the data set on. Defaults to 
"entrezgene".
}
 \item{species}{
A character string indicating the species of the GO annotations. Defaults
to "Human". This is rather here to complete the class of the object than
it has a specific use in the function.
}
}
\details{
The function relies only on one function of the biomaRt package:\code{getBM}
and was written with the purpose of no longer needing the org.Hs.eg.db package
on which the MLP package relies since this package has been reported to fail
occasionally. The result of the \code{getBM} function is transformed to the
 correct data structure to fit in the \code{MLP} function of the MLP package.
 }
\value{
The returned output is a list with one element per GO term. Each GO term
contains the genes it was annotated to. The value has class "geneSetMLP"
to fit in the \code{MLP} function.
}
\author{
Marijke Van Moerbeke
}
\examples{
data(GeneInfo)
\dontrun{
ListGO<-AnnotateEntrezIDtoGO(GeneInfo[,1],database=c("ensembl","hsapiens_gene_ensembl"),
attributes=c("entrezgene","go_id","description"),filters="entrezgene",species="Human")
}
}
\keyword{ Annotation to GO }
\keyword{ Pathway Analysis }
