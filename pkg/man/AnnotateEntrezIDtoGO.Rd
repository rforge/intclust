\name{ADC}
\alias{ADC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate Genes to GO Locations
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
AnnotateEntrezIDtoGO(entrez,database,attributes=c("entrezgene","go_id","description"),filters="entrezgene",species="Human")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{entrez}{
%%     ~~Describe \code{List} here~~
}
  \item{database}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{attributes}{
%%     ~~Describe \code{clust} here~~
}
  \item{filters}{
%%     ~~Describe \code{linkage} here~~
}
 \item{species}{
%%     ~~Describe \code{linkage} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(entrez,database,attributes=c("entrezgene","go_id","description"),filters="entrezgene",species="Human"){
	mart=useMart(database[1])
	mart=useDataset(database[2],mart=mart)
	
	goids = getBM(attributes=attributes, filters=filters, values=entrez, mart=mart)

	transformgoid<-function(goid){
		uniqueGO<-unique(goid[,2])
		GONames<-unique(goid[,3])
		listGO<-list()
	
		collect<-function(name,goid){
			genes<-goids[which(goid[,2]==name),1]
			return(genes)
		}

		ListGO<-lapply(uniqueGO,collect,goid)
		names(ListGO)=uniqueGO
		attr(ListGO,'descriptions') <- goid[,3]
		attr(attr(ListGO,"descriptions"),'names')<-goid[,2]
		return(ListGO)
	}
	
	ListGO=transformgoid(goids)
	attr(ListGO,'species')<-species	
	attr(ListGO,'geneSetSource')<-"GODB"
	class(ListGO)<-"geneSetMLP"
	return(ListGO)
}	
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
