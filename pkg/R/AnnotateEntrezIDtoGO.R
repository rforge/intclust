AnnotateEntrezIDtoGO <-
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
