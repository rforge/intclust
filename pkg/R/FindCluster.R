FindCluster<-function(List,nrclusters=NULL,select=c(1,4),fusionsLog=TRUE, WeightClust=TRUE,names=NULL){
	Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
	methodnr=select[1]
	clusternr=select[2]
	Comps=names(which(Matrix[methodnr,]==clusternr))
	return(Comps)
}	