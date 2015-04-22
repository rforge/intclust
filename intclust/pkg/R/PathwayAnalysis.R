PathwayAnalysis<-function(List,Selection=NULL,GeneExpr=geneMat,nrclusters=NULL,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,niter=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	
	Pathways=PathwaysIter(List,Selection,GeneExpr,nrclusters,method,ENTREZID,geneSetSource,top,topG,GENESET,sign,niter,fusionsLog,WeightClust,names)
	
	Intersection=Geneset.intersect(Pathways,Selection,sign,names,seperatetables,separatepvals)
	
	return(Intersection)
	
}