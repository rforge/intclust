PathwaysIter<-function(List,Selection=NULL,GeneExpr=geneMat,nrclusters=NULL,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,niter=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	
	PathwaysOutput = list() 
	for (i in 1:niter){
		message(paste("Iteration",i,sep=" "))
		mlp = Pathways(List,Selection,GeneExpr,nrclusters,method,ENTREZID,geneSetSource,top,topG,GENESET,sign=0.05,fusionsLog,WeightClust,names)
		PathwaysOutput [[length(PathwaysOutput )+1]] = mlp
		names(PathwaysOutput )[i]=paste("Iteration",i,sep=" ")
	 }
		
	
	return(PathwaysOutput )
}
