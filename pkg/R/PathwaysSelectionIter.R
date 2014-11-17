PathwaysSelectionIter<-function(List,Selection,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,topG=NULL,GENESET=ListGO,sign=0.05,niter=10,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	ListNew=list()
	element=0
	for(i in 1:length(List)){
		if(class(List[[i]]) != "CEC" & class(List[[i]]) != "Weighted"){
			element=element+1
			ListNew[[element]]=List[[i]]
		}
		else if(class(List[[i]])=="CEC" | class(List[[i]])=="Weighted"){
			ResultsClust=list()
			if(WeightClust==TRUE){
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]$Clust
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
			}			
			else{
				for (j in 1:length(List[[i]]$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					element=element+1					
					ListNew[[element]]=ResultsClust[[j]]
				}		
			}		
		}	
	}
	List=ListNew
	if(is.null(names)){
		for(j in 1:length(List)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	if(length(List)==1){
		List.output = list()
		for (i in 1:niter){
			print(paste("Iteration",i,sep=" "))
			comps.limma.mlp = Pathways.2(List[[1]],GeneExpr,nrclusters,method,ENTREZID,geneSetSource,top,topG,GENESET,sign=0.05,fusionsLog,WeightClust)
			List.output[[length(List.output)+1]] = comps.limma.mlp
			names(List.output)[i]=paste("Iteration",i,sep=" ")
			#List.cutpoints[[length(list.cutpoints)+1]] = sign		 
		}
	}
	
	else{
		List.output = list() 
		for (i in 1:niter){
			print(paste("Iteration",i,sep=" "))
			comps.limma.mlp = PathwaysSelection(List,Selection,GeneExpr,nrclusters,method,ENTREZID,geneSetSource,top,topG,GENESET,sign=0.05,fusionsLog,WeightClust,names)
			List.output[[length(List.output)+1]] = comps.limma.mlp
			names(List.output)[i]=paste("Iteration",i,sep=" ")
			#List.cutpoints[[length(list.cutpoints)+1]] = sign		 
		}
	}
	return(List.output)
}
