Geneset.intersectSelection <-
function(list.output,sign,names=NULL,seperatetables=FALSE,separatepvals=FALSE){
	if(is.null(names)){
		for(j in 1:length(list.output$"Iteration 1")){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	#put all of same method together:preparation of lists
	subsets=list()
	nmethods=length(list.output$"Iteration 1") 
	for(i in 1:nmethods){
		subsets[[i]]=list()
		
	}
	names(subsets)=names
	
	#put all of same method together: go through list.output
	for(j in 1:length(list.output)){
		name1=names(list.output)[j]
		for(k in 1:nmethods){
			name2=k
			subsets[[name2]][[name1]]=list.output[[name1]][[name2]]
			
		}
		
	}	
	
	#for every subset (= every method) take intersection over the interations per cluster
	Intersect=list()
	
	for(i in 1:length(subsets)){
		Method=subsets[[i]]
		Clusters=list()
		nclus=1
		for(j in 1:length(Method)){
			name3=paste("Iteration",j,sep=" ")
			Clusters[[name3]]=Method[[name3]]			
		}
		
		IntersectM=list()
		
		result.out=list()
		result.name = c()
		for(a in 1:length(Clusters)){ #per cluster
			if(a==1){
				Compounds=Clusters[[a]]$Compounds
				Genes=Clusters[[a]]$Genes
				
			}				
			cut = Clusters[[a]]$Pathways$ranked.genesets.table[  Clusters[[a]]$Pathways$ranked.genesets.table[,2]<=sign,]
			colnames(cut)[2] = paste("values.",a,sep="")
			cut = cut[,c(1,3,2)]
			#print(paste("For geneset table ",i,": ",dim(cut)[1]," of the ", dim(list.output[[i]]$ranked.genesets.table)[1]," remain."))
			result.out[[a]] = cut
			result.name = c(result.name,paste("genesettable",a,sep=""))
		}
		
		
		
		names(result.out) = result.name
		
		genesets.table.intersect = join_all(result.out,by=c("genesets","descriptions"),type="inner")
		genesets.table.intersect$mean_p.value=apply(genesets.table.intersect[,3:ncol(genesets.table.intersect)],1,mean)
		result.out$genesets.table.intersect = genesets.table.intersect
		
		
		
		if(separatepvals==FALSE){
			result.out$genesets.table.intersect=genesets.table.intersect[,c(1,2,ncol(genesets.table.intersect))]
		}
		
		#print(paste("All tables share ",dim(genesets.table.intersect)[1]," genesets in total."))
		
		if(seperatetables==FALSE){
			result.out=result.out$genesets.table.intersect
		}
		
		newresult=list(Compounds=Compounds,Genes=Genes,Pathways=result.out)
		
		
		#IntersectM[[a]]=
		#names(IntersectM)[a]=names(Clusters)[[a]]
		Intersect[[i]]=newresult	
		
	}
	names(Intersect)=names
	return(Intersect)
}
