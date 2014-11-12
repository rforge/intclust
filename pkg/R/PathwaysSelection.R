PathwaysSelection=function(List,Selection,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,GENESET=ListGO,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	
	ListNew=list()
	element=0
	for(i in 1:length(List)){
		if(class(List[[i]]) != "CEC" & class(List[[i]]) != "Weighted" & class(List[[i]]) != "WeightedSim"){
			element=element+1
			ListNew[[element]]=List[[i]]
		}
		else if(class(List[[i]])=="CEC" | class(List[[i]])=="Weighted" | class(List[[i]]) == "WeightedSim"){
			ResultsClust=list()
			if(WeightClust==TRUE){
				ResultsClust[[1]]=list()
				ResultsClust[[1]][[1]]=List[[i]]$Clust
				names(ResultsClust[[1]])[1]="Clust"
				element=element+1					
				ListNew[[element]]=ResultsClust[[1]]
				names(ListNew)[length(ListNew)]=c("WeightClust")
			}			
			else{
				for (j in 1:length(List[[i]]$Results)){
					ResultsClust[[j]]=list()
					ResultsClust[[j]][[1]]=List[[i]]$Results[[j]]
					names(ResultsClust[[j]])[1]="Clust"
					element=element+1					
					ListNew[[element]]=ResultsClust[[j]]
					names(ListNew)[length(ListNew)]=names(List[[i]]$Results)[j]
				}		
			}		
		}	
	}
	List=ListNew
	
	if(length(List)==1){
		ResultMLP=Pathways.2(List[[1]],GeneExpr,nrclusters,method,ENTREZID,geneSetSource,top,GENESET,sign)
	}
	else{
		method.test = function(sign.method,path.method){
			method.choice = FALSE
			
			if( sign.method=="limma"  & path.method=="MLP"  ){
				method.choice = TRUE
			}
			if(method.choice==TRUE){
				return(list(sign.method=sign.method,path.method=path.method))
			}	
			else{
				stop("Incorrect choice of method.")
			}
			
		}
		
		method.out = method.test(method[1],method[2])
		
		sign.method = method.out$sign.method
		path.method = method.out$path.method
		
		if(length(ENTREZID)==1){
			ENTREZID = colnames(GeneExpr)
		}
		
		# Determining the genesets if they were not given with the function input
		if((class(GENESET)=="geneSetMLP")[1] ){
			geneSet <- GENESET
		}
		else{
			geneSet <- AnnotateEntrezIDtoGO(GeneInfo[,1],database=c("ensembl","hsapiens_gene_ensembl"),attributes=c("entrezgene","go_id","description"),filters="entrezgene",species="Human")
		}
		
		if(is.null(top)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}
		
		
		if(is.null(names)){
			for(j in 1:length(List)){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		Matrix=MatrixFunction(List,nrclusters,fusionsLog,WeightClust,names)
		
		
		ResultMLP=list()
		for (k in 1:dim(Matrix)[1]){
			print(k)
			maxcluster=names(which(table(Matrix[k,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[k,which(colnames(Matrix)%in%Selection)]))))
			
			DataPrepared<-try_default(PreparePathway(List[[k]],GeneExpr,topG,sign),NULL,quiet=TRUE)
			if(is.null(DataPrepared)){
				Temp=List[[k]]
				Compounds=list()
				Compounds$LeadCpds=colnames(Matrix)[which(Matrix[k,]==maxcluster)]
				Compounds$OrderedCpds=as.hclust(Temp$Clust)$labels[as.hclust(Temp$Clust)$order]
				Temp[[length(Temp)+1]]=list(Compounds=Compounds)
				names(Temp)[length(Temp)]=paste("Cluster")
				
				DataPrepared<-PreparePathway(Temp,GeneExpr,topG,sign)
			}
			
			
			temp=list()
			temp[[1]]=DataPrepared$Compounds[[1]] #names of the compounds
			temp[[2]]=DataPrepared$Genes[[1]]
			
			
			if(path.method=="MLP"){
				## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
				p.values=DataPrepared$pvalsgenes[[1]]
				names(p.values) = ENTREZID
				
				out.mlp <- MLP(
						geneSet = geneSet,
						geneStatistic = p.values,
						minGenes = 5,
						maxGenes = 100,
						rowPermutations = TRUE,
						nPermutations = 100,
						smoothPValues = TRUE,
						probabilityVector = c(0.5, 0.9, 0.95, 0.99, 0.999, 0.9999, 0.99999),df = 9)
				
				output = list()
				#output$gene.p.values = p.adjust(p.values,method="fdr")
				
				ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)),p.values = as.numeric(out.mlp$geneSetPValue),descriptions = out.mlp$geneSetDescription)
				ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
				ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
				
				output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values<=0.05,]
				
				nr.genesets = c( dim(ranked.genesets.table)[1]  ,  length(geneSet) )
				names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
				output$nr.genesets = nr.genesets
				
				#output$object = out.mlp
				#output$method = "MLP"
				
				temp[[3]]=output				
			}
			names(temp)=c("Compounds","Genes","Pathways")			
			ResultMLP[[k]]=temp
			
		}
		names(ResultMLP)=names
		
		
	}
	
	return(ResultMLP)	
}
