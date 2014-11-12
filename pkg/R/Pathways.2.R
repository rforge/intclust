Pathways.2<-function(Data,GeneExpr=geneMat,nrclusters=7,method=c("limma", "MLP"),ENTREZID=GeneInfo[,1],geneSetSource = "GOBP",top=NULL,GENESET=ListGO,topGsign=0.05){
	
	DataPrepared<-try_default(PreparePathway(Data,GeneExpr,topG,sign),NULL,quiet=TRUE)
	
	if(is.null(DataPrepared)){
		clusters=cutree(Data$Clust,nrclusters)
		names(clusters)=colnames(GeneExpr)
		Temp=list(Data)
		#names of clusters???
		for(i in unique(clusters)){
			Compounds=list()
			Compounds$LeadCpds=names(clusters)[which(clusters==i)] 
			Compounds$OrderedCpds=as.hclust(Data$Clust)$labels[as.hclust(Data$Clust)$order]
			Temp[[i+1]]=list(Compounds=Compounds)
			names(Temp)[i+1]=paste("Cluster",i,sep=" ")
		}
		DataPrepared<-PreparePathway(Temp,GeneExpr,topG,sign)
	}
	
	
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
	
	
	ResultMLP=list()
	
	#clusters=cutree(Data,nrclusters)
	
	for (i in 1:length(DataPrepared$pvalsgenes)){
		print(i)
		temp=list()
		temp[[1]]=DataPrepared$Compounds[[i]] #names of the compounds
		temp[[2]]=DataPrepared$Genes[[i]]				
		pvalscluster=DataPrepared$pvalsgenes[[i]]
		
		if(path.method=="MLP"){
			## WE WILL USE THE RAW P-VALUES TO PUT IN MLP -> LESS GRANULAR
			
			names(pvalscluster) = ENTREZID
			
			out.mlp <- MLP(
					geneSet = geneSet,
					geneStatistic = pvalscluster,
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
			
			nr.genesets = c( dim(ranked.genesets.table)[1],length(geneSet) )
			names(nr.genesets) = c("used.nr.genesets","total.nr.genesets")
			output$nr.genesets = nr.genesets
			
			#output$object = out.mlp
			#output$method = "MLP"
			
			temp[[3]]=output				
		}
		names(temp)=c("Compounds","Genes","Pathways")
		ResultMLP[[i]]=temp
		names(ResultMLP)[i]=paste("Cluster",i,sep=" ")
	}
	return(ResultMLP)	
}
