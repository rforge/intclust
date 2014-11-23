ChooseFeatures<-function(Interactive=TRUE,LeadCpds=NULL,ClusterResult,ClusterColors=NULL,BinData,Datanames=c("FP"),GeneExpr,topChar = 20, topG = 20,sign=0.05,nrclusters=7,cols=Colors2,N=1){
	
	OrInteractive=Interactive
	for(a in 1:length(BinData)){
		if((all(rownames(BinData[[a]])%in%ClusterResult$Clust$order.lab))){
			BinData[[a]]=t(BinData[[a]])
		}
		
	}
	if(Interactive==TRUE){
		ClusterPlot(ClusterResult,ClusterColors,nrclusters,cols)
		hc1<-as.hclust(ClusterResult$Clust)
		ClusterSpecs<-list()
		ClusterSpecs=identify(hc1, N=N, MAXCLUSTER = ncol(BinData[[1]]), function(j) ChooseFeatures(Interactive=FALSE,LeadCpds=colnames(BinData[[1]][,j]),ClusterResult,ClusterColors=NULL,BinData,Datanames,GeneExpr,topChar,topG,sign,nrclusters,cols))		
		
		names(ClusterSpecs)<-sapply(seq(1,N),FUN=function(x) paste("Choice",x,sep=" "))
		
	}
	else{
		
		#BinData can be a list of different binary matrices :: different sources of information
		if(class(BinData)!="list"){
			stop("The binary data matrices must be put into a list")
		}
		for(i in 1:length(BinData)){
			BinData[[i]]<-BinData[[i]][which(rowSums(BinData[[i]]) != 0 & rowSums(BinData[[i]]) != ncol(BinData[[i]])),]
		}
		
		cpdSet <- colnames(BinData[[1]])
		
		DistW<-ClusterResult$DistW
		Clust<-ClusterResult$Clust
		
		hc <- as.hclust(Clust)
		OrderedCpds <- hc$labels[hc$order]
		
		if(class(LeadCpds)=="character"){
			LeadCpds=list(LeadCpds)
		}
		
		Specs=list()
		for(i in 1:length(LeadCpds)){
			Compounds=list(LeadCpds[[i]],OrderedCpds)
			names(Compounds)=c("LeadCpds","OrderedCpds")
			
			group <- factor(ifelse(cpdSet %in% LeadCpds[[i]], 1, 0)) #identify the group of interest
			
			#Determine characteristic features for the compounds: fishers exact test
			Characteristics=list()
			for(j in 1: length(BinData)){
				binMat=BinData[[j]]
				
				pFish <- apply(binMat, 1, function(x) fisher.test(table(x, group))$p.value)
				
				pFish <- sort(pFish)
				
				if(is.null(topChar)){
					topChar=length(which(pFish<0.05))
				}
				Characteristics[[j]]<- names(pFish[0:topChar])
				names(Characteristics)[j]=Datanames[j]
				
			}
			
			#Determine DE Genes with limma --> make difference between "regular" data matrix and "expression set"
			
			if(class(GeneExpr)[1]=="ExpressionSet"){
				GeneExpr$LeadCmpds<-group		
				DElead <- limmaTwoLevels(GeneExpr,"LeadCpds")
				
				#allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), resort.by = "logFC",sort.by="p")
				allDE <- a4Core::topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL),sort.by="p")
				if(is.null(allDE$ID)){
					allDE$Genes <- rownames(allDE)
				}
				else
				{
					allDE$Genes=allDE$ID
				}
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE <- allDE[0:topG, ]
				#TopAdjPval<-TopDE$adj.P.Val
				#TopRawPval<-TopDE$P.Value
				
				#RawpVal<-allDE$P.Value
				#AdjpVal <- allDE$adj.P.Val
				#genesEntrezId <- allDE$ENTREZID
				
				Genes<-list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")
				#Genes <- list(TopDE$SYMBOL,TopAdjPval,TopRawPval,genesEntrezId,RawpVal,AdjpVal)	
				#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
			}
			else{
				
				label.factor = factor(group)
				design = model.matrix(~label.factor)
				fit = lmFit(GeneExpr,design=design)
				fit = eBayes(fit)
				
				#allDE = topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr),resort.by = "logFC", sort.by="p")
				allDE = limma::topTable(fit,coef=2,adjust="fdr",n=nrow(GeneExpr), sort.by="p")
				if(is.null(allDE$ID)){
					allDE$Genes <- rownames(allDE)
				}
				else
				{
					allDE$Genes=allDE$ID
				}
				if(is.null(topG)){
					topG=length(which(allDE$adj.P.Val<=sign))
				}
				TopDE=allDE[0:topG,]
				#TopAdjPval<-TopDE$adj.P.Val
				#TopRawPval<-TopDE$P.Value
				
				#RawpVal<-allDE$P.Value
				#AdjpVal <- allDE$adj.P.Val
				
				Genes<-list(TopDE,allDE)
				names(Genes)<-c("TopDE","AllDE")
				#Genes <- list(TopDE[,1],TopAdjPval,TopRawPval,allDE[,1],RawpVal,AdjpVal)	
				#names(Genes)<-c("DE_Genes","DE_RawPvals","DE_AdjPvals", "All_Genes", "All_RawPvals","All_AdjPvals")
				
			}
			
			out=list(Compounds,Characteristics,Genes)
			names(out)=c("Compounds","Characteristics","Genes")
			Specs[[i]]=out
			names(Specs)[i]=paste("Choice",i,sep=" ")
		}		
		if(OrInteractive==TRUE|length(Specs)==1){
			return(out)
		}
		else{
			return(Specs)
		}
	}
	return(ClusterSpecs)
}
