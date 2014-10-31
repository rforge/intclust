DiffGenes.2 <-
function(Data,GeneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,top=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 
	
	DataClust=Data$Clust
	
	clusters=cutree(DataClust,nrclusters)
	
	Genes=list()
	
	if(is.null(top)){
		top1=FALSE
	}
	else{
		top1=TRUE
	}
	
	for (i in unique(clusters)){
		
		hc<-as.hclust(DataClust)
		OrderedCpds <- hc$labels[hc$order]
		
		temp=list()
		LeadCpds=colnames(GeneExpr)[which(clusters==i)] 

		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		label = rep(0,dim(GeneExpr)[2])
		label[which(clusters==i)] = 1
		
		label.factor = factor(label)
		
		if(class(GeneExpr)[1]=="ExpressionSet"){
			GeneExpr$LeadCmpds<-label.factor		
			DElead <- limmaTwoLevels(GeneExpr,"LeadCpds")
			
			allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), resort.by = "logFC",sort.by="p")
			
			if(top1==TRUE){
				result = list(allDE[1:top,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			else if(top1==FALSE){
				top=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[1:top,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			
			
		}
		else{

			design = model.matrix(~label.factor)
			fit = lmFit(GeneExpr,design=design)
			fit = eBayes(fit)
			
			allDE=topTable(fit,coef=2,adjust="fdr",resort.by = "logFC",sort.by="p")
			
			if(top1==TRUE){
				result = list(allDE[1:top,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			else if(top1==FALSE){
				top=length(which(allDE$adj.P.Val<=sign))
				result = list(allDE[1:top,],allDE)
				names(result)=c("TopDE","AllDE")
				
			}
			
			# p.values = result$adj.P.Val
		}
		temp[[2]]=result
		
		names(temp)=c("Compounds","Genes")
		
		Genes[[i]]=temp
		names(Genes)[i]=paste("Cluster",i,sep=" ")
	}
	
	return(Genes)
}
