DiffGenesSelection<-
function(List,Selection,GeneExpr=geneMat,nrclusters=7,method="limma",sign=0.05,top=NULL,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	if(method != "limma"){
		stop("Only the limma method is implemented to find differentially expressed genes")
	} 
	
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
		ResultLimma=DiffGenes.2(List[[1]],GeneExpr,nrclusters,method,sign,top)
	}
	
	else{
		if(is.null(top)){
			top1=FALSE
		}
		else{
			top1=TRUE
		}	
		
		
		Matrix=MatrixFunction(List,nrclusters,fusionsLog,WeightClust,names)
		
		ResultLimma=list()
		for(k in 1:dim(Matrix)[1]){	
			maxcluster=names(which(table(Matrix[k,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[k,which(colnames(Matrix)%in%Selection)]))))
			
			hc<-as.hclust(List[[k]]$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==maxcluster)] #names of the compounds
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			label = rep(0,dim(Matrix)[2])
			label[which(Matrix[k,]==maxcluster)] = 1
			label.factor = factor(label)
			
			GeneExpr.2=GeneExpr[,colnames(Matrix)]
			
			if(class(GeneExpr.2)[1]=="ExpressionSet"){
				GeneExpr.2$LeadCmpds<-label.factor 
				DElead <- limmaTwoLevels(GeneExpr.2,"LeadCpds")
				
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
				fit = lmFit(GeneExpr.2,design=design)
				fit = eBayes(fit)
				
				allDE=topTable(fit,n=dim(GeneExpr)[1],coef=2,adjust="fdr",resort.by = "logFC",sort.by="p")
				
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
			temp[[2]]=result
			
			names(temp)=c("Compounds","Genes")
			ResultLimma[[k]]=temp
			names(ResultLimma)[k]=paste(names[k],k,": Cluster", maxcluster, sep=" ")
		}		
	}
	
	return(ResultLimma)
	
}
