FeatSelection<-function(List,Selection=NULL,BinData,Datanames=NULL,nrclusters=NULL,top=NULL,sign=0.05,fusionsLog=TRUE,WeightClust=TRUE,names=NULL){
	
	
	if(is.null(Datanames)){
		for(j in 1:length(BinData)){
			names[j]=paste("Data",j,sep=" ")	
		}
	}
	
	if(class(Selection)=="character"){
		
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSet <- rownames(BinData[[1]])
		
		ResultFeat=list()
		Characteristics=list()
		temp=list()
		
		LeadCpds=Selection #names of the compounds
		OrderedCpds=cpdSet
		temp[[1]]=list(LeadCpds,OrderedCpds)
		names(temp[[1]])=c("LeadCpds","OrderedCpds")
		
		group <- factor(ifelse(cpdSet %in% LeadCpds, 1, 0)) #identify the group of interest
		
		#Determine characteristic features for the compounds: fishers exact test
		result=list()
		for(j in 1: length(BinData)){
			binMat=BinData[[j]]
			
			pFish <- apply(binMat, 2, function(x) fisher.test(table(x, group))$p.value)
			
			pFish <- sort(pFish)
			adjpFish<-p.adjust(pFish, method = "fdr")
			
			AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
			
			if(is.null(top)){
				topChar=length(which(pFish<0.05))
			}
			
			TopFeat=data.frame(Names=as.character(names(pFish[0:topChar])),P.Value=pFish[0:topChar],adj.P.Val=adjpFish[0:topChar])
			temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
			result[[j]]<-temp1
			names(result)[j]=Datanames[j]
			
		}
		
		
		temp[[2]]=result
		
		names(temp)=c("Compounds","Characteristics")
		
		ResultFeat[[1]]=temp
		names(ResultFeat)="Selection"
		
		
		
	}
	
	else if(class(Selection)=="numeric" & !(is.null(List))){
		
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
		
	
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		
		for(i in 1:length(BinData)){
			BinData[[i]]=BinData[[i]]+0
			BinData[[i]]<-BinData[[i]][,which(colSums(BinData[[i]]) != 0 & colSums(BinData[[i]]) != nrow(BinData[[i]]))]
		}
		
		cpdSet <- rownames(BinData[[1]])
		
		ResultFeat=list()
		for(k in 1:dim(Matrix)[1]){	
			cluster=Selection
			
			hc<-as.hclust(List[[k]]$Clust)
			OrderedCpds <- hc$labels[hc$order]
			
			Genes=list()
			temp=list()
			LeadCpds=colnames(Matrix)[which(Matrix[k,]==cluster)] #names of the compounds
			temp[[1]]=list(LeadCpds,OrderedCpds)
			names(temp[[1]])=c("LeadCpds","OrderedCpds")
			
			group <- factor(ifelse(cpdSet %in% LeadCpds[[i]], 1, 0)) #identify the group of interest
			
			#Determine characteristic features for the compounds: fishers exact test
			result=list()
			for(j in 1: length(BinData)){
				binMat=BinData[[j]]
				
				pFish <- apply(binMat, 2, function(x) fisher.test(table(x, group))$p.value)
				
				pFish <- sort(pFish)
				adjpFish<-p.adjust(pFish, method = "fdr")
				
				AllFeat=data.frame(Names=as.character(names(pFish)),P.Value=pFish,adj.P.Val=adjpFish)
				
				if(is.null(top)){
					topChar=length(which(pFish<0.05))
				}
				
				TopFeat=data.frame(Names=as.character(names(pFish[0:topChar])),P.Value=pFish[0:topChar],adj.P.Val=adjpFish[0:topChar])
				temp1=list(TopFeat=TopFeat,AllFeat=AllFeat)
				result[[j]]<-temp1
				names(result)[j]=Datanames[j]
				
			}
			
			temp[[2]]=result
			
			names(temp)=c("Compounds","Characteristics")
			ResultFeat[[k]]=temp
			names(ResultFeat)[k]=paste(names[k],": Cluster",cluster, sep=" ")
		}		
	}
	
	else{
		message("If a specific cluster is specified, clustering results must be provided in List")
	}
	return(ResultFeat)
	
	
	
}