Shared <-
function(DataLimma=NULL,DataMLP=NULL,names=NULL){  #Input=result of DiffGenes.2 and Geneset.intersect
	
	if(is.null(DataMLP)){
		ResultShared=SharedLimma(DataLimma,names)
		
	}
	else if(is.null(DataLimma)){
		ResultShared=SharedMLP(DataMLP,names)		
	}
	
	else if(is.null(DataLimma) & is.null(DataMLP)){
		
		stop("At least one Data set should be specified")
	}
	
	else{ 
		
		which=list()	
		table=c()
		
		if(length(DataLimma) != length(DataMLP)){
			stop("Unequal number of methods for limma and MLP")
		}
		else{
			nmethods=length(DataLimma)
		}
		
		if(is.null(names)){
			for(j in 1:length(DataLimma)){
				names[j]=paste("Method",j,sep=" ")	
			}
		}
		
		for (i in 1:length(DataLimma[[1]])){
			name=paste("Cluster",i,sep=" ")
			
			temp1g=c()
			temp1p=c()
			comps=c()
			
			pvalsg=c()
			pvalsp=c()
			print(name)
			
			for(j in 1:nmethods){
				if(!(is.na(DataLimma[[j]][[i]])[1]) | !(is.na(DataMLP[[j]][[i]])[1])){
					temp1g=c(temp1g,length(DataLimma[[j]][[i]]$Genes$ID))
					temp1p=c(temp1p,length(DataMLP[[j]][[i]][[3]]$descriptions))
					comps=c(comps,length(DataLimma[[j]][[i]]$Compounds))
				}
				else{
					temp1g=c(temp1g,"-")
					temp1p=c(temp1p,"-")
					comps=c(comps,"-")
				}
				
				names(temp1g)[j]=names[j]
				names(temp1p)[j]=names[j]
				names(comps)[j]=paste("Ncomps", names[j],sep=" ")
				
			}
			
			j=1		
			Continue=TRUE
			while (Continue==TRUE){
				if(!(is.na(DataLimma[[j]][[i]])[1]) | !(is.na(DataMLP[[j]][[i]])[1])){
					sharedcomps=DataLimma[[j]][[i]]$Compounds
					sharedgenes=DataLimma[[j]][[i]]$Genes$ID
					sharedpaths=DataMLP[[j]][[i]][[3]]$descriptions
					
					pvalsg=c(pvalsg,DataLimma[[j]][[i]]$Genes$adj.P.Val)
					pvalsp=c(pvalsp,DataMLP[[j]][[i]]$mean_p.value)
					
					nsharedcomps=length(DataLimma[[j]][[i]]$Compounds)
					nsharedgenes=length(DataLimma[[j]][[i]]$Genes$ID)
					nsharedpaths=length(DataMLP[[j]][[i]][[3]]$descriptions)
					names(nsharedgenes)="nshared"
					names(nsharedpaths)="nshared"
					names(nsharedcomps)="nsharedcomps"
					
					Continue=FALSE
				}
				j=j+1
			}
			
			if(nmethods>=2){
				for (j in 2:length(DataLimma)){
					if(!(is.na(DataLimma[[j]][[i]])[1]) | !(is.na(DataMLP[[j]][[i]])[1])){
						sharedcomps=intersect(sharedcomps,DataLimma[[j]][[i]]$Compounds)
						sharedgenes=intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$ID)
						sharedpaths=intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$descriptions)
						
						nsharedcomps=length(intersect(sharedcomps,DataLimma[[j]][[i]]$Compounds))
						nsharedgenes=length(intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$ID))
						nsharedpaths=length(intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$descriptions))
						names(nsharedgenes)="nshared"
						names(nsharedpaths)="nshared"
						names(nsharedcomps)="nsharedcomps"
					}
				}			
			}
			
			pvalsgenes=list()
			meanpvalsgenes=c()
			meanpvalspaths=c()
			pvalspaths=list()
			
			if(nsharedgenes != 0){
				for(c in 1:nmethods){
					pvalsg=c()
					for(g in sharedgenes){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalsg=c(pvalsg,DataLimma[[c]][[i]]$Genes$adj.P.Val[DataLimma[[c]][[i]]$Genes$ID==g])	
						}	
					}
					
					pvalsgenes[[c]]=pvalsg
					#names(pvalsgenes)[c]=paste("Method",c,sep=" ")
				}	
				
				for(g1 in 1:length(sharedgenes)){
					pvalstemp=c()			
					for(c in 1:nmethods){
						if(!(is.na(DataLimma[[c]][[i]])[1])){
							pvalstemp=c(pvalstemp,pvalsgenes[[c]][[g1]])
						}
					}			
					meanpvalsgenes=c(meanpvalsgenes,mean(pvalstemp))			
				}
				pvalsgenes[[nmethods+1]]=meanpvalsgenes	
				names(pvalsgenes)[nmethods+1]="Mean pvals genes"
			}
			else{pvalsgenes=0}
			
			if(nsharedpaths!=0){
				for(c in 1:nmethods){
					pvalsp=c()
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						for(p in sharedpaths){
							pvalsp=c(pvalsp,DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$descriptions==p,3])
						}
					}
					
					pvalspaths[[c]]=pvalsp
					
					
					#names(pvalspaths)[c]=paste("Method",c,sep=" ")
				}
				
				
				for(p1 in 1:length(sharedpaths)){
					pvalstemp1=c()
					for(c in 1:nmethods){
						if(!(is.na(DataMLP[[c]][[i]])[1])){
							pvalstemp1=c(pvalstemp1,pvalspaths[[c]][[p1]])
							
						}
						
					}			
					
					meanpvalspaths=c(meanpvalspaths,mean(pvalstemp1))
					
				}
				pvalspaths[[nmethods+1]]=meanpvalspaths	
				names(pvalspaths)[nmethods+1]="Mean pvals paths"
			}
			else{pvalpaths=0}
			
			temp=rbind(cbind(temp1g,temp1p),cbind(nsharedgenes,nsharedpaths),cbind(comps,comps),cbind(nsharedcomps,nsharedcomps))	
			
			table=cbind(table,temp)
			
			#rownames(table)=names
			
			which[[i]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes,sharedpaths=sharedpaths,pvalspaths=pvalspaths)
			names(which)[i]=paste("Cluster",i,sep=" ")
			
		}
		for (i in 1:length(seq(1,dim(table)[2],2))){
			number=seq(1,dim(table)[2],2)[i]
			colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
			colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")	
			
		}
		ResultShared=list(Table=table,Which=which)
	}
	return(ResultShared)
	
}
