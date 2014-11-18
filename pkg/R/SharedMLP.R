SharedMLP<-function(DataMLP,names=NULL){
	which=list()
	table=c()
	
	if(is.null(names)){
		for(j in 1:length(DataMLP)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nmethods=length(DataMLP)
	
	
	for (i in 1:length(DataMLP[[1]])){
		name=paste("Cluster",i,sep=" ")
		
		temp1g=c()
		temp1p=c()
		comps=c()
		
		for(j in 1:nmethods){
			if(!(is.na(DataMLP[[j]][[i]])[1])){
				temp1g=c(temp1g,length(DataMLP[[j]][[i]]$Genes$TopDE$Genes))
				temp1p=c(temp1p,length(DataMLP[[j]][[i]][[3]]$descriptions))
				comps=c(comps,length(DataMLP[[j]][[i]]$Compounds$LeadCpds))
			}
			else
			{
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
		while(Continue==TRUE){
			if(!(is.na(DataMLP[[j]][[i]])[1])){
				sharedcomps=DataMLP[[j]][[i]]$Compounds$LeadCpds
				sharedgenes=DataMLP[[j]][[i]]$Genes$TopDE$Genes
				sharedpaths=DataMLP[[j]][[i]][[3]]$descriptions
				
				nsharedcomps=length(DataMLP[[1]][[i]]$Compounds$LeadCpds)
				nsharedgenes=length(DataMLP[[1]][[i]]$Genes$TopDE$Genes)
				nsharedpaths=length(DataMLP[[1]][[i]][[3]]$descriptions)
				names(nsharedpaths)="nshared"
				names(nsharedgenes)="nshared"
				names(nsharedcomps)="nsharedcomps"
				
				Continue=FALSE
			}
			j=j+1
		}
		
		if(nmethods>=2){
			for (j in 2:length(DataMLP)){
				if(!(is.na(DataMLP[[j]][[i]])[1])){
					sharedcomps=intersect(sharedcomps,DataMLP[[j]][[i]]$Compounds$LeadCpds)
					sharedgenes=intersect(sharedgenes,DataMLP[[j]][[i]]$Genes$TopDE$Genes)
					sharedpaths=intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$descriptions)
					
					nsharedcomps=length(intersect(sharedcomps,DataMLP[[j]][[i]]$Compounds$LeadCpds))
					nsharedgenes=length(intersect(sharedgenes,DataMLP[[2]][[i]]$Genes$TopDE$Genes))
					nsharedpaths=length(intersect(sharedpaths,DataMLP[[j]][[i]][[3]]$descriptions))
					names(nsharedpaths)="nshared"
					names(nsharedgenes)="nshared"
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
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						pvalsg=c(pvalsg,DataMLP[[c]][[i]]$Genes$TopDE$adj.P.Val[DataMLP[[c]][[i]]$Genes$TopDE$Genes==g])	
					}
				}	
				print(pvalsg)
				pvalsgenes[[c]]=pvalsg
				names(pvalsgenes)[c]=paste("Method",c,sep=" ")
			}	
			
			for(g1 in 1:length(sharedgenes)){
				pvalstemp=c()			
				for(c in 1:nmethods){
					if((!(is.na(DataMLP[[c]][[i]])[1]))){
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
				for(p in sharedpaths){
					if(!(is.na(DataMLP[[c]][[i]])[1])){
						pvalsp=c(pvalsp,DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$descriptions==p,3])
					}
				}
				pvalspaths[[c]]=pvalsp
				names(pvalspaths)[c]=paste("Method",c,sep=" ")
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
		
		which[[i]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes,sharedpaths=sharedpaths,pvalspaths=pvalspaths)
		names(which)[i]=paste("Cluster",i,sep=" ")
		
		
		temp=rbind(cbind(temp1g,temp1p),cbind(nsharedgenes,nsharedpaths),cbind(comps,comps),cbind(nsharedcomps,nsharedcomps))	
		
		table=cbind(table,temp)
		
	}
	
	for (i in 1:length(seq(1,dim(table)[2],2))){
		number=seq(1,dim(table)[2],2)[i]
		colnames(table)[number]=paste("G.Cluster",i,sep=" ")	
		colnames(table)[number+1]=paste("P.Cluster",i,sep=" ")
		
	}	
	
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
	
	
}
