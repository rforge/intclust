SharedLimma <-
function(DataLimma,names=NULL){  #Input=result of DiffGenes.2
	
	which=list()
	table=c()
	
	if(is.null(names)){
		for(j in 1:length(DataLimma)){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	nmethods=length(DataLimma)
	
	
	for (i in 1:length(DataLimma[[1]])){
		name=paste("Cluster",i,sep=" ")
		
		temp1g=c()
		comps=c()
		
		for(j in 1:nmethods){
			if(!(is.na(DataLimma[[j]][[i]])[1])){
				temp1g=c(temp1g,length(DataLimma[[j]][[i]]$Genes$ID))
				comps=c(comps,length(DataLimma[[j]][[i]]$Compounds))
			}
			else
			{
				temp1g=c(temp1g,"-")
				comps=c(comps,"-")
			}			
			names(temp1g)[j]=names[j]
			names(comps)[j]=names[j]	
		}
		
		j=1
		Continue=TRUE
		while(Continue==TRUE){
			if(!(is.na(DataLimma[[j]][[i]])[1])){
				sharedcomps=DataLimma[[j]][[i]]$Compounds
				sharedgenes=DataLimma[[j]][[i]]$Genes$ID
				
				nsharedcomps=length(DataLimma[[j]][[i]]$Compounds)
				nsharedgenes=length(DataLimma[[j]][[i]]$Genes$ID)
				names(nsharedgenes)="nshared"
				names(nsharedcomps)="nsharedcomps"
				Continue=FALSE
			}
			j=j+1
		}
		
		if(nmethods>=2){
			for (j in 2:length(DataLimma)){
				if(!(is.na(DataLimma[[j]][[i]])[1])){
					sharedcomps=intersect(sharedcomps,DataLimma[[j]][[i]]$Compounds)
					sharedgenes=intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$ID)
					
					nsharedcomps=length(intersect(sharedcomps,DataLimma[[j]][[i]]$Compounds))
					nsharedgenes=length(intersect(sharedgenes,DataLimma[[j]][[i]]$Genes$ID))
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
		
		
		which[[i]]=list(sharedcomps=sharedcomps,sharedgenes=sharedgenes,pvalsgenes=pvalsgenes)
		names(which)[i]=paste("Cluster",i,sep=" ")
		
		temp=c(temp1g,nsharedgenes,comps,nsharedcomps)			
		table=cbind(table,temp)
		colnames(table)[i]=paste("G.Cluster",i,sep=" ")
		
	}
	
	ResultShared=list(Table=table,Which=which)
	return(ResultShared)
	
}
