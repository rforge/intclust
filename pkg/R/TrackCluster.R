TrackCluster <- function(List,Selection,nrclusters=NULL,followMaxComps=FALSE,followClust=TRUE,fusionsLog=TRUE,WeightClust=TRUE,names=NULL,Plot=TRUE,Table=TRUE,CompletePlot=FALSE,cols){
	FoundClusters=list()
	ClusterDistribution.2<-function(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,WeightClust){
			
		Matrix=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
		
		if(followClust==TRUE){			
			cluster.interest=NULL
			m=1
			while(m<=dim(Matrix)[1] & is.null(cluster.interest)){
				clustersfound=unique(Matrix[m,which(colnames(Matrix)%in%Selection)])
				if(length(clustersfound)==1){
					cluster.interest=clustersfound
				}
				else{
					m=m+1
				}
			}
			if(is.null(cluster.interest)){
				message("This selection is not found to be part of a cluster. FollowMaxComps will be put to true and the function will proceed")
				followMaxComps=TRUE
				followClust=FALSE
			}
			
		}	

		for(i in 1:dim(Matrix)[1]){

			temp=list()
			temp[[1]]=Selection
			names(temp)[1]="Selection"
			
			clusternumbers=unique(Matrix[i,which(colnames(Matrix)%in%Selection)])
			
			nr.clusters=length(clusternumbers)
			temp[[2]]=nr.clusters
			names(temp)[2]="nr.clusters"
			
			min.together=min(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			max.together=max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))
			
			nr.min.max.together=c(min.together,max.together)
			temp[[3]]=nr.min.max.together
			names(temp)[3]="nr.min.max.together"
			
			min.perc.together <- min.together/length(Selection) *100
			max.perc.together <- max.together/length(Selection) *100
			
			perc.min.max.together =c(min.perc.together,max.perc.together) 
			temp[[4]]=perc.min.max.together	
			names(temp)[4]="perc.min.max.together"
			
			temp[[5]]=list()
			names(temp)[5]="AllClusters"
			
			for(a in 1:length(clusternumbers)){
				temp[[5]][[a]]=list()
				names(temp[[5]])[a]=paste("Cluster",clusternumbers[a],sep=" ")
				
				temp[[5]][[a]][[1]]=clusternumbers[a]
				temp[[5]][[a]][[2]]=names(which(Matrix[i,]==clusternumbers[a])) #complete cluster
				temp[[5]][[a]][[3]]=intersect(Selection,temp[[5]][[a]][[2]]) #Objects from original selection in this cluster
				temp[[5]][[a]][[4]]=temp[[5]][[a]][[2]][which(!(temp[[5]][[a]][[2]] %in% Selection))] #Objects extra to this cluster
				names(temp[[5]][[a]])=c("clusternumber","Complete cluster","Objects from original selection in this cluster","Objects extra to this cluster")					
			}
			
			if(followMaxComps==TRUE){
				
				maxcluster=names(which(table(Matrix[i,which(colnames(Matrix)%in%Selection)])==max(table(Matrix[i,which(colnames(Matrix)%in%Selection)]))))
				temp[[6]]=maxcluster
				names(temp)[6]="Cluster with max Objects"
				complabels=names(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(maxcluster))])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(maxcluster))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
				
			}
			
			if(followClust==TRUE){
								
				temp[[6]]=cluster.interest
				names(temp)[6]="Cluster"
				complabels=names(Matrix[i,which(colnames(Matrix)%in%Selection & Matrix[i,]==as.numeric(cluster.interest))])
				temp[[7]]=complabels
				names(temp)[7]="Complabels"
				complete.new.cluster=names(Matrix[i,which(Matrix[i,]==as.numeric(cluster.interest))])
				temp[[8]]=complete.new.cluster
				names(temp)[8]="Complete.new.cluster"
				extra.new.cluster=complete.new.cluster[which(!(complete.new.cluster %in% Selection))]
				temp[[9]]=extra.new.cluster
				names(temp)[9]="Extra.new.cluster"
				
				
			}
			FoundClusters[[i]]=temp	
		}
		
		FoundClusters[[length(FoundClusters)+1]]=Matrix
		return(FoundClusters)
	}	
	
	Found=ClusterDistribution.2(List,Selection,nrclusters,followMaxComps,followClust,fusionsLog,WeightClust)
	
	Matrix=Found[[length(Found)]]
	Found=Found[-length(Found)]
	
	if(is.null(names)){
		for(j in 1:dim(Matrix)[1]){
			names[j]=paste("Method",j,1)
		}
	}
	
	if(Plot==TRUE){
	
		if(followMaxComps==TRUE){
			lab1=c("Maximum of compounds of original cluster together")
			labelcluster=c()
			for(z in 1:length(Found)){
				labelcluster=c(labelcluster,Found[[z]]$Cluster)				
			}
		}
		else{lab1=c("Number of compounds still in original cluster")}
		
		nrcluster=c()
		nrcomps=c()
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
		
		plot(type="n",x=0,y=0,xlim=c(0,length(Found)+0.5),ylim=c(0,max(nrcluster,nrcomps)+2.0),xlab="",ylab="",xaxt="n",yaxt="n",cex.lab=1.25)
		lines(x=seq(1,length(Found)),y=nrcluster,lty=1,col="red",lwd=1.5)
		points(x=seq(1,length(Found)),y=nrcluster,pch=19,col="red",cex=1.5)
		
		lines(x=c(0,1),y=c(length(Selection),nrcomps[1]),lty=1,col="black",lwd=1.5)
		points(x=0,y=length(Selection),pch=19,col="black",cex=1.5)
		
		lines(x=seq(1,length(Found)),y=nrcomps,lty=1,col="blue",lwd=1.5)
		points(x=seq(1,length(Found)),y=nrcomps,pch=19,col="blue",cex=1.5)
		
		text(0,length(Selection), "S",cex=1.5,pos=1,col="black",font=2)	
		text(seq(1,length(Found)),nrcomps, labelcluster,cex=1.5,pos=1,col="black",font=2)
			
		axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1)
		axis(2,at=seq(0,max(nrcluster,nrcomps)+2.0),labels=seq(0,max(nrcluster,nrcomps)+2.0),cex=1)
		

		legend(0,max(nrcluster,nrcomps)+2.4,lty=c(1,1,0),pch=c(19,19,0),col=c("blue","red","black"),legend=c(lab1,"Number of clusters original cluster divided amongst","Cluster number"),bty="n",cex=1.2)
		
	}
	
	if(CompletePlot==TRUE){
		dev.new()
		nrcluster=c(1)
		nrcomps=c(length(Selection))
		for(j in 1:length(Found)){
			nrcluster=c(nrcluster,Found[[j]]$nr.clusters)
			nrcomps=c(nrcomps,length(Found[[j]]$Complabels))
		}
		
		plot(type="n",x=0,y=0,xlim=c(0,length(Found)+0.5),ylim=c(0,max(nrcluster,nrcomps)+1.5),xlab="",ylab="Number of Compounds",xaxt="n",yaxt="n")
		
		xnext=c()
		ynext=c()
		colorsp=c()
		for(m in 0:length(Found)){
			if(m==0){
				howmany=1
			}
			else{
				howmany=length(Found[[m]]$AllClusters)
			}

			if(m==0){
				xnext=0.1
				ynext=c(length(Selection))
				colorsp=c("black")
				
					
				points(x=xnext,y=ynext,col=colorsp,pch=19,cex=1.25)
				labelcluster="S"
				position=3
				if(!(is.integer(ynext[length(ynext)]))){
						position=1
				}
				text(xnext,ynext,labelcluster,cex=1.5,pos=position,col="black",font=2)					
				
				L1=list()
				for(n in 1:howmany){
					L1[[n]]=Selection				
				}
				
			}
			else{
				xprev=xnext
				yprev=ynext
				ynext=c()
				colorsp=c()
				xnext=rep(seq(1,length(Found))[m],howmany)
				for(n in 1:howmany){										
					ynext=c(ynext,length(Found[[m]]$AllClusters[[n]][[3]]))
					
					colorsp=c(colorsp,cols[Found[[m]]$AllClusters[[n]]$clusternumber])
					
					if(length(ynext)>1){
						for(t in 1:(length(ynext)-1)){
							if(ynext[t]==ynext[length(ynext)]){
								ynext[length(ynext)]=ynext[length(ynext)]-0.3
							}
						}	
					}
					
					points(x=xnext[n],y=ynext[length(ynext)],col=colorsp[length(colorsp)],pch=19,cex=1.25)
					labelcluster=Found[[m]]$AllClusters[[n]]$clusternumber
					position=3
					if(!(is.integer(ynext[length(ynext)]))){
						position=1
					}
					text(xnext[n],ynext[length(ynext)],labelcluster,cex=1.5,pos=position,col="black",font=2)					
				}		
				
				
				L2=L1				
				L1=list()
				for(n in 1:howmany){
					L1[[n]]=Found[[m]]$AllClusters[[n]][[3]]					
				}
				
				for(q in 1:length(L1)){
					for(p in 1:length(L2)){
						if(length(which(L2[[p]] %in% L1[[q]])) != 0){
							segments(x0=xprev[p],y0=yprev[p],x1=xnext[q],y1=ynext[q],col=colorsp[q],lwd=2)
						}							
					}
				}
			}
		}
		axis(1,at=seq(0,length(Found)),labels=c("Selection",names),las=2,cex=1.5)
		axis(2,at=seq(0,max(nrcluster,nrcomps)),labels=seq(0,max(nrcluster,nrcomps)),cex=1.5)
		legend(0,max(nrcluster,nrcomps)+1.6,pch=c(0),col=c("black"),legend=c("Cluster number"),bty="n",cex=1.2)
		
	}
	
	if(Table==TRUE & Plot==TRUE){
		SharedComps=Selection
		Extra=list()
		temp=c()
		for(a in 1:length(Found)){
			SharedComps=intersect(SharedComps,Found[[a]]$Complabels)
		}
		for(a in 1:length(Found)){
			Extra[[a]]=Found[[a]]$Complete.new.cluster[which(!(Found[[a]]$Complete.new.cluster%in%SharedComps))]
			names(Extra)[a]=names[a]
			if(followMaxComps==TRUE){
				names(Extra)[a]=paste(names(Extra)[a],"_",labelcluster[a],sep="")
			}
			temp=c(temp,length(Extra[[a]]))
		}
		
		ExtraOr=Selection[which(!(Selection%in%SharedComps))]
		
		collength=max(length(SharedComps),length(ExtraOr),temp)
		
		if(length(SharedComps)<collength){
			spaces=collength-length(SharedComps)
			SharedComps=c(SharedComps,rep(" ",spaces))
		}
		
		if(length(ExtraOr)<collength){
			spaces=collength-length(ExtraOr)
			ExtraOr=c(ExtraOr,rep(" ",spaces))
		}
		
		for(b in 1:length(Extra)){
			if(length(Extra[[b]])<collength){
				spaces=collength-length(Extra[[b]])
				Extra[[b]]=c(Extra[[b]],rep(" ",spaces))
			}
			
		}		
		table=data.frame(ExtraOr=ExtraOr,SharedComps=SharedComps)
		for(b in 2:length(Extra)){
			table[1+b]=Extra[[b]]
			colnames(table)[1+b]=names(Extra)[b]
		}
		
	}
	
	names(Found)=names
	Found[[length(Found)+1]]=table
	names(Found)[length(Found)]="Table"
	return(Found)	
}
