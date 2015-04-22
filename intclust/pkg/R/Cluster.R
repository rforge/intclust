Cluster<-function(Data,distmeasure="tanimoto",clust="agnes",linkage="ward",gap=TRUE,maxK=50){	
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented for now. Method continues with these options")
		clust="agnes"
		linkage="ward"
	}
		
	#STEP 1: Distance Matrices

	DistM=Distance(Data,distmeasure)
	rownames(DistM)=rownames(Data)
	colnames(DistM)=rownames(Data)	
	
	#STEP 2: Hierarchical Clustering with Ward Link
	
	Clust=agnes(DistM,diss=TRUE,method=linkage)
	
	func = function(x,k){
		return(list(cluster=cutree(Clust,k=k)  )  )
	}
	
	#with optional gap statistic
	if(gap==TRUE){
		Clust_gap = clusGap(Data,FUNcluster=func,K.max=maxK,B=500)
		gapdata = as.data.frame(Clust_gap$Tab)
		
		k1 = maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstSEmax")
		k2 = maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalSEmax")
		k3 = maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"firstmax")
		k4 = maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"globalmax")
		k5 = maxSE(gapdata[-maxK,3],gapdata[-maxK,4],"Tibs2001SEmax")
		
		k = data.frame(firstSEmax=k1,globalSEmax=k2,firstmax=k3,globalmax=k4,Tibs2001SEmax=k5)
		
		out = list(DistM=DistM,Clust=Clust,Clust_gap=Clust_gap,gapdata=gapdata,k=k)
	}
	else{
		out=list(DistM=DistM,Clust=Clust)
	}
	attr(out,'method')<-'Single Clustering'
	return(out)
}
