ADC<-function(List,distmeasure="tanimoto",normalize=FALSE,method=NULL,clust="agnes",linkage="ward"){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						fused matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	#Fuse variables into 1 data Matrix
	
	OrderNames=rownames(List[[1]])
	for(i in 1:length(List)){
		List[[i]]=List[[i]][OrderNames,]
	}
	
	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}
	
	#Compute Distance Matrix on AllData
	
	AllDataDist=Distance(AllData,distmeasure,normalize,method)
	
	#Perform hierarchical clustering with ward link on distance matrix
	
	HClust = agnes(AllDataDist,diss=TRUE,method=linkage)		
	
	
	out=list(AllData=AllData,DistM=AllDataDist,Clust=HClust)
	attr(out,'method')<-'ADC'	
	return(out)
	
}

