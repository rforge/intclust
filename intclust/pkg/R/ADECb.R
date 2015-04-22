ADECb<-function(List,distmeasure="tanimoto",nrclusters=seq(5,25,1),clust="agnes",linkage="ward"){
	
	if(class(List) != "list"){
		stop("Data must be of type lists")
	}

	if(is.null(nrclusters)){
		stop("Give a number of cluters to cut the dendrogram into.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						coassociation matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	#Fuse A1 and A2 into 1 Data Matrix
	
	AllData<-NULL
	for (i in 1:length(List)){
		if(i==1){
			AllData=List[[1]]
		}
		else{
			AllData=cbind(AllData,List[[i]])
		}
	}
	
	
	#Put up Incidence matrix
	Incidence=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
	rownames(Incidence)=rownames(AllData)
	colnames(Incidence)=rownames(AllData)
	
	
	#Repeat for t iterations: not necessary in version b since no random number of features taken.
	
	
	#Step 2: apply hierarchical clustering on A1_prime and A2_prime + cut tree into nrclusters
	
	DistM=Distance(AllData,distmeasure)
	
	HClust_A=agnes(DistM,diss=TRUE,method=linkage)
	
	for(k in 1:length(nrclusters)){
		message(k)
		Temp=cutree(HClust_A,nrclusters[k])	
		MembersofClust=matrix(1,dim(List[[1]])[1],dim(List[[1]])[1])
		
		for(l in 1:length(Temp)){
			label=Temp[l]
			sameclust=which(Temp==label)
			MembersofClust[l,sameclust]=0
		}
		Incidence=Incidence+MembersofClust
	}	
	
	Clust=agnes(Incidence,diss=TRUE,method=linkage)
	
	out=list(AllData=AllData,S=Incidence,Clust=Clust)
	attr(out,'method')<-'ADEC'
	return(out)
	
}
