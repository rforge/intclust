ADECc<-function(List,distmeasure="tanimoto",t=10,r=NULL,nrclusters=seq(5,25,1),clust="agnes",linkage="ward"){
	
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
	
	
	#take random sample of features
	
	nc=ncol(AllData)
	evenn=function(x){if(x%%2!=0)x=x-1 else x}
	nc=evenn(nc)
	
	
	#Put up Incidence matrix
	Incidence=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
	rownames(Incidence)=rownames(AllData)
	colnames(Incidence)=rownames(AllData)
	
	
	#Repeat for t iterations
	
	
	for(g in 1:t){
		
		message(g)
		#if r is not fixed: changes per iteration. Need 1 value for r.
		if(is.null(r)){
			r=sample((nc/2):(ncol(AllData)-1),1)
		}
		
		#take random sample:
		
		temp1=sample(ncol(AllData),r,replace=FALSE)
		
		A_prime=AllData[,temp1]
		
		
		#Step 2: apply hierarchical clustering on A1_prime and A2_prime + cut tree into nrclusters
		
		DistM=Distance(A_prime,distmeasure)
		
		HClust_A_prime=agnes(DistM,diss=TRUE,method=linkage)
		
		for(k in 1:length(nrclusters)){
			
			Temp=cutree(HClust_A_prime,nrclusters[k])	
			MembersofClust=matrix(1,dim(List[[1]])[1],dim(List[[1]])[1])
			
			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				MembersofClust[l,sameclust]=0	
			}
			Incidence=Incidence+MembersofClust		
		}
		
	}
	
	Clust=agnes(Incidence,diss=TRUE,method=linkage)
	
	out=list(AllData=AllData,S=Incidence,Clust=Clust)
	attr(out,'method')<-'ADEC'
	return(out)
	
}
