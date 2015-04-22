SNFa=function(List,distmeasure=c("tanimoto","tanimoto"),NN=20,alpha=0.5,T=20,clust="agnes",linkage="ward"){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(alpha<0.3 | alpha >1){
		message("Warning: alpha is recommended to be between 0.3 and 1 for the SNF method. Default is 0.5.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						fused matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	
	#STEP 1: Distance Matrices
	DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i]))
	
	#STEP 2: Affinity Matrices
	
	AffM=lapply(seq(length(List)), function(i) affinityMatrix(DistM[[i]], NN, alpha))
	
	#STEP 3: Fuse Networks Into 1 Single Network
	
	SNF_FusedM=SNF(AffM, NN, T)
	rownames(SNF_FusedM)=rownames(List[[1]])
	colnames(SNF_FusedM)=rownames(List[[1]])
	Dist=1-SNF_FusedM
	
	#STEP 4: Perform Hierarchical Clustering with WARD Link
	if(clust=="agnes" & linkage=="ward"){
		HClust = agnes(Dist,diss=TRUE,method=linkage)		
	}
	
	#Output= list with the fused matrix and the performed clustering
	out=list(SNF_FusedM=SNF_FusedM,DistM=Dist,Clust=HClust)
	attr(out,'method')<-'SNF'
	return(out)
}
