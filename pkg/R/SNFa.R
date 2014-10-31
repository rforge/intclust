SNFa <-
function(List,distmeasure=c("tanimoto","tanimoto"),NN=20,alpha=0.5,T=20,clust="agnes",linkage="ward"){
	
	#Checking required data types and methods:
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	if(alpha<0.3 | alpha >1){
		print("Warning: alpha is recommended to be between 0.3 and 1 for the SNF method. Default is 0.5.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						fused matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	
	#STEP 1: Distance Matrices
		Distance=function(data,distmeasure){
		data <- data+0
		
		tanimoto = function(m){
			S = matrix(0,nrow=dim(m)[1],ncol=dim(m)[1])
			
			for(i in 1:dim(m)[1]){
				for(j in 1:i){
					N.A = sum(m[i,])
					N.B = sum(m[j,])
					N.C = sum(m[i,(m[i,]==m[j,])])
					
					if(N.A==0&N.B==0){
						coef = 1				
					}
					else{
						coef = N.C / (N.A+N.B-N.C)
					}
					S[i,j] = coef
					S[j,i] = coef
				}
				
			}
			D = 1 - S
			return(D)
		}
		
		# Computing the distance matrices
		
		if(distmeasure=="jaccard"){
			dist = dist.binary(data,method=1)
			dist = as.matrix(dist)
		}
		else if(distmeasure=="tanimoto"){
			dist = tanimoto(data)
			dist = as.matrix(dist)
			rownames(dist) <- rownames(data)
		}
		else if(distmeasure=="euclidean"){
			dist = daisy(data,metric="euclidean")
			dist = as.matrix(dist)
		}
		
		else{
			stop("Incorrect choice of distmeasure. Must be one of: tanimoto, jaccard or euclidean.")
		}
		
		return(dist)
	}
	
	DistM=vector("list",length(List))
	for (i in 1:length(List)){
		DistM[[i]]=Distance(List[[i]],distmeasure[i])
		
	}
	
	#STEP 2: Affinity Matrices
	
	AffM=vector("list",length(List))
	for(i in 1:length(List)){
		AffM[[i]]=affinityMatrix(DistM[[i]], NN, alpha)
	}
	
	#STEP 3: Fuse Networks Into 1 Single Network
	
	SNF_FusedM=SNF(AffM, NN, T)
	rownames(SNF_FusedM)=rownames(List[[1]])
	colnames(SNF_FusedM)=rownames(List[[1]])
	
	#STEP 4: Perform Hierarchical Clustering with WARD Link
	if(clust=="agnes" & linkage=="ward"){
		HClust = agnes(SNF_FusedM,diss=FALSE,method=linkage)
		
	}
	
	#Output= list with the fused matrix and the performed clustering
	out=list(SNF_FusedM=SNF_FusedM,Clust=HClust)
	return(out)
}
