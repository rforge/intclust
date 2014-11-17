Cluster<-function(Data,distmeasure="tanimoto",clust="agnes",linkage="ward",gap=TRUE,maxK=50){
	Data<-Data+0
	
	if(clust != "agnes" | linkage != "ward"){
		print("Only hierarchical clustering with WARD link is implemented for now. Method continues with these options")
		clust="agnes"
		linkage="ward"
	}
	
	if(distmeasure=="tanimoto"| distmeasure == "jaccard"){
		print("If the distance measure is tanimoto or jaccard, the values are subtracted from 1 to create distances from the similarities")
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
	
	return(out)
}
