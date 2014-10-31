WonM <-
function(List,distmeasure=c("tanimoto","tanimoto"),nrclusters=seq(5,25,1),clust="agnes",linkage="ward"){
	
	#Step 1: Distance Matrices
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
		rownames(DistM[[i]])=rownames(List[[i]])
		colnames(DistM[[i]])=rownames(List[[i]])
		
	}
	
	#Step 2: perform hierarchical clustering on both distance matrices
	HClustering=list()
	
	if(clust=="agnes" & linkage=="ward"){
		for(i in 1:length(List)){
			HClustering[[i]]=agnes(DistM[[i]],diss=TRUE,method=linkage)
			
		}
		
	}
	
	#Step 3: cut the dendrograms into a range of K values
	
	Consensus=list()
	for (i in 1:length(List)){
		Consensus[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
		rownames(Consensus[[i]])=rownames(List[[i]])
		colnames(Consensus[[i]])=rownames(List[[i]])
		
	}
	
	for(k in nrclusters){
		for(j in 1:length(List)){
			Temp=cutree(HClustering[[j]],k)
			CM=matrix(0,dim(List[[j]])[1],dim(List[[j]])[1])
			
			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				CM[l,sameclust]=1				
			}	
			
			Consensus[[j]]=Consensus[[j]]+CM
			
		}
		
		
	}
	
	OverallConsensus=matrix(0,dim(List[[1]])[1],dim(List[[1]])[1])
	for (i in 1:length(Consensus)){
		OverallConsensus=OverallConsensus+Consensus[[i]]  #Maybe weighted combination?
	}
	
	OverallClustering=agnes(OverallConsensus,diss=FALSE,method=linkage)
	
	out=list(DistanceM=DistM,ClustSep=HClustering,Consensus=OverallConsensus,Clust=OverallClustering)
	return(out)
}
