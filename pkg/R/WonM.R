WonM=function(List,distmeasure=c("tanimoto","tanimoto"),nrclusters=seq(5,25,1),clust="agnes",linkage="ward"){
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						fused matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	#Step 1: Distance Matrices
	Dist=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure[i]))
	
	#Step 2: perform hierarchical clustering on both distance matrices

	HClustering=lapply(seq(length(List)),function(i) agnes(Dist[[i]],diss=TRUE,method=linkage))

	
	#Step 3: cut the dendrograms into a range of K values
	
	#Give 0 to pair belonging together, give 1 to a pair not belonging together : ==> Distances created otherwise similarities.
	ClusterMembers<-function(HClust,nrclusters){
		Temp=lapply(seq(length(nrclusters)),function(i) cutree(HClust,nrclusters[i]))		
		CM=lapply(seq(length(nrclusters)),function(i) matrix(1,dim(List[[1]])[1],dim(List[[1]])[1]))
		
		clusters<-function(temp,cm){
			for(l in 1:length(temp)){
				label=temp[l]
				sameclust=which(temp==label)
				cm[l,sameclust]=0			
			}
			return(cm)
		}
		
		CM2=lapply(seq(length(nrclusters)),function(i) clusters(temp=Temp[[i]],cm=CM[[i]]))
		Consensus2=Reduce("+",CM2)
		return(Consensus2)
		
		
	}
	
	Consensus=lapply(seq(length(List)), function(i) ClusterMembers(HClustering[[i]],nrclusters))
	
	OverallConsensus=Reduce("+",Consensus)	
	OverallClustering=agnes(OverallConsensus,diss=TRUE,method=linkage)
	
	out=list(DistanceM=Dist,ClustSep=HClustering,Consensus=OverallConsensus,Clust=OverallClustering)
	attr(out,'method')<-'WonM'
	return(out)
}
