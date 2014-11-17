WeightedClust <- function(List,distmeasure=c("tanimoto","tanimoto"),weight=seq(1,0,-0.1),Clustweight=0.5,clust="agnes",linkage="ward"){ # weight = weight to data1
	
	#Step 1: compute distance matrices:
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
	
	for(a in 1:length(distmeasure)){
		if(distmeasure[[a]]=='euclidean'){
			stand<-function(c){
				minc=min(c)
				maxc=max(c)
				c1=(c-minc)/(maxc-minc)
				return(c1)				
			}
			List[[a]]=apply(List[[a]],2,stand)			
		}
		
	}
	DistM=vector("list",length(List))
	for (i in 1:length(List)){
		DistM[[i]]=Distance(List[[i]],distmeasure[i])
		rownames(DistM[[i]])=rownames(List[[i]])
		colnames(DistM[[i]])=rownames(List[[i]])
		
	}
	
	#Step 2: Weighted linear combination of the distance matrices:
	if(length(List) != 2){
		stop("Method implemented for 2 data matrices only.")
	}
	
	WeightedDist=list()
	for (i in 1:length(weight)){
		alpha=weight[i]
		WeightedDist[[i]]=alpha*DistM[[1]]+(1-alpha)*DistM[[2]]
	}	
	
	WeightedClust=list()
	for(i in 1:length(WeightedDist)){
		WeightedClust[[i]]=agnes(WeightedDist[[i]],diss=TRUE,method=linkage)
		names(WeightedClust)[i]=paste("Weight",weight[i],sep=" ")
		if(weight[i]==Clustweight){
			Clust=WeightedClust[[i]]			
		}
		
	}
	
	# return list with objects
	out=list(Dist=DistM,WeightedDist=WeightedDist,Results=WeightedClust,Clust=Clust)
	class(out)="Weighted"
	return(out)
	
}
