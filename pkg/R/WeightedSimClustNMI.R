WeightedSimClustNMI <-
function(List,type=c("data","clusters"),w=seq(0,1,0.01),clust="agnes",linkage="ward",distmeasure=c("euclidean","tanimoto"),gap=FALSE,maxK=50,nrclusters=NULL,names=c("B","FP"),AllClusters=FALSE){
	if(length(w)==1){
		if(type=="data"){
			#If given data matrices.
			##1: clustering on separate sources.
			##2: extract distance matrices.
			##3: determine nrclusters if not given.
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
			
			
			Clust1 <- Cluster(List[[1]], distmeasure = distmeasure[1], clust, linkage, gap, maxK)
			Clust2 <- Cluster(List[[2]], distmeasure = distmeasure[2], clust, linkage, gap, maxK)
			
			Dist1=Clust1$DistM
			Dist2=Clust2$DistM
			
			#need: gap? ncluster?
		}
		else{
			
			Dist1=List[[1]]$DistM
			Dist2=List[[2]]$DistM
			
			if((all(range(Dist1)==c(0,1)) & (range(Dist2)[1] !=0 & range(Dist2)[2]!=1)) | ((range(Dist1)[1] !=0 & range(Dist1)[2]!=1) & all(range(Dist2)==c(0,1)))){
				print('warning: the distance matrices do not have the same range. It might be that that one data set is bainary and the other is continuous. Standardzie the variables.')
				
			}
			
		}
		Weight=w
		
	}
	else{
		
		temp=ChooseWeightNMI(List,type,w,nrclusters,distmeasure,clust,linkage,gap,maxK,names)
				
		Dist1=temp[[1]]$DistM
		Dist2=temp[[2]]$DistM
		
		Weight=temp[[4]]		
	}
	
	#Weighted Distance matrix
	
	DistW=Weight*Dist1+(1-Weight)*Dist2
	
	#Clustering on weighted distance matrix
	
	if(clust!="agnes"|linkage!="ward"){
		print("Only agglomerative hierarchical clustering with ward linkage implemented.")
	}
	
	WeightedSimCluster=agnes(DistW,diss=TRUE,method=linkage)
	
	out=list(Dist1, Dist2, Weight,DistW,Clust=WeightedSimCluster)	
	names(out)=c("Dist1","Dist2","Weight","DistW","Clust")
	
	AllCluster=list()
	if(AllClusters==TRUE){
		for(a in 1:length(w)){
			W1=w[a]*Dist1+(1-w[a])*Dist2
			AllCluster[[a]]=agnes(W1,diss=TRUE,method=linkage)
			names(AllCluster)[a]=paste("Weight_",w[a],sep="")
			
		}
		
		out=list(Dist1, Dist2, Weight,DistW,Clust=WeightedSimCluster,AllCluster)	
		names(out)=c("Dist1","Dist2","Weight","DistW","Clust","Results")
	}	
	class(out)="WeightedSim"
	return(out)
}
