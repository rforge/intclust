CECa<-function(List,distmeasure=c("tanimoto","tanimoto"),t=10,r=NULL,nrclusters=NULL,weight=NULL,clust="agnes",linkage="ward",Clustweight=0.5){
	
	if(class(List) != "list"){
		stop("Data must be of type lists")
	}
	
	if(length(List) != 2){
		stop("Method is just implemented for 2 data modalities.")
	}
	
	if(is.null(nrclusters)){
		stop("Give a number of cluters to cut the dendrogram into for each data modality.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						coassociation matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	#Step 1: Take a random subset of features from A1 and A2
	#Notation facility:
	A1=List[[1]]
	A2=List[[2]]
	A=list(A1,A2)
	
	nc1=ncol(A1)
	nc2=ncol(A2)
	nc=c(nc1,nc2)
	evenn=function(x){if(x%%2!=0)x=x-1 else x}
	nc=lapply(nc,FUN=evenn)
	nc=unlist(nc)
	
	#Put up Incidence matrix for each data modality
	Incidence=list()
	for (i in 1:length(List)){
		Incidence[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
		rownames(Incidence[[i]])=rownames(List[[i]])
		colnames(Incidence[[i]])=rownames(List[[i]])
	}
	
	#Repeat for t iterations
	for(g in 1:t){
		
		
		#if r is not fixed: changes per iteration. Need 2 values of r since 2 data modalities.
		if(is.null(r)){
			r1=sample((nc[1]/2):(nc1-1),1)
			r2=sample((nc[2]/2):(nc2-1),1)
			r=c(r1,r2)
		}
		
		#take random sample:
		
		temp1=sample(ncol(A1),r[1],replace=FALSE)
		temp2=sample(ncol(A2),r[2],replace=FALSE)
		
		A1_prime=A1[,temp1]
		A2_prime=A2[,temp2]
		
		A_prime=list(A1_prime,A2_prime)
		
		#Step 2: apply hierarchical clustering on A1_prime and A2_prime + cut tree into nrclusters
		
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
		
		DistM=list()
		for (i in 1:length(A_prime)){
			DistM[[i]]=Distance(A_prime[[i]],distmeasure[i])
			
		}
		
		HClust_A_prime=list()
		for (i in 1:length(DistM)){
			HClust_A_prime[[i]]=agnes(DistM[[i]],diss=TRUE,method=linkage)
		}
		
		MembersofClust=list()
		for (i in 1:length(HClust_A_prime)){
			Temp=cutree(HClust_A_prime[[i]],nrclusters[i])	
			MembersofClust[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
			
			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				MembersofClust[[i]][l,sameclust]=1	
			}
		}
		
		for (i in 1:length(Incidence)){
			Incidence[[i]]=Incidence[[i]]+MembersofClust[[i]]			
		}
		
		
	}
	
	if(is.null(weight)){
		weight=seq(1,0,-0.1)	
	}
	
	IncidenceComb=list()
	for (i in 1:length(weight)){
		IncidenceComb[[i]]=weight[i]*Incidence[[1]]+(1-weight[i])*Incidence[[2]]
	}
	
	
	CEC=list()
	for (i in 1:length(IncidenceComb)){
		CEC[[i]]=agnes(IncidenceComb[[i]],diss=FALSE,method=linkage)
		names(CEC)[i]=paste("Weight",weight[i],sep="")
		if(weight[i]==Clustweight){
			Clust=CEC[[i]]
			
		}
	}
	
	out=list(Incidence=Incidence,IncidenceComb=IncidenceComb,Results=CEC,Clust=Clust)
	class(out)="CEC"
	return(out)	
}
