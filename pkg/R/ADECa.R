ADECa <-
function(List,distmeasure="tanimoto",t=10,r=NULL,nrclusters=NULL,clust="agnes",linkage="ward"){
	
	if(class(List) != "list"){
		stop("Data must be of type lists")
	}
	
	if(length(List) != 2){
		stop("Method is just implemented for 2 data modalities.")
	}
	
	if(is.null(nrclusters)){
		stop("Give a number of cluters to cut the dendrogram into.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
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
		print(g)
		
		#if r is not fixed: changes per iteration. Need 1 value for r.
		if(is.null(r)){
			r=sample((nc/2):(ncol(AllData)-1),1)
		}
		
		#take random sample:
		
		temp1=sample(ncol(AllData),r,replace=FALSE)
		
		A_prime=AllData[,temp1]
		
		
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
	
		DistM=Distance(A_prime,distmeasure)
		
		HClust_A_prime=agnes(DistM,diss=TRUE,method=linkage)
		
		
		Temp=cutree(HClust_A_prime,nrclusters)	
		MembersofClust=matrix(0,dim(AllData)[1],dim(AllData)[1])
		
		for(l in 1:length(Temp)){
			label=Temp[l]
			sameclust=which(Temp==label)
			MembersofClust[l,sameclust]=1	
			
		}
		
		Incidence=Incidence+MembersofClust		
	}
	
	Clust=agnes(Incidence,diss=FALSE,method=linkage)
	
	out=list(AllData=AllData,Se=Incidence,Clust=Clust)
	return(out)
	
}
