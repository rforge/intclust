SNFb<-function(List,distmeasure=c("tanimoto","tanimoto"),NN=20,alpha=0.5,T=20,clust="agnes",linkage="ward"){
	
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
	
	AffinityMatrix.2=function (Diff, K = 20, sigma = 0.5) 
	{
		N = nrow(Diff)
		#Diff = (Diff + t(Diff))/2 #Why is this here?
		diag(Diff) = 0
		sortedColumns = as.matrix(t(apply(Diff, 2, sort)))
		finiteMean <- function(x) {
			mean(x[is.finite(x)])
		}
		means = apply(sortedColumns[, 1:K + 1], 1, finiteMean) + 
				.Machine$double.eps
		avg <- function(x, y) ((x + y)/2)
		Sig = outer(means, means, avg)/3 * 2 + Diff/3 + .Machine$double.eps
		Sig[Sig <= .Machine$double.eps] = .Machine$double.eps
		densities = dnorm(Diff, 0, sigma * Sig, log = FALSE)
		densities=densities*(0.5*Sig*sqrt(2*pi))  #Rescale them back to have 1 on the diagonal
		W = (densities + t(densities))/2 #Why is this here?
		W=densities
		return(W)
	}
	
	AffM=vector("list",length(List))
	for(i in 1:length(List)){
		AffM[[i]]=AffinityMatrix.2(DistM[[i]], NN, alpha)
	}
	
	#STEP 3: Fuse Networks Into 1 Single Network
	
	snf.2=function (Wall, K = 20, t = 20) 
	{
		LW = length(Wall)
		
		normalize <- function(X) {
			NMatrix=matrix(0,dim(X)[1],dim(X)[2])
			for(i in 1:dim(X)[1]){
				row=X[i,]
				row[i]=0
				D=sum(row)
				for(j in 1:dim(X)[2]){
					N=X[i,j]
					
					NMatrix[i,j]=N/(2*D)
					
					if(i==j){
						NMatrix[i,j]=1/2
					}
				}	
				
			}
			return(NMatrix)
		}	
		
		PMatrix<-vector("list",LW)
		SMatrix <- vector("list", LW)
		nextW <- vector("list", LW)
		
		
		for (i in 1:LW) {
			PMatrix[[i]] = normalize(Wall[[i]])
			Wall[[i]] = (Wall[[i]] + t(Wall[[i]]))/2
		}
		
		for (i in 1:LW) {
			
			zero <- function(x) {  #After affinityMatrix: the closest obs have the highest values.
				s = sort(x, index.return = TRUE)
				x[s$ix[1:(length(x) - K)]] = 0
				return(x)
			}
			
			
			SMatrix[[i]] = matrix(0,dim(Wall[[i]])[1],dim(Wall[[i]])[2])
			for (k in 1:nrow(SMatrix[[i]])) {
				SMatrix[[i]][k, ] = zero(Wall[[i]][k, ])
			}	
			SMatrix[[i]]=normalize(SMatrix[[i]])*2
			SMatrix[[i]]=SMatrix[[i]]-diag(nrow(SMatrix[[i]])) #Diagonal elements should be put to 0 or 1??
			
		}	
		
		for (i in 1:t) {
			
			for (j in 1:LW) {
				sumWJ = matrix(0, dim(PMatrix[[j]])[1], dim(PMatrix[[j]])[2])
				for (k in 1:LW) {
					if (k != j) {
						sumWJ = sumWJ + PMatrix[[k]]
					}
				}
				nextW[[j]] = SMatrix[[j]] %*% (sumWJ/(LW - 1)) %*% t(SMatrix[[j]]) #update PMatrix
			}
			for (j in 1:LW) {
				PMatrix[[j]] = nextW[[j]] #+ diag(nrow(Wall[[j]]))  #why +1 on diagonals?
				Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2           #why is this necessary?
				PMatrix[[j]]=normalize(PMatrix[[j]])   #normalize after every iteration
			}
		}
		
		W = matrix(0, nrow(PMatrix[[1]]), ncol(PMatrix[[1]]))
		for (i in 1:LW) {
			W = W + PMatrix[[i]]
		}
		W = W/LW
		W = normalize(W) #again normalization?
		W = (W + t(W) + diag(nrow(W)))/2
		return(W)
	}
	
	SNF_FusedM=snf.2(AffM, NN, T)
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
