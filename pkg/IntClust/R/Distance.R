Distance=function(Data,distmeasure=c("tanimoto","jaccard","euclidean","hamming","cont tanimoto"), normalize=FALSE,method=NULL){
	Data <- Data+0
	distmeasure=match.arg(distmeasure)
	if((distmeasure=="euclidean") & normalize==TRUE){
		Data=Normalization(Data,method)				
	}
	
	tanimoto = function(m){
		S = matrix(0,nrow=dim(m)[1],ncol=dim(m)[1])
		
#		for(i in 1:dim(m)[1]){
#			for(j in 1:i){
#				N.A = sum(m[i,])
#				N.B = sum(m[j,])
#				N.C = sum(m[i,(m[i,]==m[j,])])
#				
#				if(N.A==0&N.B==0){
#					coef = 1				
#				}
#				else{
#					coef = N.C / (N.A+N.B-N.C)
#				}
#				S[i,j] = coef
#				S[j,i] = coef
#			}
#			
#		}
		#via matrix multiplication
		m=as.matrix(m)
		N.C=m %*% t(m)
		N.A=m %*% (1-t(m))
		N.B=(1-m) %*% t(m)
		S=N.C/(N.A+N.B+N.C)
		D = 1 - S
		return(D)
	}
	
	# Computing the distance matrices
	
	if(distmeasure=="jaccard"){
		dist = dist.binary(Data,method=1)
		dist = as.matrix(dist)
	}
	else if(distmeasure=="tanimoto"){
		dist = tanimoto(Data)
		dist = as.matrix(dist)
		rownames(dist) <- rownames(Data)
	}
	else if(distmeasure=="euclidean"){
		dist = daisy(Data,metric="euclidean")
		dist = as.matrix(dist)
	}
	else if(distmeasure=="hamming"){
		dist=hamming.distance(Data)
		dist=as.matrix(dist)
	}
	else if(distmeasure=="cont tanimoto"){
		if(normalize==TRUE){
			Data=Normalization(Data,method)				
		}		
#		S=matrix(0,nrow=nrow(Data),ncol=nrow(Data))
#		for(i in 1:nrow(Data)){
#          for(j in 1:i){
#
#				N.AB= sum(Data[i,]*Data[j,])
#				N.A = sum(Data[i,]*Data[i,])
#				N.B = sum(Data[j,]*Data[j,])
#				
#			
#				coef = N.AB / (N.A+N.B-N.AB)
#				S[i,j] = coef
#				S[j,i] = coef
#			}
#			
#		}
		#matrix manipulation (generally faster)
		m=as.matrix(Data)
		X_A_X_B=m %*% t(m)
		temp=diag(X_A_X_B)
		X_A=matrix(rep(temp,nrow(Data)),nrow=nrow(Data),ncol=nrow(Data),byrow=FALSE)
		X_B=t(X_A)
		
		Denom=X_A+X_B-X_A_X_B
		
		S=X_A_X_B/Denom
		
		dist=1-S
		
			
	}

	else{
		stop("Incorrect choice of distmeasure. Must be one of: tanimoto, jaccard or euclidean.")
	}
	colnames(dist)=rownames(dist)
	
	return(dist)
}
