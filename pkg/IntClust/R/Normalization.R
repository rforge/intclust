Normalization<-function(Data,method=c("Quantile","Fisher-Yates","Standardize","Range","Q","q","F","f","S","s","R","r")){
		
		method=substring(method[1],1,1)  #Function also accepts begin letter of each method	
		method=match.arg(method)
		if(method=="S"|method=="s"){
			Data1<-stdize(Data)#center and divide by standard deviation
			return(Data1)
		}
		
		else if(method=="R"|method=="r"){
			rangenorm<-function(x){
				minc=min(x)
				maxc=max(x)
				rx=(x-minc)/(maxc-minc)
			}
			
			DataN=apply(Data,2,rangenorm)
			Data=DataN
			return(Data)

		}
		else if(method=="Q"|method=="q"){
			DataColSorted <- apply(Data, 2, sort)
			NValues=apply(DataColSorted,1,median,na.rm=T)
		}
		else{
			DataColSorted <- apply(Data, 2, sort)
			NValues=qnorm(1:nrow(Data)/(nrow(Data)+1))
		}
		
	DataRanked=c(apply(Data,2,rank))
		
	array(approx(1:nrow(Data),NValues,DataRanked)$y,dim(Data),dimnames(Data))
	t=approx(1:nrow(Data),NValues,DataRanked)
	return(Data)
	
}

