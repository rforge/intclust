CECb<-function(List,distmeasure=c("tanimoto","tanimoto"),nrclusters=seq(5,25,1),weight=NULL,clust="agnes",linkage="ward",WeightClust=0.5){
	
	if(class(List) != "list"){
		stop("Data must be of type list")
	}
	
	
	if(is.null(nrclusters)){
		stop("Give a number of clusters to cut the dendrogram into for each data modality.")
	}
	
	if(clust != "agnes" | linkage != "ward"){
		message("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting
						coassociation matrix.")
		clust="agnes"
		linkage="ward"
	}
	
	#Step 1: Take all features from A1 and A2
	#Notation facility:

	
	#Put up Incidence matrix for each data modality
	Incidence=list()
	for (i in 1:length(List)){
		Incidence[[i]]=matrix(0,dim(List[[i]])[1],dim(List[[i]])[1])
		rownames(Incidence[[i]])=rownames(List[[i]])
		colnames(Incidence[[i]])=rownames(List[[i]])
	}
	
	#Repeat for t iterations: not necessary here since only thing that changes is the number of clusters the tree is cut into
	
	
	#Step 2: apply hierarchical clustering on A1_prime and A2_prime + cut tree into nrclusters
	
	DistM=lapply(seq(length(List)),function(i) Distance(List[[i]],distmeasure=distmeasure[i]))
	
	HClust_A=lapply(seq(length(DistM)),function(i) agnes(DistM[[i]],diss=TRUE,method=linkage))
	
	for(k in 1:length(nrclusters)){
		message(k)
		MembersofClust=list()
		for (i in 1:length(HClust_A)){
			Temp=cutree(HClust_A[[i]],nrclusters[k])	
			MembersofClust[[i]]=matrix(1,dim(List[[i]])[1],dim(List[[i]])[1])
			
			for(l in 1:length(Temp)){
				label=Temp[l]
				sameclust=which(Temp==label)
				MembersofClust[[i]][l,sameclust]=0
			}
			Incidence[[i]]=Incidence[[i]]+MembersofClust[[i]]
		}	
	}
	
	
	if(is.null(weight)){
		weight=seq(1,0,-0.1)	
	}
	else if(class(weight)=='list' & length(weight[[1]])!=length(List)){
		stop("Give a weight for each data matrix or specify a sequence of weights")
	}
	else{
		message('The weights are considered to be a sequence, each situation is investigated')
	}
	
	if(class(weight)!="list"){
		condition<-function(l){		
			l=as.numeric(l)
			if( sum(l)==1 ){  #working with characters since with the numeric values of comb or permutations something goes not the way is should: 0.999999999<0.7+0.3<1??
				#return(row.match(l,t1))
				return(l)
			}
			else(return(0))
		}
		
		t1=permutations(n=length(weight),r=length(List),v=as.character(weight),repeats.allowed = TRUE)
		t2=lapply(seq_len(nrow(t1)), function(i) if(sum(as.numeric(t1[i,]))==1) return(as.numeric(t1[i,])) else return(0)) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
		t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
		t4=t2[which(t3!=0)]
		weight=lapply(seq(length(t4)),function(i) rev(t4[[i]]))
		
	}
	
	if(class(weight)=="list" & "x" %in% weight[[1]]){ #x indicates a free weight
		for(i in 1:length(weight)){
			w=weight[[i]]
			weightsfordata=which(w!="x") #position of the provided weight = position of the data to which the weight is given
			givenweights=as.numeric(w[weightsfordata])
			
			stilltodistribute=1-sum(givenweights)
			
			newweights=seq(stilltodistribute,0,-0.1)
			
			t1=permutations(n=length(newweights),r=length(List)-length(weightsfordata),v=as.character(newweights),repeats.allowed = TRUE)
			Input1=as.list(seq_len(nrow(t1)))
			Input2=lapply(seq(length(Input1)),function(i) {Input1[[i]][length(Input1[[i]])+1]=stilltodistribute
						return(Input1[[i]])})
			t2=lapply(seq(length(Input2)), FUN=function(i){if(sum(as.numeric(t1[Input2[[i]][1],])+0.00000000000000002775)==Input2[[i]][2]) return(as.numeric(t1[i,])) else return(0)}) #make this faster: lapply on a list or adapt permutations function itself: first perform combinations under restriction then perform permutations
			t3=sapply(seq(length(t2)),function(i) if(!all(t2[[i]]==0)) return (i) else return(0))
			weightsforotherdata=t2[which(t3!=0)]
			
			new=list()
			for(i in 1:length(weightsforotherdata)){
				w1=weightsforotherdata[[i]]
				new[[i]]=rep(0,length(List))
				new[[i]][weightsfordata]=givenweights
				new[[i]][which(new[[i]]==0)]=w1
			}
			
			weight=new
		}
	}
	
	weightedcomb<-function(w,Dist){
		temp=lapply(seq_len(length(Dist)),function(i) w[i]*Dist[[i]])
		temp=Reduce("+",temp)	
		return(temp)
	}
	IncidenceComb=lapply(weight,weightedcomb,Incidence)
	
	
	CEC=list()
	for (i in 1:length(IncidenceComb)){
		CEC[[i]]=agnes(IncidenceComb[[i]],diss=TRUE,method=linkage)
		names(CEC)[i]=paste("Weight",weight[i],sep="")
		if(all(weight[[i]]==WeightClust)){
			Clust=CEC[[i]]
		}
	}
	
	out=list(Incidence=Incidence,IncidenceComb=IncidenceComb,Results=CEC,Clust=Clust)
	attr(out,'method')<-'CEC'
	return(out)	
	
}
