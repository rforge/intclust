BioassayPlot<-function(LeadCpds,OrderedCpds,Data,nrclusters,Color=NULL,cols=Colors2,ylab="bio-assays"){
	temp=OrderedCpds[which(!(OrderedCpds%in%LeadCpds))]
	AllCpds=c(LeadCpds,temp)
	
	if(!(all(AllCpds%in%rownames(Data)))){
		Data=t(Data)
	}
	
	plot(x=0,y=0,xlim=c(0,(nrow(Data)+3)),ylim=c(min(Data)-0.5,max(Data)+0.5),type="n",ylab=ylab,xlab='',xaxt='n')
	for(i in c(1:ncol(Data))){	
		lines(x=seq(1,length(LeadCpds)),y=Data[which(rownames(Data)%in%LeadCpds),i],col=i)
	}
	for(i in c(1:ncol(Data))){	
		lines(x=seq(length(LeadCpds)+4,(nrow(Data)+3)),y=Data[which(!(rownames(Data)%in%LeadCpds)),i],col=i)
	}

	if(!(is.null(Color))){
		Data1 <- Color$Clust
		ClustData1=cutree(Data1,nrclusters) 
	
		ordercolors=ClustData1[Data1$order]
		names(ordercolors)=Data1$order.lab
	
		ClustData1=ClustData1[Data1$order]	
	
	
		order=seq(1,nrclusters)
	
		for (k in 1:length(unique(ClustData1))){
			select=which(ClustData1==unique(ClustData1)[k])
			ordercolors[select]=order[k]
		}
	
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
	}
	else{
		colors<-rep("black",length(AllCpds))
		names(colors)<-AllCpds
	}
	axis(1, labels=FALSE)
	mtext(c(LeadCpds,"","","",temp), side = 1,  at=c(seq(1,length(LeadCpds)),seq(length(LeadCpds),(nrow(Data)+3))), line=0, las=2, cex=0.70,col=colors[AllCpds])
}