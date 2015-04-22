ContFeaturesPlot<-function(LeadCpds,OrderedCpds,Data,nrclusters,Color=NULL,cols=NULL,ylab="bio-assays",AddLegend=TRUE,margins=c(5.5,3.5,0.5,8.7),plottype="new",location=NULL){
	
	temp=OrderedCpds[which(!(OrderedCpds%in%LeadCpds))]
	AllCpds=c(LeadCpds,temp)
	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			dev.off()
		}
	}
	
	if(!(all(AllCpds%in%rownames(Data)))){
		Data=t(Data)
	}
	plottypein(plottype,location)
	par(mar=margins)
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
	
	mtext(LeadCpds,side=1,at=seq(1,length(LeadCpds)),line=0,las=2,cex=0.70,col=colors[LeadCpds])
	mtext(temp, side = 1, at=c(seq(length(LeadCpds)+4,(nrow(Data)+3))), line=0, las=2, cex=0.70,col=colors[temp])
	if(AddLegend==TRUE){
		
		labels=colnames(Data)
		colslegend=seq(1,length(colnames(Data)))
		
		par(xpd=T,mar=margins)
		legend(nrow(Data)+5,mean(c(min(c(min(Data)-0.5,max(Data)+0.5)),max(c(min(Data)-0.5,max(Data)+0.5)))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)
		
	}
	plottypeout(plottype)
}