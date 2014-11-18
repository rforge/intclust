FeaturesPlot<-function(LeadCpds,OrderedCpds,Features,Data,Color,nrclusters=7,cols=Colors2,name=c("FP")){
	par(mar=c(6,3,0,6))
	
	temp=OrderedCpds[which(!(OrderedCpds%in%LeadCpds))]
	AllCpds=c(LeadCpds,temp)
	
	if(all(AllCpds%in%rownames(Data))){
		Data=t(Data)
	}
	x<-c(1:length(AllCpds))
	y<-c(1:length(Features))
	PlotData<-t(Data[Features,AllCpds])
	image(x,y,PlotData,col=c('gray90','blue'),axes=FALSE)
	image(x[1:length(LeadCpds)],y,PlotData[1:length(LeadCpds),],col=c('gray90','green'),add=TRUE)
	
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
	
	
	mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0, las=2,cex=0.8)
	mtext(name, side = 2,  line=1, las=0, cex=1)
	mtext(rownames(PlotData), side = 1, at= c(1:nrow(PlotData)), line=0, las=2, cex=0.8,col=colors[AllCpds])
	
}
