BinFeaturesPlot<-function(LeadCpds,OrderedCpds,Features,Data,Color,nrclusters=NULL,cols=NULL,name=c("FP"),margins=c(5.5,3.5,0.5,5.5),plottype="new",location=NULL){
	
	temp=OrderedCpds[which(!(OrderedCpds%in%LeadCpds))]
	AllCpds=c(LeadCpds,temp)
	
	if(all(AllCpds%in%rownames(Data))){
		Data=t(Data)
	}
	
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
	
	x<-c(1:length(AllCpds)) #x=comps
	y<-c(1:length(Features)) #y=feat
	PlotData<-t(Data[as.character(Features),AllCpds])
	plottypein(plottype,location)
	par(mar=margins)
	image(x,y,PlotData,col=c('gray90','blue'),xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	image(x[1:length(LeadCpds)],y,PlotData[1:length(LeadCpds),],col=c('gray90','green'),add=TRUE,xlab="",axes=FALSE,ann=FALSE,xaxt='n')
	
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
	
	#axis(1, labels=FALSE)
	mtext(colnames(PlotData), side = 4, at= c(1:ncol(PlotData)), line=0.2, las=2,cex=0.8)
	mtext(name, side = 2,  line=1, las=0, cex=1)
	mtext(rownames(PlotData), side = 1, at= c(1:nrow(PlotData)), line=0.2, las=2, cex=0.8,col=colors[AllCpds])
	plottypeout(plottype)
}
