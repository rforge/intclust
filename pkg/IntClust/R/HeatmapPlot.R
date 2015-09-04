HeatmapPlot<-function(Data1,Data2,names=NULL,nrclusters=NULL,cols=NULL,plottype="new",location=NULL){
	data1=Data1$Clust
	data2=Data2$Clust
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
	DistM=.distanceheatmaps(data1,data2,names,nrclusters)
	plottypein(plottype,location)
	heatmap.2(DistM,Rowv =as.dendrogram(data1), Colv=as.dendrogram(data2),trace="none",col=cols,key=FALSE)
	plottypeout(plottype)
}