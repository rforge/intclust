PlotPathways<-function(Pathways,nRow=5,main=NULL,plottype="new",location=NULL){	
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
	#preparing data structure for the plotGOGraph
	colnames(Pathways)[3:ncol(Pathways)]=sub("mean_","",colnames(Pathways)[3:ncol(Pathways)])	
	#plot GOgraph
	plottypein(plottype,location)
	PlotGOGraph(Pathways,nRow=nRow,main=main)
	plottypeout(plottype)
	
}


