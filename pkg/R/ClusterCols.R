ClusterCols <- function(x,Data,nrclusters=7,cols=my_palette2) {
	if(length(cols)<nrclusters){
		stop("Not for every cluster a color is specified")
	}	
	
	Clustdata=cutree(Data,nrclusters)
	Clustdata=Clustdata[Data$order]
	
	ordercolors=Clustdata
	order=seq(1,nrclusters)
	
	for (k in 1:length(unique(Clustdata))){
		select=which(Clustdata==unique(Clustdata)[k])
		ordercolors[select]=order[k]
	}
	
	colfunc=function(x,cols){
		indextemp=which(attr(Data$diss,"Labels")==x)
		index1=which(Data$order==indextemp)	
		
		index2=ordercolors[index1]
		
		color=cols[index2]
		return(color)	
	}
	
	if (is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to clustercolor
		attr(x, "nodePar") <- list(pch=NA,lab.col=colfunc(label,cols),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=colfunc(label,cols))
	}
	return(x)
}
