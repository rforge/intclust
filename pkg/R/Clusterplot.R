ClusterPlot <-
function(Data1,Data2,nrclusters=7,cols=my_palette2,...){
	x=Data1$Clust
	Data=Data2$Clust
	
	d_temp<- dendrapply(as.dendrogram(x,hang=0.02),ClusterCols,Data,nrclusters,cols)
	
	plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2,...)
	axis(side = 2, lwd = 2)	
}
