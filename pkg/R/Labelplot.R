LabelPlot <-
function(Data,Sel1,Sel2=NULL,col1,col2=NULL,...){
	x=Data$Clust
	
	d_temp <- dendrapply(as.dendrogram(x,hang=0.02),LabelCols,Sel1,Sel2,col1,col2)
	
	plot(d_temp,nodePar=list(pch=NA),edgePar=list(lwd=2),ylab="Height",font.axis=2,font.lab=2,font=2,...)
	axis(side = 2, lwd = 2)
}
