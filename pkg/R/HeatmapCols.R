HeatmapCols<-function(Data1,Data2,names=rownames(fingerprintMat),nrclusters=7,cols=my_palette){
	data1=Data1$Clust
	data2=Data2$Clust
	
	DistM=.distanceheatmaps(data1,data2,names,nrclusters)
	
	heatmap.2(DistM,Rowv =as.dendrogram(data1), Colv=as.dendrogram(data2),trace="none",col=cols,key=FALSE)
}