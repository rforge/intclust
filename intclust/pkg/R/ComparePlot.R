ComparePlot<-function(List,nrclusters=NULL,cols=Colors2,fusionsLog=FALSE,WeightClust=FALSE,names=NULL,margins=c(8.1,3.1,3.1,4.1),...){
	
	MatrixColors=ReorderToReference(List,nrclusters,fusionsLog,WeightClust,names)
	
	Names=ColorsNames(MatrixColors,cols)
	
	nobs=dim(MatrixColors)[2]
	nmethods=dim(MatrixColors)[1]
	
	if(is.null(names)){
		for(j in 1:nmethods){
			names[j]=paste("Method",j,sep=" ")	
		}
	}
	
	similar=round(SimilarityMeasure(MatrixColors),2)
	
	par(mar=margins)
	color2D.matplot(MatrixColors,cellcolors=Names,show.values=FALSE,axes=FALSE,xlab="",ylab="",...)
	axis(1,at=seq(0.5,(nobs-0.5)),labels=colnames(MatrixColors),las=2,cex.axis=0.70)
	axis(2,at=seq(0.5,(nmethods-0.5)),labels=rev(names),cex.axis=0.65,las=2)
	axis(4,at=seq(0.5,(nmethods-0.5)),labels=rev(similar),cex.axis=0.65,las=2)
	
}
