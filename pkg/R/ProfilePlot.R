ProfilePlot <-
function(Gene,Comps=Comps,GeneExpr=geneMat,Clusters=NULL,cols=NULL,AddLegend=TRUE,names=NULL,margins=c(8.1,4.1,1.1,6.5),extra,...){
	GenesComps=as.vector(GeneExpr[which(rownames(GeneExpr)==Gene),colnames(GeneExpr)%in%Comps])
	Others=as.vector(GeneExpr[which(rownames(GeneExpr)==Gene),!(colnames(GeneExpr)%in%Comps)])
	if(length(Others)==0){
		Continue=FALSE
	}else{Continue=TRUE}
	
	ylims=c(GenesComps,Others)
	
	par(mar=margins,xpd=TRUE)
	#par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
	plot(type="n",x=0,y=0,xlim=c(0,ncol(GeneExpr)),ylim=c(min(ylims),max(ylims)),ylab="Gene Expression",xlab=" ",xaxt="n",...)
	points(x=seq(1,length(GenesComps)),y=GenesComps,pch=19,col="red")
	segments(x0=1,y0=mean(ylims[1:length(GenesComps)]),x1=length(GenesComps),y1=mean(ylims[1:length(GenesComps)]),lwd=3)
	
	Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)),labels=Indices,las=2,cex.axis=0.70,xlab=" ")
	axis(2,ylim=c(min(ylims),max(ylims)),ylab="Gene Expression")
	
	if(Continue==TRUE){
		points(x=seq(length(GenesComps)+1,ncol(GeneExpr)),y=Others,pch=19,col="blue")
		segments(x0=length(GenesComps)+1,y0=mean(ylims[(length(GenesComps)+1):ncol(GeneExpr)]),x1=ncol(GeneExpr),y1=mean(ylims[(length(GenesComps)+1):ncol(GeneExpr)]),lwd=3)
		
	}
	
	if(!(is.null(Clusters))){
		
		for(i in 1:length(Clusters)){
			temp=Clusters[[i]][!(Clusters[[i]] %in% Comps)]
			index=which(Indices%in%temp)
			
			values=as.vector(GeneExpr[which(rownames(GeneExpr)==Gene),colnames(GeneExpr)%in%temp])
			
			points(x=index,y=values,pch=19,col=cols[i])
			
		}
		
	}
	
	if(AddLegend==TRUE){
		if(is.null(names)){
			stop("Provide the names of the used methods to draw a legend")
		}
		labels=c("Shared")
		colslegend=c("red")
		if(Continue==TRUE){
			colslegend=c(colslegend,"blue")
			labels=c(labels,"Others")
		}
		
		if(!(is.null(Clusters))){
			colslegend=c(colslegend,cols)
			for(j in 1:length(Clusters)){
				labels=c(labels,paste("Extra in",names[j],sep=" "))
			}
		}
		
		par(xpd=T,mar=margins)
		legend(length(ylims)+extra,mean(c(min(ylims),max(ylims))),legend=labels,col=colslegend,pch=19)
	}
}
