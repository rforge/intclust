ProfilePlotAll<-function(Genes=Genes,Comps=Comps,GeneExpr=geneMat,Raw=FALSE,Order=NULL,Color=NULL,nrclusters=7,Clusters=NULL,cols=NULL,AddLegend=TRUE,margins=c(8.1,4.1,1.1,6.5),extra=5,...){
	par(mar=margins,xpd=TRUE)
	
	if(class(GeneExpr)[1]=="ExpressionSet"){
		GeneExpr <- exprs(GeneExpr)
		
	}
	
	if(!is.null(Order)){
		if(class(Order)=="character"){
			orderlabs=Order
		}
		else{
			orderlabs=Order$Clust$order.lab
			GeneExpr=GeneExpr[,match(orderlabs,colnames(GeneExpr))]
		}
	}
	
	if(!is.null(Color)){
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
		
		if(!is.null(Order)){
			if(class(Order)=="character"){
				ordernames=Order
			}
			else{
				ordernames=Order$Clust$order.lab
			}	
			ordercolors=ordercolors[ordernames]
		}
		
		colors<- cols[ordercolors]
		names(colors) <-names(ordercolors)	
		
	}
	else{
		colors<-rep(1,ncol(GeneExpr))
	}
	
	#yvalues=c()
	#allvalues=c()
	#for(i in 1:length(Genes)){
	#	yvalues=as.vector(GeneExpr[which(rownames(GeneExpr)==Genes[i]),])
	#	allvalues=c(allvalues,yvalues-mean(yvalues))
	#}	
	yvalues=GeneExpr[Genes,]
	if(Raw==FALSE){
		allvalues=as.vector(apply(yvalues,1,function(c) c-mean(c)))
	}
	else{
		allvalues=as.vector(yvalues)
	}	
	ylims=c(min(allvalues)-0.1,max(allvalues)+0.1)
	
	plot(type="n",x=0,y=0,xlim=c(0,ncol(GeneExpr)),ylim=ylims,ylab=expression(log[2] ~ paste("fold ", "change")),xlab=" ",xaxt="n")
	#ylims=c()
	Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	
	for(i in 1:length(Genes)){
		GenesComps=as.numeric(GeneExpr[which(rownames(GeneExpr)==Genes[i]),colnames(GeneExpr)%in%Comps])
		Others=as.numeric(GeneExpr[which(rownames(GeneExpr)==Genes[i]),!(colnames(GeneExpr)%in%Comps)])
		if(length(Others)==0){
			Continue=FALSE
		}else{Continue=TRUE}
		
		#ylims=c(ylims,c(GenesComps,Others)-mean(c(GenesComps,Others)))	
		if(Raw==FALSE){
			yvalues1=GenesComps-mean(c(GenesComps,Others))	
		}
		else{
			yvalues1=GenesComps
		}
		
		lines(x=seq(1,length(GenesComps)),y=yvalues1,lty=1,col=i,lwd=1.6)
		#points(x=seq(1,length(GenesComps)),y=yvalues1,pch=19,col=i)
		segments(x0=1,y0=mean(yvalues1[1:length(GenesComps)]),x1=length(GenesComps),y1=mean(yvalues1[1:length(GenesComps)]),lwd=1.5,col=i)
		
		
		if(Continue==TRUE){
			if(Raw==FALSE){
				yvalues2=Others-mean(c(GenesComps,Others))	
			}
			else{
				yvalues2=Others
			}
			
			
			lines(x=seq(length(GenesComps)+1,ncol(GeneExpr)),y=yvalues2,lty=1,col=i,lwd=1.6)
			segments(x0=length(GenesComps)+1,y0=mean(yvalues2[1:length(Others)]),x1=ncol(GeneExpr),y1=mean(yvalues2[1:length(Others)]),lwd=1.5,col=i)
			
		}
		
		if(!(is.null(Clusters))){
			
			for(j in 1:length(Clusters)){
				temp=Clusters[[j]][!(Clusters[[j]] %in% Comps)]
				index=which(Indices%in%temp)
				
				values=as.vector(GeneExpr[which(rownames(GeneExpr)==Genes[i]),colnames(GeneExpr)%in%temp])
				values=values-mean(c(GenesComps,Others))	
				points(x=index,y=values,pch=19,col=cols[j])
				
			}
			
		}
		
	}	
	#Indices=c(colnames(GeneExpr)[which(colnames(GeneExpr)%in%Comps)],colnames(GeneExpr)[which(!colnames(GeneExpr)%in%Comps)])
	if(!is.null(Color)){
		axis(1, labels=FALSE)
		box("outer")
		mtext(substr(Indices,1,15), side = 1,  at=seq(0.5,(ncol(GeneExpr)-0.5)), line=0, las=2, cex=0.70,col=colors[Indices])
		
	}
	else{
		axis(1,at=seq(0.5,(ncol(GeneExpr)-0.5)),labels=Indices,las=2,cex.axis=0.70,xlab=" ")
	}
	axis(2,ylab=expression(log[2] ~ paste("fold ", "change")))
	
	if(AddLegend==TRUE){
		
		labels=Genes
		colslegend=seq(1,length(Genes))
		
		par(xpd=T,mar=margins)
		legend(ncol(GeneExpr)+extra,mean(c(min(ylims),max(ylims))),legend=c(labels),col=c(colslegend),lty=1,lwd=3,cex=0.8)
		
	}
}
