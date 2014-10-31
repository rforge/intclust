BoxPlotDist <-
function(Data1,Data2,lab1,lab2,limits1,limits2,type=1){ #Data input should be from cluster or integrative cluster function (such that DistM element is available).
#type: 1=plot 1, 2=plot 2, 3=both plots
	Dist1<-Data1$DistM
	Dist2<-Data2$DistM
	
	Dist1lower <- Dist1[lower.tri(Dist1)]
	Dist2lower <- Dist2[lower.tri(Dist2)]
	
	Categorize<-function(Distlower,limits){
		Cat=c(rep(0,length(Distlower)))
		for(j in 1:(length(limits)+1)){
			if(j==1){
				Cat[Distlower<=limits[j]]=j
			}
			else if(j<=length(limits)){
				Cat[Distlower>limits[j-1] & Distlower<=limits[j]]=j
			}	
			else{
				Cat[Distlower>limits[j-1]]=j
			}
		}
		Cat<-factor(Cat)
		return(Cat)
		
	}
	
	Dist1cat<-Categorize(Dist1lower,limits1)
	Dist2cat<-Categorize(Dist2lower,limits2)
	
	#plot1
	
	dataBox1<-data.frame(D1=Dist1lower,C2=Dist2cat)
	p1<-ggplot(dataBox1,aes(factor(C2),D1)) #x,y
	p1<-p1+geom_boxplot(outlier.shape=NA)+geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+
			xlab(lab2)+ylab(lab1)
	
	#plot2
	
	dataBox2<-data.frame(D2=Dist2lower,C1=Dist1cat)
	p2<-ggplot(dataBox2,aes(factor(C1),D2)) #x,y
	p2<-p2+geom_boxplot(outlier.shape=NA)+geom_point(color="blue",size=2,shape=19,position="jitter",cex=1.5)+
			xlab(lab1)+ylab(lab2)
	
	if(type==3){
		dev.new(width=14,height=7)
		grid.arrange(p1, p2, ncol=2,nrow=1)		
		#dev.new()
		#print(p1)
		#dev.new()
		#print(p2)	
	}
	else if(type==1){
		dev.new()
		print(p1)
		
	}
	else if(type==2){
		dev.new()
		print(p2)		
	}
	
}
