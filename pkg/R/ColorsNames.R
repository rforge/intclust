ColorsNames <-
function(MatrixColors,cols=my_palette2){
	Names=c()
	for (i in 1:dim(MatrixColors)[2]){
		for(j in 1:dim(MatrixColors)[1]){
			temp=MatrixColors[j,i]	
			Color=cols[temp]
			Names=c(Names,Color)
		}
	}
	return(Names)
}
