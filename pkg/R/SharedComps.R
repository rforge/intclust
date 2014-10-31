SharedComps <-
function(List){
	comps=List[[1]]$Compounds	
	for(i in 2:length(List)){
		comps=intersect(comps,List[[i]]$Compounds)
	}	
	return(comps)
}
