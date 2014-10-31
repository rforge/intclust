labelCol <-
function(x,Sel) {
	if (is.leaf(x)) {
		## fetch label
		label <- attr(x, "label") 
		## set label color to red for SelF, to black otherwise
		attr(x, "nodePar") <- list(pch=NA,lab.col=ifelse(label %in% Sel, "red", "black"),lab.cex=0.9,font=2)
		attr(x, "edgePar") <- list(lwd=2,col=ifelse(label %in% Sel, "red", "black"));
	}
	return(x)
}
