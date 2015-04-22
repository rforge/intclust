PlotPathways<-function(Pathways,nRow=5,main=NULL,plottype="new",location=NULL){	
	plottypein<-function(plottype,location){
		if(plottype=="pdf" & !(is.null(location))){
			pdf(paste(location,".pdf",sep=""))
		}
		if(plottype=="new"){
			dev.new()
		}
		if(plottype=="sweave"){
			
		}
	}
	plottypeout<-function(plottype){
		if(plottype=="pdf"){
			dev.off()
		}
	}
	#preparing data structure for the plotGOGraph
	colnames(Pathways)[3:ncol(Pathways)]=sub("mean_","",colnames(Pathways)[3:ncol(Pathways)])	
	#plot GOgraph
	plottypein(plottype,location)
	PlotGOGraph_AdjustLegend(Pathways,nRow=nRow,main=main)
	plottypeout(plottype)
	
}



PlotGOGraph_AdjustLegend<-function (object, nRow = 5, main = NULL) 
{
	if (!inherits(object, "MLP")) 
		stop("The 'object' argument should be an object of class 'MLP' as produced by the MLP function")
	if (is.data.frame(attributes(object)$geneSetSource)) 
		stop("Plotting a GO graph is only possible for MLP results based om geneSetSource 'GOBP', 'GOMF', or 'GOCC'")
	main <- if (is.null(main)) 
				"Go graph"
			else main
	require(GO.db)
	require(Rgraphviz)
	require(gplots)
	require(gmodels)
	require(gdata)
	require(gtools)
	require(GOstats)
	require(annotate)
	goids <- rownames(object)[1:nRow]
	ontology <- sub("GO", "", attributes(object)$geneSetSource)
	basicGraph <- GOGraph(goids, get(paste("GO", ontology, "PARENTS", 
							sep = "")))
	basicGraph <- removeNode("all", basicGraph)
	basicGraph <- removeNode(setdiff(nodes(basicGraph), rownames(object)), 
			basicGraph)
	basicGraph <- layoutGraph(basicGraph)
	pvalues <- object[nodes(basicGraph), "geneSetPValue"]
	names(pvalues) <- nodes(basicGraph)
	pvalues <- pvalues[!is.na(pvalues)]
	pvalues[pvalues == 0] <- min(pvalues[pvalues != 0])/10
	scores <- -log10(pvalues)
	scores[scores <= 0.1] <- 0.1
	nColors <- round(max(scores) * 10)
	gocolors <- colorpanel(nColors, low = "lightyellow", high = "olivedrab")
	nodeFillColor <- rep("white", length(nodes(basicGraph)))
	names(nodeFillColor) <- nodes(basicGraph)
	nodeFillColor[names(scores)] <- gocolors[trunc(scores * 10)]
	nodeRenderInfo(basicGraph) <- list(fill = nodeFillColor)
	inbin <- object$totalGeneSetSize
	names(inbin) <- rownames(object)
	onchip <- object$testedGeneSetSize
	names(onchip) <- rownames(object)
	allcounts <- matrix(ncol = 2, nrow = length(inbin), dimnames = list(c(names(inbin)), 
					c("onchip", "inbin")))
	allcounts[names(onchip), "onchip"] <- onchip
	allcounts[names(inbin), "inbin"] <- inbin
	counts <- apply(allcounts, 1, function(x) {
				paste(x[1], x[2], sep = " - ")
			})
	terms <- getGOTerm(nodes(basicGraph))
	goTerm <- terms[[1]]
	goTerm <- sapply(goTerm, function(x) {
				paste(substr(x, start = 1, stop = 30), substr(x, start = 31, 
								stop = 60), sep = "\\\n")
			})
	nodeLabel <- paste(nodes(basicGraph), goTerm[nodes(basicGraph)], 
			counts[nodes(basicGraph)], sep = "\\\n")
	names(nodeLabel) <- nodes(basicGraph)
	nodeRenderInfo(basicGraph) <- list(label = nodeLabel)
	edgeRenderInfo(basicGraph)$arrowhead <- "none"
	nodeRenderInfo(basicGraph) <- list(shape = "ellipse")
	nodeRenderInfo(basicGraph) <- list(cex = 0.5)
	nodeRenderInfo(basicGraph) <- list(lWidth = 60)
	nodeRenderInfo(basicGraph) <- list(labelJust = "c")
	graphRenderInfo(basicGraph)$main <- main
	renderGraph(basicGraph)
	legend(300, 30, legend = paste(c(" least", 
							" medium", " most"), " (scores ", round(max(scores)) * 
							c(2, 5, 8)/10, ")", sep = ""), fill = gocolors[round(max(scores)) * 
							c(2, 5, 8)], cex = 0.7)
}