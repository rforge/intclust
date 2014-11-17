\name{Pathways}
\alias{Pathways}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pathways for multiple functions
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Pathways(List, GeneExpr = geneMat, nrclusters = 7, method = c("limma", "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, fusionsLog = TRUE, WeightClust = TRUE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{GeneExpr}{
%%     ~~Describe \code{GeneExpr} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{ENTREZID}{
%%     ~~Describe \code{ENTREZID} here~~
}
  \item{geneSetSource}{
%%     ~~Describe \code{geneSetSource} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{topG}{
%%     ~~Describe \code{topG} here~~
}
  \item{GENESET}{
%%     ~~Describe \code{GENESET} here~~
}
  \item{sign}{
%%     ~~Describe \code{sign} here~~
}
  \item{fusionsLog}{
%%     ~~Describe \code{fusionsLog} here~~
}
  \item{WeightClust}{
%%     ~~Describe \code{WeightClust} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (List, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
    "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", 
    top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, fusionsLog = TRUE, 
    WeightClust = TRUE, names = NULL) 
{
    ListNew = list()
    element = 0
    for (i in 1:length(List)) {
        if (class(List[[i]]) != "CEC" & class(List[[i]]) != "Weighted" & 
            class(List[[i]]) != "WeightedSim") {
            element = element + 1
            ListNew[[element]] = List[[i]]
        }
        else if (class(List[[i]]) == "CEC" | class(List[[i]]) == 
            "Weighted" | class(List[[i]]) == "WeightedSim") {
            ResultsClust = list()
            if (WeightClust == TRUE) {
                ResultsClust[[1]] = list()
                ResultsClust[[1]][[1]] = List[[i]]$Clust
                names(ResultsClust[[1]])[1] = "Clust"
                element = element + 1
                ListNew[[element]] = ResultsClust[[1]]
                names(ListNew)[length(ListNew)] = c("WeightClust")
            }
            else {
                for (j in 1:length(List[[i]]$Results)) {
                  ResultsClust[[j]] = list()
                  ResultsClust[[j]][[1]] = List[[i]]$Results[[j]]
                  names(ResultsClust[[j]])[1] = "Clust"
                  element = element + 1
                  ListNew[[element]] = ResultsClust[[j]]
                  names(ListNew)[length(ListNew)] = names(List[[i]]$Results)[j]
                }
            }
        }
    }
    List = ListNew
    if (length(List) == 1) {
        ResultMLP = Pathways.2(List[[1]], GeneExpr, nrclusters, 
            method, ENTREZID, geneSetSource, top, GENESET, sign)
    }
    else {
        method.test = function(sign.method, path.method) {
            method.choice = FALSE
            if (sign.method == "limma" & path.method == "MLP") {
                method.choice = TRUE
            }
            if (method.choice == TRUE) {
                return(list(sign.method = sign.method, path.method = path.method))
            }
            else {
                stop("Incorrect choice of method.")
            }
        }
        method.out = method.test(method[1], method[2])
        sign.method = method.out$sign.method
        path.method = method.out$path.method
        if (length(ENTREZID) == 1) {
            ENTREZID = colnames(GeneExpr)
        }
        if ((class(GENESET) == "geneSetMLP")[1]) {
            geneSet <- GENESET
        }
        else {
            geneSet <- AnnotateEntrezIDtoGO(GeneInfo[, 1], database = c("ensembl", 
                "hsapiens_gene_ensembl"), attributes = c("entrezgene", 
                "go_id", "description"), filters = "entrezgene", 
                species = "Human")
        }
        if (is.null(top)) {
            top1 = FALSE
        }
        else {
            top1 = TRUE
        }
        if (is.null(names)) {
            for (j in 1:length(List)) {
                names[j] = paste("Method", j, sep = " ")
            }
        }
        MatrixClusters = MatrixFunction(List, nrclusters, fusionsLog, 
            WeightClust, names)
        ResultMLP = list()
        maxclus = 0
        for (k in 1:dim(MatrixClusters)[1]) {
            print(k)
            clusters = MatrixClusters[k, ]
            if (max(clusters) > maxclus) {
                maxclus = max(clusters)
            }
            DataPrepared <- try_default(PreparePathway(List[[k]], 
                GeneExpr, topG, sign), NULL, quiet = TRUE)
            if (is.null(DataPrepared)) {
                Temp = List[[k]]
                for (i in unique(clusters)) {
                  Compounds = list()
                  Compounds$LeadCpds = names(clusters)[which(clusters == 
                    i)]
                  Compounds$OrderedCpds = as.hclust(List[[k]]$Clust)$labels[as.hclust(List[[k]]$Clust)$order]
                  Temp[[i + 1]] = list(Compounds = Compounds)
                  names(Temp)[i + 1] = paste("Cluster", i, sep = " ")
                }
                DataPrepared <- PreparePathway(Temp, GeneExpr, 
                  topG, sign)
            }
            PathwaysResults = list()
            for (i in 1:length(DataPrepared$pvalsgenes)) {
                print(paste(k, i, sep = "."))
                temp = list()
                temp[[1]] = DataPrepared$Compounds[[i]]
                temp[[2]] = DataPrepared$Genes[[i]]
                pvalscluster = DataPrepared$pvalsgenes[[i]]
                if (path.method == "MLP") {
                  names(pvalscluster) = ENTREZID
                  out.mlp <- MLP(geneSet = geneSet, geneStatistic = pvalscluster, 
                    minGenes = 5, maxGenes = 100, rowPermutations = TRUE, 
                    nPermutations = 100, smoothPValues = TRUE, 
                    probabilityVector = c(0.5, 0.9, 0.95, 0.99, 
                      0.999, 0.9999, 0.99999), df = 9)
                  output = list()
                  ranked.genesets.table = data.frame(genesets = (rownames(out.mlp)), 
                    p.values = as.numeric(out.mlp$geneSetPValue), 
                    descriptions = out.mlp$geneSetDescription)
                  ranked.genesets.table$genesets = as.character(ranked.genesets.table$genesets)
                  ranked.genesets.table$descriptions = as.character(ranked.genesets.table$descriptions)
                  output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values <= 
                    0.05, ]
                  nr.genesets = c(dim(ranked.genesets.table)[1], 
                    length(geneSet))
                  names(nr.genesets) = c("used.nr.genesets", 
                    "total.nr.genesets")
                  output$nr.genesets = nr.genesets
                  temp[[3]] = output
                }
                names(temp) = c("Compounds", "Genes", "Pathways")
                PathwaysResults[[i]] = temp
                names(PathwaysResults)[i] = paste("Cluster", 
                  i, sep = " ")
            }
            ResultMLP[[k]] = PathwaysResults
        }
        names(ResultMLP) = names
        for (i in 1:length(ResultMLP)) {
            for (k in 1:length(ResultMLP[[i]])) {
                if (is.null(ResultMLP[[i]][[k]])[1]) {
                  ResultMLP[[i]][[k]] = NA
                  names(ResultMLP[[i]])[k] = paste("Cluster", 
                    k, sep = " ")
                }
            }
            if (length(ResultMLP[[i]]) != maxclus) {
                extra = maxclus - length(ResultMLP[[i]])
                for (j in 1:extra) {
                  ResultMLP[[i]][[length(ResultMLP[[i]]) + j]] = NA
                  names(ResultMLP[[i]])[length(ResultMLP[[i]])] = paste("Cluster", 
                    length(ResultMLP[[i]]), sep = " ")
                }
            }
        }
    }
    return(ResultMLP)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
