\name{PathwaysCluster}
\alias{PathwaysCluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pathways for a specific cluster
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PathwaysCluster(Object, GeneExpr, topG, topP = NULL, method = c("limma", "MLP"), ENTREZIDs = NULL, geneSetSource = "GOBP", GENESET = ListGO, sign = 0.05, niter = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Object}{
%%     ~~Describe \code{Object} here~~
}
  \item{GeneExpr}{
%%     ~~Describe \code{GeneExpr} here~~
}
  \item{topG}{
%%     ~~Describe \code{topG} here~~
}
  \item{topP}{
%%     ~~Describe \code{topP} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{ENTREZIDs}{
%%     ~~Describe \code{ENTREZIDs} here~~
}
  \item{geneSetSource}{
%%     ~~Describe \code{geneSetSource} here~~
}
  \item{GENESET}{
%%     ~~Describe \code{GENESET} here~~
}
  \item{sign}{
%%     ~~Describe \code{sign} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
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
function (Object, GeneExpr, topG, topP = NULL, method = c("limma", 
    "MLP"), ENTREZIDs = NULL, geneSetSource = "GOBP", GENESET = ListGO, 
    sign = 0.05, niter = 10) 
{
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
    if (is.null(ENTREZIDs)) {
        ENTREZIDs = FoundGenes$Genes_1$All_Genes
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
    PreparedData = PreparePathway(Object, GeneExpr, topG, sign)
    pvalsgenes = PreparedData[[1]]
    Compounds = PreparedData[[2]]
    Pathways = list()
    for (k in 1:length(pvalsgenes)) {
        print(k)
        pvalsgenesRaw = as.vector(pvalsgenes[[k]])
        temp = list()
        for (j in 1:niter) {
            print(paste(k, ".", j, sep = ""))
            Paths = function(pvaluesRaw, path.method) {
                if (path.method == "MLP") {
                  names(pvaluesRaw) = ENTREZIDs
                  out.mlp <- MLP(geneSet = geneSet, geneStatistic = pvaluesRaw, 
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
                  if (is.null(topP)) {
                    output$ranked.genesets.table = ranked.genesets.table[ranked.genesets.table$p.values <= 
                      sign, ]
                  }
                  else {
                    output$ranked.genesets.table = ranked.genesets.table[1:topP, 
                      ]
                  }
                  nr.genesets = c(dim(ranked.genesets.table)[1], 
                    length(geneSet))
                  names(nr.genesets) = c("used.nr.genesets", 
                    "total.nr.genesets")
                  output$nr.genesets = nr.genesets
                }
                return(output)
            }
            temp[[j]] = Paths(pvalsgenesRaw, path.method)
            names(temp)[j] = paste("Iteration", j, sep = " ")
        }
        IntersectGenesets = function(list.output) {
            result.out = list()
            for (j in 1:length(list.output)) {
                cut = list.output[[j]]$ranked.genesets.table[list.output[[j]]$ranked.genesets.table[, 
                  2] <= sign, ]
                colnames(cut)[2] = paste("values.", j, sep = "")
                cut = cut[, c(1, 3, 2)]
                result.out[[j]] = cut
            }
            genesets.table.intersect = join_all(result.out, by = c("genesets", 
                "descriptions"), type = "inner")
            genesets.table.intersect$mean_p.value = apply(genesets.table.intersect[, 
                3:ncol(genesets.table.intersect)], 1, mean)
            result.out$genesets.table.intersect = genesets.table.intersect
            Paths = result.out$genesets.table.intersec
            return(Paths)
        }
        Pathways[[k]] = IntersectGenesets(temp)
        names(Pathways)[k] = paste("Pathways_", k, sep = "")
    }
    out = list(CompsP, FoundGenes, Pathways)
    names(out) = c("Compounds", "Genes", "Pathways")
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
