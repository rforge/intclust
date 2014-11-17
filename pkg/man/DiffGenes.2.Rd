\name{DiffGenes.2}
\alias{DiffGenes.2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DE for one result
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DiffGenes.2(Data, GeneExpr = geneMat, nrclusters = 7, method = "limma", sign = 0.05, top = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data}{
%%     ~~Describe \code{Data} here~~
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
  \item{sign}{
%%     ~~Describe \code{sign} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
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
function (Data, GeneExpr = geneMat, nrclusters = 7, method = "limma", 
    sign = 0.05, top = NULL) 
{
    if (method != "limma") {
        stop("Only the limma method is implemented to find differentially expressed genes")
    }
    DataClust = Data$Clust
    clusters = cutree(DataClust, nrclusters)
    Genes = list()
    if (is.null(top)) {
        top1 = FALSE
    }
    else {
        top1 = TRUE
    }
    for (i in unique(clusters)) {
        hc <- as.hclust(DataClust)
        OrderedCpds <- hc$labels[hc$order]
        temp = list()
        LeadCpds = colnames(GeneExpr)[which(clusters == i)]
        temp[[1]] = list(LeadCpds, OrderedCpds)
        names(temp[[1]]) = c("LeadCpds", "OrderedCpds")
        label = rep(0, dim(GeneExpr)[2])
        label[which(clusters == i)] = 1
        label.factor = factor(label)
        if (class(GeneExpr)[1] == "ExpressionSet") {
            GeneExpr$LeadCmpds <- label.factor
            DElead <- limmaTwoLevels(GeneExpr, "LeadCpds")
            allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), 
                resort.by = "logFC", sort.by = "p")
            if (top1 == TRUE) {
                result = list(allDE[1:top, ], allDE)
                names(result) = c("TopDE", "AllDE")
            }
            else if (top1 == FALSE) {
                top = length(which(allDE$adj.P.Val <= sign))
                result = list(allDE[1:top, ], allDE)
                names(result) = c("TopDE", "AllDE")
            }
        }
        else {
            design = model.matrix(~label.factor)
            fit = lmFit(GeneExpr, design = design)
            fit = eBayes(fit)
            allDE = topTable(fit, coef = 2, adjust = "fdr", resort.by = "logFC", 
                sort.by = "p")
            if (top1 == TRUE) {
                result = list(allDE[1:top, ], allDE)
                names(result) = c("TopDE", "AllDE")
            }
            else if (top1 == FALSE) {
                top = length(which(allDE$adj.P.Val <= sign))
                result = list(allDE[1:top, ], allDE)
                names(result) = c("TopDE", "AllDE")
            }
        }
        temp[[2]] = result
        names(temp) = c("Compounds", "Genes")
        Genes[[i]] = temp
        names(Genes)[i] = paste("Cluster", i, sep = " ")
    }
    return(Genes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
