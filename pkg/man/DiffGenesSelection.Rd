\name{DiffGenesSelection}
\alias{DiffGenesSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
DE for a selection of compounds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DiffGenesSelection(List, Selection, GeneExpr = geneMat, nrclusters = 7, method = "limma", sign = 0.05, top = NULL, fusionsLog = TRUE, WeightClust = TRUE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{Selection}{
%%     ~~Describe \code{Selection} here~~
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
function (List, Selection, GeneExpr = geneMat, nrclusters = 7, 
    method = "limma", sign = 0.05, top = NULL, fusionsLog = TRUE, 
    WeightClust = TRUE, names = NULL) 
{
    if (method != "limma") {
        stop("Only the limma method is implemented to find differentially expressed genes")
    }
    ListNew = list()
    element = 0
    for (i in 1:length(List)) {
        if (class(List[[i]]) != "CEC" & class(List[[i]]) != "Weighted") {
            element = element + 1
            ListNew[[element]] = List[[i]]
        }
        else if (class(List[[i]]) == "CEC" | class(List[[i]]) == 
            "Weighted") {
            ResultsClust = list()
            if (WeightClust == TRUE) {
                ResultsClust[[1]] = list()
                ResultsClust[[1]][[1]] = List[[i]]$Clust
                names(ResultsClust[[1]])[1] = "Clust"
                element = element + 1
                ListNew[[element]] = ResultsClust[[1]]
            }
            else {
                for (j in 1:length(List[[i]]$Results)) {
                  ResultsClust[[j]] = list()
                  ResultsClust[[j]][[1]] = List[[i]]$Results[[j]]
                  names(ResultsClust[[j]])[1] = "Clust"
                  element = element + 1
                  ListNew[[element]] = ResultsClust[[j]]
                }
            }
        }
    }
    List = ListNew
    if (is.null(names)) {
        for (j in 1:length(List)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    if (length(List) == 1) {
        ResultLimma = DiffGenes.2(List[[1]], GeneExpr, nrclusters, 
            method, sign, top)
    }
    else {
        if (is.null(top)) {
            top1 = FALSE
        }
        else {
            top1 = TRUE
        }
        Matrix = MatrixFunction(List, nrclusters, fusionsLog, 
            WeightClust, names)
        ResultLimma = list()
        for (k in 1:dim(Matrix)[1]) {
            maxcluster = names(which(table(Matrix[k, which(colnames(Matrix) \%in\% 
                Selection)]) == max(table(Matrix[k, which(colnames(Matrix) \%in\% 
                Selection)]))))
            hc <- as.hclust(List[[k]]$Clust)
            OrderedCpds <- hc$labels[hc$order]
            Genes = list()
            temp = list()
            LeadCpds = colnames(Matrix)[which(Matrix[k, ] == 
                maxcluster)]
            temp[[1]] = list(LeadCpds, OrderedCpds)
            names(temp[[1]]) = c("LeadCpds", "OrderedCpds")
            label = rep(0, dim(Matrix)[2])
            label[which(Matrix[k, ] == maxcluster)] = 1
            label.factor = factor(label)
            GeneExpr.2 = GeneExpr[, colnames(Matrix)]
            if (class(GeneExpr.2)[1] == "ExpressionSet") {
                GeneExpr.2$LeadCmpds <- label.factor
                DElead <- limmaTwoLevels(GeneExpr.2, "LeadCpds")
                allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), 
                  resort.by = "logFC", sort.by = "p")
                if (is.null(allDE$ID)) {
                  allDE$Genes <- rownames(allDE)
                }
                else {
                  allDE$Genes = allDE$ID
                }
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
                fit = lmFit(GeneExpr.2, design = design)
                fit = eBayes(fit)
                allDE = topTable(fit, coef = 2, n = dim(GeneExpr)[1], 
                  adjust = "fdr", resort.by = "logFC", sort.by = "p")
                if (is.null(allDE$ID)) {
                  allDE$Genes <- rownames(allDE)
                }
                else {
                  allDE$Genes = allDE$ID
                }
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
            ResultLimma[[k]] = temp
            names(ResultLimma)[k] = paste(names[k], k, ": Cluster", 
                maxcluster, sep = " ")
        }
    }
    return(ResultLimma)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
