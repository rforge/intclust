\name{PreparePathway}
\alias{PreparePathway}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Preparing a data set for pathway analysis if this is not the case
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PreparePathway(Object, GeneExpr, topG, sign)
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
  \item{sign}{
%%     ~~Describe \code{sign} here~~
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
function (Object, GeneExpr, topG, sign) 
{
    FoundGenes = NULL
    FoundComps = NULL
    FoundGenes = FindElement("Genes", Object)
    if (is.null(FoundGenes) | (is.list(FoundGenes) & length(FoundGenes) == 
        0)) {
        FoundComps = FindElement("Compounds", Object)
        if (is.null(FoundComps) | (is.list(FoundComps) & length(FoundComps) == 
            0)) {
            stop("Specify either the p-values of the genes or a selection of compounds to test for DE genes.")
        }
        pvalsgenes = list()
        FoundGenes = list()
        CompsP = list()
        TopDEP = list()
        for (i in 1:length(FoundComps)) {
            LeadCpds = FoundComps[[i]]$LeadCpds
            CompsP[[i]] = FoundComps[[i]]
            names(CompsP)[[i]] = paste("Compounds_", i, sep = "")
            if (is.null(LeadCpds)) {
                stop("In the Compounds element, specify an element LeadCpds")
            }
            group <- factor(ifelse(colnames(GeneExpr) \%in\% LeadCpds, 
                1, 0))
            if (class(GeneExpr)[1] == "ExpressionSet") {
                GeneExpr$LeadCmpds <- group
                DElead <- limmaTwoLevels(GeneExpr, "LeadCmpds")
                allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), 
                  resort.by = "logFC", sort.by = "p")
                if (is.null(allDE$ID)) {
                  allDE$Genes <- rownames(allDE)
                }
                else {
                  allDE$Genes = allDE$ID
                }
                if (is.null(topG)) {
                  topG = length(which(allDE$adj.P.Val <= sign))
                }
                TopDE <- allDE[1:topG, ]
                Genes <- list(TopDE, allDE)
                names(Genes) <- c("TopDE", "AllDE")
            }
            else {
                label.factor = factor(group)
                design = model.matrix(~label.factor)
                fit = lmFit(GeneExpr, design = design)
                fit = eBayes(fit)
                allDE = topTable(fit, coef = 2, adjust = "fdr", 
                  n = dim(GeneExpr)[1], resort.by = "logFC", 
                  sort.by = "p")
                if (is.null(allDE$ID)) {
                  allDE$Genes <- rownames(allDE)
                }
                else {
                  allDE$Genes = allDE$ID
                }
                if (is.null(topG)) {
                  topG = length(which(allDE$adj.P.Val <= sign))
                }
                TopDE <- allDE[1:topG, ]
                Genes <- list(TopDE, allDE)
                names(Genes) <- c("TopDE", "AllDE")
            }
            FoundGenes[[i]] = Genes
            TopDEP[[i]] = FoundGenes[[i]]
            names(TopDEP)[i] = paste("genes_", i, sep = "")
            names(FoundGenes) = paste("Genes_", i, sep = "")
            pvalsgenes[[i]] = Genes$AllDE$P.Value
            names(pvalsgenes)[i] = paste("pvals_", i, sep = "")
        }
    }
    else {
        pvalsgenes = list()
        TopDEP = list()
        for (i in 1:length(FoundGenes)) {
            TopDEP[[i]] = FoundGenes[[i]]
            names(TopDEP)[i] = paste("genes_", i, sep = "")
            pvalsgenes[[i]] = FoundGenes[[i]]$AllDE$P.Value
            names(pvalsgenes)[i] = paste("pvals_", i, sep = "")
        }
        FoundComps = FindElement("Compounds", Object)
        CompsP = list()
        for (i in 1:length(FoundComps)) {
            if (is.null(FoundComps[[i]]$LeadCpds)) {
                CompsP[[i]] = "No LeadCpds specified"
            }
            else {
                CompsP[[i]] = FoundComps[[i]]
                names(CompsP)[[i]] = paste("Compounds_", i, sep = "")
            }
        }
    }
    return(list(pvalsgenes = pvalsgenes, Compounds = CompsP, 
        Genes = TopDEP))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
