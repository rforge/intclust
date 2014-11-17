\name{ChooseFeatures}
\alias{ChooseFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Interactive plot for DE Genes and features
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ChooseFeatures(Interactive = TRUE, LeadCpds, ClusterResult, ClusterColors = NULL, BinData, Datanames = c("FP"), GeneExpr, topChar = 20, topG = 20, sign = 0.05, nrclusters = 7, cols = my_palette2, N = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Interactive}{
%%     ~~Describe \code{Interactive} here~~
}
  \item{LeadCpds}{
%%     ~~Describe \code{LeadCpds} here~~
}
  \item{ClusterResult}{
%%     ~~Describe \code{ClusterResult} here~~
}
  \item{ClusterColors}{
%%     ~~Describe \code{ClusterColors} here~~
}
  \item{BinData}{
%%     ~~Describe \code{BinData} here~~
}
  \item{Datanames}{
%%     ~~Describe \code{Datanames} here~~
}
  \item{GeneExpr}{
%%     ~~Describe \code{GeneExpr} here~~
}
  \item{topChar}{
%%     ~~Describe \code{topChar} here~~
}
  \item{topG}{
%%     ~~Describe \code{topG} here~~
}
  \item{sign}{
%%     ~~Describe \code{sign} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (Interactive = TRUE, LeadCpds, ClusterResult, ClusterColors = NULL, 
    BinData, Datanames = c("FP"), GeneExpr, topChar = 20, topG = 20, 
    sign = 0.05, nrclusters = 7, cols = my_palette2, N = 1) 
{
    OrInteractive = Interactive
    for (a in 1:length(BinData)) {
        if ((all(rownames(BinData[[a]]) \%in\% ClusterResult$Clust$order.lab))) {
            BinData[[a]] = t(BinData[[a]])
        }
    }
    if (Interactive == TRUE) {
        windows()
        Clusterplot(ClusterResult, ClusterColors, nrclusters, 
            cols)
        hc1 <- as.hclust(ClusterResult$Clust)
        ClusterSpecs <- list()
        ClusterSpecs = identify(hc1, N = N, MAXCLUSTER = ncol(BinData[[1]]), 
            function(j) ChooseFeatures(Interactive = FALSE, LeadCpds = colnames(BinData[[1]][, 
                j]), ClusterResult, ClusterColors = NULL, BinData, 
                Datanames, GeneExpr, topChar, topG, sign, nrclusters, 
                cols))
        names(ClusterSpecs) <- sapply(seq(1, N), FUN = function(x) paste("Choice", 
            x, sep = " "))
    }
    else {
        if (class(BinData) != "list") {
            stop("The binary data matrices must be put into a list")
        }
        for (i in 1:length(BinData)) {
            BinData[[i]] <- BinData[[i]][which(rowSums(BinData[[i]]) != 
                0 & rowSums(BinData[[i]]) != ncol(BinData[[i]])), 
                ]
        }
        cpdSet <- colnames(BinData[[1]])
        DistW <- ClusterResult$DistW
        Clust <- ClusterResult$Clust
        hc <- as.hclust(Clust)
        OrderedCpds <- hc$labels[hc$order]
        if (class(LeadCpds) == "character") {
            LeadCpds = list(LeadCpds)
        }
        Specs = list()
        for (i in 1:length(LeadCpds)) {
            Compounds = list(LeadCpds[[i]], OrderedCpds)
            names(Compounds) = c("LeadCpds", "OrderedCpds")
            group <- factor(ifelse(cpdSet \%in\% LeadCpds[[i]], 
                1, 0))
            Characteristics = list()
            for (j in 1:length(BinData)) {
                binMat = BinData[[j]]
                pFish <- apply(binMat, 1, function(x) fisher.test(table(x, 
                  group))$p.value)
                pFish <- sort(pFish)
                if (is.null(topChar)) {
                  topChar = length(which(pFish < 0.05))
                }
                Characteristics[[j]] <- names(pFish[0:topChar])
                names(Characteristics)[j] = Datanames[j]
            }
            if (class(GeneExpr)[1] == "ExpressionSet") {
                GeneExpr$LeadCmpds <- group
                DElead <- limmaTwoLevels(GeneExpr, "LeadCpds")
                allDE <- topTable(DElead, n = length(DElead@MArrayLM$genes$SYMBOL), 
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
                TopDE <- allDE[0:topG, ]
                Genes <- list(TopDE, allDE)
                names(Genes) <- c("TopDE", "AllDE")
            }
            else {
                label.factor = factor(group)
                design = model.matrix(~label.factor)
                fit = lmFit(GeneExpr, design = design)
                fit = eBayes(fit)
                allDE = topTable(fit, coef = 2, adjust = "fdr", 
                  n = nrow(GeneExpr), sort.by = "p")
                if (is.null(allDE$ID)) {
                  allDE$Genes <- rownames(allDE)
                }
                else {
                  allDE$Genes = allDE$ID
                }
                if (is.null(topG)) {
                  topG = length(which(allDE$adj.P.Val <= sign))
                }
                TopDE = allDE[0:topG, ]
                Genes <- list(TopDE, allDE)
                names(Genes) <- c("TopDE", "AllDE")
            }
            out = list(Compounds, Characteristics, Genes)
            names(out) = c("Compounds", "Characteristics", "Genes")
            Specs[[i]] = out
            names(Specs)[i] = paste("Choice", i, sep = " ")
        }
        if (OrInteractive == TRUE | length(Specs) == 1) {
            return(out)
        }
        else {
            return(Specs)
        }
    }
    return(ClusterSpecs)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
