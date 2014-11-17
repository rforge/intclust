\name{SharedSelection}
\alias{SharedSelection}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersection over multiple methods over genes and pathways for a selection of compounds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SharedSelection(DataLimma = NULL, DataMLP = NULL, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataLimma}{
%%     ~~Describe \code{DataLimma} here~~
}
  \item{DataMLP}{
%%     ~~Describe \code{DataMLP} here~~
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
function (DataLimma = NULL, DataMLP = NULL, names = NULL) 
{
    if (is.null(DataMLP)) {
        ResultShared = SharedSelectionLimma(DataLimma, names)
    }
    which = list()
    table = c()
    if (length(DataLimma) != length(DataMLP)) {
        stop("Unequal number of methods for limma and MLP")
    }
    else {
        nmethods = length(DataLimma)
    }
    if (is.null(names)) {
        for (j in 1:length(DataLimma)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    temp1g = c()
    temp1p = c()
    comps = c()
    pvalsg = c()
    pvalsp = c()
    for (i in 1:nmethods) {
        temp1g = c(temp1g, length(DataLimma[[i]]$Genes$ID))
        temp1p = c(temp1p, length(DataMLP[[i]][[3]]$descriptions))
        comps = c(comps, length(DataLimma[[i]]$Compounds))
        names(temp1g)[i] = names[i]
        names(temp1p)[i] = names[i]
        names(comps)[i] = paste("Ncomps", names[i], i, sep = " ")
        if (i == 1) {
            if (!(is.na(DataLimma[[i]])[1]) | !(is.na(DataMLP[[i]])[1])) {
                sharedcomps = DataLimma[[i]]$Compounds
                sharedgenes = DataLimma[[i]]$Genes$ID
                sharedpaths = DataMLP[[i]][[3]]$descriptions
                pvalsg = c(pvalsg, DataLimma[[i]]$Genes$adj.P.Val)
                pvalsp = c(pvalsp, DataMLP[[i]]$mean_p.value)
                nsharedcomps = length(DataLimma[[i]]$Compounds)
                nsharedgenes = length(DataLimma[[i]]$Genes$ID)
                nsharedpaths = length(DataMLP[[i]][[3]]$descriptions)
                names(nsharedgenes) = "nshared"
                names(nsharedpaths) = "nshared"
                names(nsharedcomps) = "nsharedcomps"
            }
        }
        else {
            sharedcomps = intersect(sharedcomps, DataLimma[[i]]$Compounds)
            sharedgenes = intersect(sharedgenes, DataLimma[[i]]$Genes$ID)
            sharedpaths = intersect(sharedpaths, DataMLP[[i]][[3]]$descriptions)
            nsharedcomps = length(intersect(sharedcomps, DataLimma[[i]]$Compounds))
            nsharedgenes = length(intersect(sharedgenes, DataLimma[[i]]$Genes$ID))
            nsharedpaths = length(intersect(sharedpaths, DataMLP[[i]][[3]]$descriptions))
            names(nsharedgenes) = "nshared"
            names(nsharedpaths) = "nshared"
            names(nsharedcomps) = "nsharedcomps"
        }
    }
    pvalsgenes = list()
    meanpvalsgenes = c()
    meanpvalspaths = c()
    pvalspaths = list()
    if (nsharedgenes != 0) {
        for (c in 1:nmethods) {
            pvalsg = c()
            for (g in sharedgenes) {
                if (!(is.na(DataLimma[[c]])[1])) {
                  pvalsg = c(pvalsg, DataLimma[[c]]$Genes$adj.P.Val[DataLimma[[c]]$Genes$ID == 
                    g])
                }
            }
            pvalsgenes[[c]] = pvalsg
        }
        for (g1 in 1:length(sharedgenes)) {
            pvalstemp = c()
            for (c in 1:nmethods) {
                if (!(is.na(DataLimma[[c]])[1])) {
                  pvalstemp = c(pvalstemp, pvalsgenes[[c]][[g1]])
                }
            }
            meanpvalsgenes = c(meanpvalsgenes, mean(pvalstemp))
        }
        pvalsgenes[[nmethods + 1]] = meanpvalsgenes
        names(pvalsgenes)[nmethods + 1] = "Mean pvals genes"
    }
    else {
        pvalsgenes = 0
    }
    if (nsharedpaths != 0) {
        for (c in 1:nmethods) {
            pvalsp = c()
            if (!(is.na(DataMLP[[c]])[1])) {
                for (p in sharedpaths) {
                  pvalsp = c(pvalsp, DataMLP[[c]][[3]][DataMLP[[c]][[3]]$descriptions == 
                    p, 3])
                }
            }
            pvalspaths[[c]] = pvalsp
        }
        for (p1 in 1:length(sharedpaths)) {
            pvalstemp1 = c()
            for (c in 1:nmethods) {
                if (!(is.na(DataMLP[[c]])[1])) {
                  pvalstemp1 = c(pvalstemp1, pvalspaths[[c]][[p1]])
                }
            }
            meanpvalspaths = c(meanpvalspaths, mean(pvalstemp1))
        }
        pvalspaths[[nmethods + 1]] = meanpvalspaths
        names(pvalspaths)[nmethods + 1] = "Mean pvals paths"
    }
    else {
        pvalpaths = 0
    }
    temp = rbind(cbind(temp1g, temp1p), cbind(nsharedgenes, nsharedpaths), 
        cbind(comps, comps), cbind(nsharedcomps, nsharedcomps))
    table = cbind(table, temp)
    which[[1]] = list(sharedcomps = sharedcomps, sharedgenes = sharedgenes, 
        pvalsgenes = pvalsgenes, sharedpaths = sharedpaths, pvalspaths = pvalspaths)
    for (i in 1:length(seq(1, dim(table)[2], 2))) {
        number = seq(1, dim(table)[2], 2)[i]
        colnames(table)[number] = c("G.Cluster")
        colnames(table)[number + 1] = c("P.Cluster")
    }
    ResultShared = list(Table = table, Which = which)
    return(ResultShared)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
