\name{SharedSelectionLimma}
\alias{SharedSelectionLimma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersection over multiple methods over pathways for a selection of compounds
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SharedSelectionLimma(DataMLP = NULL, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function (DataMLP = NULL, names = NULL) 
{
    which = list()
    table = c()
    if (is.null(names)) {
        for (j in 1:length(DataMLP)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    nmethods = length(DataMLP)
    temp1g = c()
    comps = c()
    pvalsg = c()
    for (i in 1:nmethods) {
        temp1g = c(temp1g, length(DataMLP[[i]]$Genes$ID))
        comps = c(comps, length(DataMLP[[i]]$Compounds))
        names(temp1g)[i] = names[i]
        names(comps)[i] = paste("Ncomps", names[i], i, sep = " ")
        if (i == 1) {
            if (!(is.na(DataMLP[[i]])[1])) {
                sharedcomps = DataMLP[[i]]$Compounds
                sharedgenes = DataMLP[[i]]$Genes$ID
                pvalsg = c(pvalsg, DataMLP[[i]]$Genes$adj.P.Val)
                nsharedcomps = length(DataMLP[[i]]$Compounds)
                nsharedgenes = length(DataMLP[[i]]$Genes$ID)
                names(nsharedgenes) = "nshared"
                names(nsharedcomps) = "nsharedcomps"
            }
        }
        else {
            sharedcomps = intersect(sharedcomps, DataMLP[[i]]$Compounds)
            sharedgenes = intersect(sharedgenes, DataMLP[[i]]$Genes$ID)
            nsharedcomps = length(intersect(sharedcomps, DataMLP[[i]]$Compounds))
            nsharedgenes = length(intersect(sharedgenes, DataMLP[[i]]$Genes$ID))
            names(nsharedgenes) = "nshared"
            names(nsharedcomps) = "nsharedcomps"
        }
    }
    pvalsgenes = list()
    meanpvalsgenes = c()
    if (nsharedgenes != 0) {
        for (c in 1:nmethods) {
            pvalsg = c()
            for (g in sharedgenes) {
                if (!(is.na(DataMLP[[c]])[1])) {
                  pvalsg = c(pvalsg, DataMLP[[c]]$Genes$adj.P.Val[DataMLP[[c]]$Genes$ID == 
                    g])
                }
            }
            pvalsgenes[[c]] = pvalsg
        }
        for (g1 in 1:length(sharedgenes)) {
            pvalstemp = c()
            for (c in 1:nmethods) {
                if (!(is.na(DataMLP[[c]])[1])) {
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
    temp = rbind(temp1g, nsharedgenes, comps, nsharedcomps)
    table = cbind(table, temp)
    which[[1]] = list(sharedcomps = sharedcomps, sharedgenes = sharedgenes, 
        pvalsgenes = pvalsgenes)
    ResultShared = list(Table = table, Which = which)
    return(ResultShared)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


