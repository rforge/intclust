\name{SharedLimma}
\alias{SharedLimma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersection over multiple methods over genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SharedLimma(DataLimma, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataLimma}{
%%     ~~Describe \code{DataLimma} here~~
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
function (DataLimma, names = NULL) 
{
    which = list()
    table = c()
    if (is.null(names)) {
        for (j in 1:length(DataLimma)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    nmethods = length(DataLimma)
    for (i in 1:length(DataLimma[[1]])) {
        name = paste("Cluster", i, sep = " ")
        temp1g = c()
        comps = c()
        for (j in 1:nmethods) {
            if (!(is.na(DataLimma[[j]][[i]])[1])) {
                temp1g = c(temp1g, length(DataLimma[[j]][[i]]$Genes$ID))
                comps = c(comps, length(DataLimma[[j]][[i]]$Compounds))
            }
            else {
                temp1g = c(temp1g, "-")
                comps = c(comps, "-")
            }
            names(temp1g)[j] = names[j]
            names(comps)[j] = names[j]
        }
        j = 1
        Continue = TRUE
        while (Continue == TRUE) {
            if (!(is.na(DataLimma[[j]][[i]])[1])) {
                sharedcomps = DataLimma[[j]][[i]]$Compounds
                sharedgenes = DataLimma[[j]][[i]]$Genes$ID
                nsharedcomps = length(DataLimma[[j]][[i]]$Compounds)
                nsharedgenes = length(DataLimma[[j]][[i]]$Genes$ID)
                names(nsharedgenes) = "nshared"
                names(nsharedcomps) = "nsharedcomps"
                Continue = FALSE
            }
            j = j + 1
        }
        if (nmethods >= 2) {
            for (j in 2:length(DataLimma)) {
                if (!(is.na(DataLimma[[j]][[i]])[1])) {
                  sharedcomps = intersect(sharedcomps, DataLimma[[j]][[i]]$Compounds)
                  sharedgenes = intersect(sharedgenes, DataLimma[[j]][[i]]$Genes$ID)
                  nsharedcomps = length(intersect(sharedcomps, 
                    DataLimma[[j]][[i]]$Compounds))
                  nsharedgenes = length(intersect(sharedgenes, 
                    DataLimma[[j]][[i]]$Genes$ID))
                  names(nsharedgenes) = "nshared"
                  names(nsharedcomps) = "nsharedcomps"
                }
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
                  if (!(is.na(DataLimma[[c]][[i]])[1])) {
                    pvalsg = c(pvalsg, DataLimma[[c]][[i]]$Genes$adj.P.Val[DataLimma[[c]][[i]]$Genes$ID == 
                      g])
                  }
                }
                pvalsgenes[[c]] = pvalsg
            }
            for (g1 in 1:length(sharedgenes)) {
                pvalstemp = c()
                for (c in 1:nmethods) {
                  if (!(is.na(DataLimma[[c]][[i]])[1])) {
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
        which[[i]] = list(sharedcomps = sharedcomps, sharedgenes = sharedgenes, 
            pvalsgenes = pvalsgenes)
        names(which)[i] = paste("Cluster", i, sep = " ")
        temp = c(temp1g, nsharedgenes, comps, nsharedcomps)
        table = cbind(table, temp)
        colnames(table)[i] = paste("G.Cluster", i, sep = " ")
    }
    ResultShared = list(Table = table, Which = which)
    return(ResultShared)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
