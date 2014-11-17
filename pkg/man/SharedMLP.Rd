\name{SharedMLP}
\alias{SharedMLP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Intersection over multiple methods over pathways
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SharedMLP(DataMLP)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DataMLP}{
%%     ~~Describe \code{DataMLP} here~~
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
function (DataMLP) 
{
    which = list()
    table = c()
    if (is.null(names)) {
        for (j in 1:length(DataMLP)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    nmethods = length(DataMLP)
    for (i in 1:length(DataMLP[[1]])) {
        name = paste("Cluster", i, sep = " ")
        temp1g = c()
        temp1p = c()
        comps = c()
        for (j in 1:nmethods) {
            if (!(is.na(DataMLP[[j]][[i]])[1])) {
                temp1g = c(temp1g, length(DataMLP[[j]][[i]]$Genes$ID))
                temp1p = c(temp1p, length(DataMLP[[j]][[i]][[3]]$descriptions))
                comps = c(comps, length(DataMLP[[j]][[i]]$Compounds))
            }
            else {
                temp1g = c(temp1g, "-")
                temp1p = c(temp1p, "-")
                comps = c(comps, "-")
            }
            names(temp1g)[j] = names[j]
            names(temp1p)[j] = names[j]
            names(comps)[j] = paste("Ncomps", names[j], sep = " ")
        }
        j = 1
        Continue = TRUE
        while (Continue == TRUE) {
            if (!(is.na(DataMLP[[j]][[i]])[1])) {
                sharedcomps = DataMLP[[j]][[i]]$Compounds
                sharedgenes = DataMLP[[j]][[i]]$Genes$ID
                sharedpaths = DataMLP[[j]][[i]][[3]]$descriptions
                nsharedcomps = length(DataMLP[[1]][[i]]$Compounds)
                nsharedgenes = length(DataMLP[[1]][[i]]$Genes$ID)
                nsharedpaths = length(DataMLP[[1]][[i]][[3]]$descriptions)
                names(nsharedpaths) = "nshared"
                names(nsharedgenes) = "nshared"
                names(nsharedcomps) = "nsharedcomps"
                Continue = FALSE
            }
            j = j + 1
        }
        if (nmethods >= 2) {
            for (j in 2:length(DataLimma)) {
                if (!(is.na(DataMLP[[j]][[i]])[1])) {
                  sharedcomps = intersect(sharedcomps, DataMLP[[j]][[i]]$Compounds)
                  sharedgenes = intersect(sharedgenes, DataMLP[[j]][[i]]$Genes$ID)
                  sharedpaths = intersect(sharedpaths, DataMLP[[j]][[i]]$descriptions)
                  nsharedcomps = length(intersect(sharedcomp, 
                    DataMLP[[j]][[i]]$Compounds))
                  nsharedgenes = length(intersect(sharedgenes, 
                    DataMLP[[2]][[i]]$Genes$ID))
                  nsharedpaths = length(intersect(sharedgenes, 
                    DataMLP[[j]][[i]]$descriptions))
                  names(nsharedpaths) = "nshared"
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
                  if (!(is.na(DataMLP[[c]][[i]])[1])) {
                    pvalsg = c(pvalsg, DataMLP[[c]][[i]]$Genes$adj.P.Val[DataMLP[[c]][[i]]$Genes$ID == 
                      g])
                  }
                }
                print(pvalsg)
                pvalsgenes[[c]] = pvalsg
            }
            for (g1 in 1:length(sharedgenes)) {
                pvalstemp = c()
                for (c in 1:nmethods) {
                  if ((!(is.na(DataMLP[[c]][[i]])[1]))) {
                    pvalstemp = c(pvalstemp, pvalsgenes[[c]][[g1]])
                  }
                }
                meanpvalsgenes = c(meanpvalsgenes, mean(pvalstemp))
            }
            pvalsgenes[[nmethods + 1]] = meanpvalsgenes
        }
        else {
            pvalsgenes = 0
        }
        if (nsharedpaths != 0) {
            for (c in 1:nmethods) {
                pvalsp = c()
                for (p in sharedpaths) {
                  if (!(is.na(DataMLP[[c]][[i]])[1])) {
                    pvalsp = c(pvalsp, DataMLP[[c]][[i]][[3]][DataMLP[[c]][[i]][[3]]$descriptions == 
                      p, 3])
                  }
                }
                pvalspaths[[c]] = pvalsp
                names(pvalspaths)[c] = paste("Method", c, sep = " ")
            }
            for (p1 in 1:length(sharedpaths)) {
                pvalstemp1 = c()
                for (c in 1:nmethods) {
                  if (!(is.na(DataMLP[[c]][[i]])[1])) {
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
        which[[i]] = list(sharedcomps = sharedcomps, sharedgenes = sharedgenes, 
            pvalsgenes = pvalsgenes, sharedpaths = sharedpaths, 
            pvalspaths = pvalspaths)
        names(which)[i] = paste("Cluster", i, sep = " ")
        temp = rbind(cbind(temp1g, temp1p), cbind(nsharedgenes, 
            nsharedpaths), cbind(comps, comps), cbind(nsharedcomps, 
            nsharedcomps))
        table = cbind(table, temp)
    }
    for (i in 1:length(seq(1, dim(table)[2], 2))) {
        number = seq(1, dim(table)[2], 2)[i]
        colnames(table)[number] = paste("G.Cluster", i, sep = " ")
        colnames(table)[number + 1] = paste("P.Cluster", i, sep = " ")
    }
    ResultShared = list(Table = table, Which = which)
    return(ResultShared)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
