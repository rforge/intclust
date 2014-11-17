\name{Geneset.intersect}
\alias{Geneset.intersect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the intersection of PathwaysIter function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Geneset.intersect(list.output, sign, names = NULL, seperatetables = FALSE, separatepvals = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list.output}{
%%     ~~Describe \code{list.output} here~~
}
  \item{sign}{
%%     ~~Describe \code{sign} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{seperatetables}{
%%     ~~Describe \code{seperatetables} here~~
}
  \item{separatepvals}{
%%     ~~Describe \code{separatepvals} here~~
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
function (list.output, sign, names = NULL, seperatetables = FALSE, 
    separatepvals = FALSE) 
{
    if (is.null(names)) {
        for (j in 1:length(list.output$"Iteration 1")) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    subsets = list()
    nmethods = length(list.output$"Iteration 1")
    for (i in 1:nmethods) {
        subsets[[i]] = list()
    }
    names(subsets) = names
    for (j in 1:length(list.output)) {
        name1 = names(list.output)[j]
        for (k in 1:nmethods) {
            name2 = names[k]
            subsets[[name2]][[name1]] = list.output[[name1]][[name2]]
        }
    }
    Intersect = list()
    for (i in 1:length(subsets)) {
        Method = subsets[[i]]
        Clusters = list()
        nclus = length(Method[[1]])
        for (j in 1:length(Method)) {
            name3 = paste("Iteration", j, sep = " ")
            for (k in 1:nclus) {
                name4 = paste("Cluster", k, sep = " ")
                if (!(is.na(Method[[name3]][[name4]])[1])) {
                  Clusters[[name4]][[name3]] = Method[[name3]][[name4]]
                }
                else {
                  Clusters[[name4]][[name3]] = NA
                }
            }
        }
        IntersectM = list()
        for (a in 1:length(Clusters)) {
            if (!(is.na(Clusters[[a]])[1])) {
                result.out = list()
                result.name = c()
                for (b in 1:length(Clusters[[a]])) {
                  if (b == 1) {
                    Compounds = Clusters[[a]][[1]]$Compounds
                    Genes = Clusters[[a]][[1]]$Genes
                  }
                  cut = Clusters[[a]][[b]]$Pathways$ranked.genesets.table[Clusters[[a]][[b]]$Pathways$ranked.genesets.table[, 
                    2] <= sign, ]
                  colnames(cut)[2] = paste("values.", b, sep = "")
                  cut = cut[, c(1, 3, 2)]
                  result.out[[b]] = cut
                  result.name = c(result.name, paste("genesettable", 
                    b, sep = ""))
                }
                names(result.out) = result.name
                genesets.table.intersect = join_all(result.out, 
                  by = c("genesets", "descriptions"), type = "inner")
                genesets.table.intersect$mean_p.value = apply(genesets.table.intersect[, 
                  3:ncol(genesets.table.intersect)], 1, mean)
                result.out$genesets.table.intersect = genesets.table.intersect
                if (separatepvals == FALSE) {
                  result.out$genesets.table.intersect = genesets.table.intersect[, 
                    c(1, 2, ncol(genesets.table.intersect))]
                }
                if (seperatetables == FALSE) {
                  result.out = result.out$genesets.table.intersect
                }
                newresult = list(Compounds = Compounds, Genes = Genes, 
                  Pathways = result.out)
                IntersectM[[a]] = newresult
                names(IntersectM)[a] = names(Clusters)[[a]]
            }
            else {
                IntersectM[[a]] = NA
                names(IntersectM)[a] = names(Clusters)[[a]]
            }
        }
        Intersect[[i]] = IntersectM
    }
    names(Intersect) = names
    return(Intersect)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
