\name{PathwaysIter}
\alias{PathwaysIter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Iterations of the Pathways function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PathwaysIter(List, GeneExpr = geneMat, nrclusters = 7, method = c("limma", "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, niter = 10, fusionsLog = TRUE, WeightClust = TRUE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
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
  \item{ENTREZID}{
%%     ~~Describe \code{ENTREZID} here~~
}
  \item{geneSetSource}{
%%     ~~Describe \code{geneSetSource} here~~
}
  \item{top}{
%%     ~~Describe \code{top} here~~
}
  \item{topG}{
%%     ~~Describe \code{topG} here~~
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
function (List, GeneExpr = geneMat, nrclusters = 7, method = c("limma", 
    "MLP"), ENTREZID = GeneInfo[, 1], geneSetSource = "GOBP", 
    top = NULL, topG = NULL, GENESET = ListGO, sign = 0.05, niter = 10, 
    fusionsLog = TRUE, WeightClust = TRUE, names = NULL) 
{
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
        List.output = list()
        for (i in 1:niter) {
            print(paste("Iteration", i, sep = " "))
            comps.limma.mlp = Pathways.2(List[[1]], GeneExpr, 
                nrclusters, method, ENTREZID, geneSetSource, 
                top, topG, GENESET, sign = 0.05)
            List.output[[length(List.output) + 1]] = comps.limma.mlp
            names(List.output)[i] = paste("Iteration", i, sep = " ")
        }
    }
    else {
        List.output = list()
        for (i in 1:niter) {
            print(paste("Iteration", i, sep = " "))
            comps.limma.mlp = Pathways(List, GeneExpr, nrclusters, 
                method, ENTREZID, geneSetSource, top, topG, GENESET, 
                sign = 0.05, fusionsLog, WeightClust, names)
            List.output[[length(List.output) + 1]] = comps.limma.mlp
            names(List.output)[i] = paste("Iteration", i, sep = " ")
        }
    }
    return(List.output)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
