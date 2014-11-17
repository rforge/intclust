\name{FindGenes}
\alias{FindGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find Genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FindGenes(DataLimma, names = NULL)
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
    FoundGenes = list()
    if (is.null(names)) {
        for (j in 1:length(DataLimma)) {
            names[j] = paste("Method", j, sep = " ")
        }
    }
    nrclusters = length(DataLimma[[1]])
    for (j in 1:nrclusters) {
        FoundGenes[[j]] = list()
    }
    for (i in 1:length(DataLimma)) {
        for (j in 1:nrclusters) {
            if (!(is.na(DataLimma[[i]][[j]]))[1]) {
                tempgenes = DataLimma[[i]][[j]]$Genes$TopDE$Genes[1:5]
                if (!(is.null(tempgenes))) {
                  tempgenes = tempgenes[!(is.na(tempgenes))]
                  for (k in 1:length(tempgenes)) {
                    if (!(tempgenes[k] \%in\% names(FoundGenes[[j]]))) {
                      FoundGenes[[j]][[length(FoundGenes[[j]]) + 
                        1]] = names[i]
                      names(FoundGenes[[j]])[length(FoundGenes[[j]])] = tempgenes[k]
                    }
                    else if (tempgenes[k] \%in\% names(FoundGenes[[j]])) {
                      found = which(names(FoundGenes[[j]]) == 
                        tempgenes[k])
                      FoundGenes[[j]][[found]] = c(FoundGenes[[j]][[found]], 
                        names[i])
                    }
                  }
                }
            }
        }
    }
    for (l in 1:nrclusters) {
        namesl = names(FoundGenes[[l]])
        if (!(is.null(namesl))) {
            for (m in l:nrclusters) {
                namesm = names(FoundGenes[[m]])
                Templist = FoundGenes[[m]]
                if (length(namesm) != 0) {
                  if (l != m) {
                    for (k in 1:length(namesm)) {
                      if (namesm[k] \%in\% namesl) {
                        found = which(namesl == namesm[k])
                        methods = Templist[[k]]
                        del = which(names(FoundGenes[[m]]) == 
                          namesm[k])
                        FoundGenes[[m]][[del]] = c()
                        for (a in 1:length(methods)) {
                          methods[a] = paste(methods[a], "_", 
                            m, sep = "")
                        }
                        FoundGenes[[l]][[found]] = c(FoundGenes[[l]][[found]], 
                          methods)
                      }
                    }
                  }
                }
            }
        }
    }
    for (i in 1:length(FoundGenes)) {
        names(FoundGenes)[i] = paste("Cluster", i, sep = " ")
    }
    return(FoundGenes)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
