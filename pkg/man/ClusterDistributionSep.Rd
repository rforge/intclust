\name{ClusterDistributionSep}
\alias{ClusterDistributionSep}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Follow a selection of compounds over the different methods
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ClusterDistributionSep(ClusterData, Selection, K, Height = NULL, extraObjects = FALSE, OriginalExtra = FALSE, MaxExtra = FALSE, ref = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ClusterData}{
%%     ~~Describe \code{ClusterData} here~~
}
  \item{Selection}{
%%     ~~Describe \code{Selection} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{Height}{
%%     ~~Describe \code{Height} here~~
}
  \item{extraObjects}{
%%     ~~Describe \code{extraObjects} here~~
}
  \item{OriginalExtra}{
%%     ~~Describe \code{OriginalExtra} here~~
}
  \item{MaxExtra}{
%%     ~~Describe \code{MaxExtra} here~~
}
  \item{ref}{
%%     ~~Describe \code{ref} here~~
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
function (ClusterData, Selection, K, Height = NULL, extraObjects = FALSE, 
    OriginalExtra = FALSE, MaxExtra = FALSE, ref = FALSE) 
{
    if (!is.null(Height)) {
        K = NULL
    }
    data <- ClusterData$Clust
    if (class(Selection) == "character") {
        Selection.names = Selection
    }
    else if (class(Selection) == "numeric") {
        Selection.names <- attr(data$diss, "Labels")[Selection]
    }
    else {
        stop("Selection must be of type character of numeric")
    }
    data.hclust <- as.hclust(data)
    labels <- cutree(data.hclust, k = K, h = Height)
    if (class(Selection) == "character") {
        labels.interest = labels[which(names(labels) \%in\% Selection)]
    }
    else if (class(Selection) == "numeric") {
        labels.interest <- labels[Selection]
    }
    nr.clusters <- length(unique(labels.interest))
    min.together <- min(table(labels.interest))
    max.together <- max(table(labels.interest))
    min.perc.together <- min.together/length(Selection) * 100
    max.perc.together <- max.together/length(Selection) * 100
    clusterlabel.max <- as.numeric(names(table(labels.interest))[which(table(labels.interest) == 
        max.together)])
    complabels.max <- list()
    for (i in 1:length(clusterlabel.max)) {
        complabels.max[[i]] <- names(labels.interest)[labels.interest == 
            clusterlabel.max[i]]
    }
    ObjectsinClust = list()
    for (i in 1:length(unique(labels.interest))) {
        ObjectsinClust[[i]] = list()
        temp = ObjectsinClust[[i]]
        label = unique(labels.interest)[i]
        temp[[1]] = label
        names(temp)[1] = c("Cluster number:")
        temp[[2]] = names(labels.interest[which(labels.interest == 
            label)])
        names(temp)[2] = c("Objects from original selection in this cluster:")
        temp[[3]] = names(labels[which(labels == label)])
        names(temp)[3] = c("Complete cluster:")
        if (extraObjects == TRUE) {
            temp[[4]] = c(temp[[3]][!(temp[[3]] \%in\% temp[[2]])])
            names(temp)[4] = c("Objects extra to this cluster:")
        }
        ObjectsinClust[[i]] = temp
        names(ObjectsinClust)[i] = paste("Cluster", label, sep = " ")
    }
    new.cluster.includesmax <- list()
    for (i in 1:length(clusterlabel.max)) {
        new.cluster.includesmax[[i]] <- which(labels == clusterlabel.max[i])
    }
    if (is.null(ClusterData$weight)) {
        ClusterData$weight = 0
    }
    if (ClusterData$weight == 1 | OriginalExtra == TRUE & ref == 
        TRUE) {
        additional.in.cluster <- list()
        for (i in 1:length(unique(labels.interest))) {
            original.clusterlabel <- unique(labels.interest)[i]
            full.cluster <- names(labels[labels == original.clusterlabel])
            additional.in.cluster[[i]] <- full.cluster[!(full.cluster \%in\% 
                Selection.names)]
            names(additional.in.cluster)[i] <- paste("Cluster", 
                original.clusterlabel, sep = " ")
        }
    }
    else {
        additional.in.cluster <- NA
    }
    if (MaxExtra == TRUE) {
        extra.in.max.fused.cluster <- list()
        for (i in 1:length(clusterlabel.max)) {
            total.fused.cluster <- names(new.cluster.includesmax[[i]])
            extra.in.max.fused.cluster[[i]] <- total.fused.cluster[!(total.fused.cluster \%in\% 
                complabels.max[[i]])]
            names(extra.in.max.fused.cluster)[i] <- paste("Cluster", 
                clusterlabel.max[i], sep = " ")
        }
    }
    else {
        extra.in.max.fused.cluster <- NA
    }
    out = list(Selection = Selection.names, nr.clusters = nr.clusters, 
        nr.min.max.together = c(min.together, max.together), 
        perc.min.max.together = c(min.perc.together, max.perc.together), 
        AllClusters = ObjectsinClust, complabels.max = complabels.max, 
        new.cluster.includesmax = new.cluster.includesmax, additional.in.originalcluster = additional.in.cluster, 
        extra.in.max.fused.cluster = extra.in.max.fused.cluster)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
