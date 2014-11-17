\name{ClusterDistribution}
\alias{ClusterDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Follow a cluster over the different methods
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ClusterDistribution(List, Selection, nrclusters, followMaxComps = FALSE, followClust = TRUE, fusionsLog = TRUE, WeightClust = TRUE, names = NULL, reverse = FALSE, Plot = TRUE, Table = TRUE, CompletePlot = FALSE, cols)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{Selection}{
%%     ~~Describe \code{Selection} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{followMaxComps}{
%%     ~~Describe \code{followMaxComps} here~~
}
  \item{followClust}{
%%     ~~Describe \code{followClust} here~~
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
  \item{reverse}{
%%     ~~Describe \code{reverse} here~~
}
  \item{Plot}{
%%     ~~Describe \code{Plot} here~~
}
  \item{Table}{
%%     ~~Describe \code{Table} here~~
}
  \item{CompletePlot}{
%%     ~~Describe \code{CompletePlot} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
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
function (List, Selection, nrclusters, followMaxComps = FALSE, 
    followClust = TRUE, fusionsLog = TRUE, WeightClust = TRUE, 
    names = NULL, reverse = FALSE, Plot = TRUE, Table = TRUE, 
    CompletePlot = FALSE, cols) 
{
    FoundClusters = list()
    ClusterDistribution.2 <- function(List, Selection, nrclusters, 
        followMaxComps, followClust, fusionsLog, WeightClust, 
        reverse) {
        if ((length(List) == 1 & (class(List[[1]]) == "Weighted" | 
            class(List[[1]]) == "CEC" | class(List[[1]]) == "WeightedSim")) & 
            reverse == TRUE) {
            List[[1]]$Results = rev(List[[1]]$Results)
        }
        else if (length(List) == 2) {
            if ((class(List[[2]]) == "Weighted" | class(List[[2]]) == 
                "CEC" | class(List[[2]]) == "WeightedSim") & 
                reverse == TRUE) {
                List[[2]]$Results = rev(List[[2]]$Results)
            }
        }
        Matrix = MatrixFunction(List, nrclusters, fusionsLog, 
            WeightClust, names)
        label.interest = unique(Matrix[1, which(colnames(Matrix) \%in\% 
            Selection)])
        for (i in 1:dim(Matrix)[1]) {
            temp = list()
            temp[[1]] = Selection
            names(temp)[1] = "Selection"
            clusternumbers = unique(Matrix[i, which(colnames(Matrix) \%in\% 
                Selection)])
            nr.clusters = length(clusternumbers)
            temp[[2]] = nr.clusters
            names(temp)[2] = "nr.clusters"
            min.together = min(table(Matrix[i, which(colnames(Matrix) \%in\% 
                Selection)]))
            max.together = max(table(Matrix[i, which(colnames(Matrix) \%in\% 
                Selection)]))
            nr.min.max.together = c(min.together, max.together)
            temp[[3]] = nr.min.max.together
            names(temp)[3] = "nr.min.max.together"
            min.perc.together <- min.together/length(Selection) * 
                100
            max.perc.together <- max.together/length(Selection) * 
                100
            perc.min.max.together = c(min.perc.together, max.perc.together)
            temp[[4]] = perc.min.max.together
            names(temp)[4] = "perc.min.max.together"
            temp[[5]] = list()
            names(temp)[5] = "AllClusters"
            for (a in 1:length(clusternumbers)) {
                temp[[5]][[a]] = list()
                names(temp[[5]])[a] = paste("Cluster", clusternumbers[a], 
                  sep = " ")
                temp[[5]][[a]][[1]] = clusternumbers[a]
                temp[[5]][[a]][[2]] = names(which(Matrix[i, ] == 
                  clusternumbers[a]))
                temp[[5]][[a]][[3]] = intersect(Selection, temp[[5]][[a]][[2]])
                temp[[5]][[a]][[4]] = temp[[5]][[a]][[2]][which(!(temp[[5]][[a]][[2]] \%in\% 
                  Selection))]
                names(temp[[5]][[a]]) = c("clusternumber", "Complete cluster", 
                  "Objects from original selection in this cluster", 
                  "Objects extra to this cluster")
            }
            if (followMaxComps == TRUE) {
                maxcluster = names(which(table(Matrix[i, which(colnames(Matrix) \%in\% 
                  Selection)]) == max(table(Matrix[i, which(colnames(Matrix) \%in\% 
                  Selection)]))))
                temp[[6]] = maxcluster
                names(temp)[6] = "Cluster with max Objects"
                complabels = names(Matrix[i, which(colnames(Matrix) \%in\% 
                  Selection & Matrix[i, ] == as.numeric(maxcluster))])
                temp[[7]] = complabels
                names(temp)[7] = "Complabels"
                complete.new.cluster = names(Matrix[i, which(Matrix[i, 
                  ] == as.numeric(maxcluster))])
                temp[[8]] = complete.new.cluster
                names(temp)[8] = "Complete.new.cluster"
                extra.new.cluster = complete.new.cluster[which(!(complete.new.cluster \%in\% 
                  Selection))]
                temp[[9]] = extra.new.cluster
                names(temp)[9] = "Extra.new.cluster"
            }
            if (followClust == TRUE) {
                temp[[6]] = label.interest
                names(temp)[6] = "Cluster"
                complabels = names(Matrix[i, which(colnames(Matrix) \%in\% 
                  Selection & Matrix[i, ] == as.numeric(label.interest))])
                temp[[7]] = complabels
                names(temp)[7] = "Complabels"
                complete.new.cluster = names(Matrix[i, which(Matrix[i, 
                  ] == as.numeric(label.interest))])
                temp[[8]] = complete.new.cluster
                names(temp)[8] = "Complete.new.cluster"
                extra.new.cluster = complete.new.cluster[which(!(complete.new.cluster \%in\% 
                  Selection))]
                temp[[9]] = extra.new.cluster
                names(temp)[9] = "Extra.new.cluster"
            }
            FoundClusters[[i]] = temp
        }
        FoundClusters[[length(FoundClusters) + 1]] = Matrix
        return(FoundClusters)
    }
    Found = ClusterDistribution.2(List, Selection, nrclusters, 
        followMaxComps, followClust, fusionsLog, WeightClust, 
        reverse)
    Matrix = Found[[length(Found)]]
    Found = Found[-length(Found)]
    if (is.null(names)) {
        for (j in 1:dim(Matrix)[1]) {
            names[j] = paste("Method", j, 1)
        }
    }
    if (Plot == TRUE) {
        if (followMaxComps == TRUE) {
            lab1 = c("Maximum of compounds of original cluster together")
            labelcluster = c()
            for (z in 1:length(Found)) {
                labelcluster = c(labelcluster, Found[[z]]$Cluster)
            }
        }
        else {
            lab1 = c("Number of compounds still in original cluster")
        }
        nrcluster = c()
        nrcomps = c()
        for (j in 1:length(Found)) {
            nrcluster = c(nrcluster, Found[[j]]$nr.clusters)
            nrcomps = c(nrcomps, length(Found[[j]]$Complabels))
        }
        plot(type = "n", x = 0, y = 0, xlim = c(1, length(Found)), 
            ylim = c(0, max(nrcluster, nrcomps) + 2), xlab = "", 
            ylab = "", xaxt = "n", yaxt = "n", cex.lab = 1.25)
        lines(x = seq(1, length(Found)), y = nrcluster, lty = 1, 
            col = "red", lwd = 1.5)
        points(x = seq(1, length(Found)), y = nrcluster, pch = 19, 
            col = "red", cex = 1.5)
        lines(x = seq(1, length(Found)), y = nrcomps, lty = 1, 
            col = "blue", lwd = 1.5)
        points(x = seq(1, length(Found)), y = nrcomps, pch = 19, 
            col = "blue", cex = 1.5)
        if (followMaxComps == TRUE) {
            text(seq(1, length(Found)), nrcomps, labelcluster, 
                cex = 1.5, pos = 1, col = "black", font = 2)
        }
        axis(1, at = seq(1, length(Found)), labels = names, las = 2, 
            cex = 1)
        axis(2, at = seq(0, max(nrcluster, nrcomps) + 2), labels = seq(0, 
            max(nrcluster, nrcomps) + 2), cex = 1)
        if (followMaxComps == TRUE) {
            legend(1, max(nrcluster, nrcomps) + 2.7, lty = c(1, 
                1, 0), pch = c(19, 19, 0), col = c("blue", "red", 
                "black"), legend = c(lab1, "Number of clusters original cluster divided amongst", 
                "Cluster number"), bty = "n", cex = 1.2)
        }
        else {
            legend(1, max(nrcluster, nrcomps) + 2.7, lty = c(1, 
                1), pch = c(19, 19), col = c("blue", "red"), 
                legend = c(lab1, "Number of clusters original cluster divided amongst"), 
                bty = "n", cex = 1.2)
        }
    }
    if (CompletePlot == TRUE) {
        windows()
        nrcluster = c()
        nrcomps = c()
        for (j in 1:length(Found)) {
            nrcluster = c(nrcluster, Found[[j]]$nr.clusters)
            nrcomps = c(nrcomps, length(Found[[j]]$Complabels))
        }
        plot(type = "n", x = 0, y = 0, xlim = c(1, length(Found)), 
            ylim = c(0, max(nrcluster, nrcomps) + 0.5), xlab = "", 
            ylab = "Number of Compounds", xaxt = "n", yaxt = "n")
        xnext = c()
        ynext = c()
        colorsp = c()
        for (m in 1:length(Found)) {
            howmany = length(Found[[m]]$AllClusters)
            if (m == 1) {
                xnext = rep(seq(1, length(Found))[m], howmany)
                for (g in 1:howmany) {
                  ynext = c(ynext, length(Found[[m]]$AllClusters[[g]][[3]]))
                  colorsp = c(colorsp, cols[Found[[m]]$AllClusters[[g]]$clusternumber])
                  if (length(ynext) > 1) {
                    for (t in 1:(length(ynext) - 1)) {
                      if (ynext[t] == ynext[length(ynext)]) {
                        ynext[length(ynext)] = ynext[length(ynext)] - 
                          0.3
                      }
                    }
                  }
                  points(x = xnext[g], y = ynext[length(ynext)], 
                    col = colorsp[length(colorsp)], pch = 19, 
                    cex = 1.25)
                  labelcluster = Found[[m]]$AllClusters[[g]]$clusternumber
                  position = 3
                  if (!(is.integer(ynext[length(ynext)]))) {
                    position = 1
                  }
                  text(xnext[g], ynext[length(ynext)], labelcluster, 
                    cex = 1.5, pos = position, col = "black", 
                    font = 2)
                }
                L1 = list()
                for (n in 1:howmany) {
                  L1[[n]] = Found[[m]]$AllClusters[[n]][[3]]
                }
            }
            else {
                xprev = xnext
                yprev = ynext
                ynext = c()
                colorsp = c()
                xnext = rep(seq(1, length(Found))[m], howmany)
                for (n in 1:howmany) {
                  ynext = c(ynext, length(Found[[m]]$AllClusters[[n]][[3]]))
                  colorsp = c(colorsp, cols[Found[[m]]$AllClusters[[n]]$clusternumber])
                  if (length(ynext) > 1) {
                    for (t in 1:(length(ynext) - 1)) {
                      if (ynext[t] == ynext[length(ynext)]) {
                        ynext[length(ynext)] = ynext[length(ynext)] - 
                          0.3
                      }
                    }
                  }
                  points(x = xnext[n], y = ynext[length(ynext)], 
                    col = colorsp[length(colorsp)], pch = 19, 
                    cex = 1.25)
                  labelcluster = Found[[m]]$AllClusters[[n]]$clusternumber
                  position = 3
                  if (!(is.integer(ynext[length(ynext)]))) {
                    position = 1
                  }
                  text(xnext[n], ynext[length(ynext)], labelcluster, 
                    cex = 1.5, pos = position, col = "black", 
                    font = 2)
                }
                L2 = L1
                L1 = list()
                for (n in 1:howmany) {
                  L1[[n]] = Found[[m]]$AllClusters[[n]][[3]]
                }
                for (q in 1:length(L1)) {
                  for (p in 1:length(L2)) {
                    if (length(which(L2[[p]] \%in\% L1[[q]])) != 
                      0) {
                      segments(x0 = xprev[p], y0 = yprev[p], 
                        x1 = xnext[q], y1 = ynext[q], col = colorsp[q], 
                        lwd = 2)
                    }
                  }
                }
            }
        }
        axis(1, at = seq(1, length(Found)), labels = names, las = 2, 
            cex = 1.5)
        axis(2, at = seq(0, max(nrcluster, nrcomps)), labels = seq(0, 
            max(nrcluster, nrcomps)), cex = 1.5)
        legend(1, max(nrcluster, nrcomps) + 1, pch = c(0), col = c("black"), 
            legend = c("Cluster number"), bty = "n", cex = 1.2)
    }
    if (Table == TRUE & Plot == TRUE) {
        SharedComps = Selection
        Extra = list()
        temp = c()
        for (a in 1:length(Found)) {
            SharedComps = intersect(SharedComps, Found[[a]]$Complabels)
        }
        for (a in 1:length(Found)) {
            Extra[[a]] = Found[[a]]$Complete.new.cluster[which(!(Found[[a]]$Complete.new.cluster \%in\% 
                SharedComps))]
            names(Extra)[a] = names[a]
            if (followMaxComps == TRUE) {
                names(Extra)[a] = paste(names(Extra)[a], "_", 
                  labelcluster[a], sep = "")
            }
            temp = c(temp, length(Extra[[a]]))
        }
        ExtraOr = Selection[which(!(Selection \%in\% SharedComps))]
        collength = max(length(SharedComps), length(ExtraOr), 
            temp)
        if (length(SharedComps) < collength) {
            spaces = collength - length(SharedComps)
            SharedComps = c(SharedComps, rep(" ", spaces))
        }
        if (length(ExtraOr) < collength) {
            spaces = collength - length(ExtraOr)
            ExtraOr = c(ExtraOr, rep(" ", spaces))
        }
        for (b in 1:length(Extra)) {
            if (length(Extra[[b]]) < collength) {
                spaces = collength - length(Extra[[b]])
                Extra[[b]] = c(Extra[[b]], rep(" ", spaces))
            }
        }
        table = data.frame(ExtraOr = ExtraOr, SharedComps = SharedComps)
        for (b in 2:length(Extra)) {
            table[1 + b] = Extra[[b]]
            colnames(table)[1 + b] = names(Extra)[b]
        }
    }
    names(Found) = names
    Found[[length(Found) + 1]] = table
    names(Found)[length(Found)] = "Table"
    return(Found)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
