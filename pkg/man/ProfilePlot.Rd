\name{ProfilePlot}
\alias{ProfilePlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of genes
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ProfilePlot(Gene, Comps = Comps, GeneExpr = geneMat, Clusters = NULL, cols = NULL, AddLegend = TRUE, names = NULL, margins = c(8.1, 4.1, 1.1, 6.5), extra, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Gene}{
%%     ~~Describe \code{Gene} here~~
}
  \item{Comps}{
%%     ~~Describe \code{Comps} here~~
}
  \item{GeneExpr}{
%%     ~~Describe \code{GeneExpr} here~~
}
  \item{Clusters}{
%%     ~~Describe \code{Clusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{AddLegend}{
%%     ~~Describe \code{AddLegend} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{margins}{
%%     ~~Describe \code{margins} here~~
}
  \item{extra}{
%%     ~~Describe \code{extra} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
function (Gene, Comps = Comps, GeneExpr = geneMat, Clusters = NULL, 
    cols = NULL, AddLegend = TRUE, names = NULL, margins = c(8.1, 
        4.1, 1.1, 6.5), extra, ...) 
{
    GenesComps = as.vector(GeneExpr[which(rownames(GeneExpr) == 
        Gene), colnames(GeneExpr) \%in\% Comps])
    Others = as.vector(GeneExpr[which(rownames(GeneExpr) == Gene), 
        !(colnames(GeneExpr) \%in\% Comps)])
    if (length(Others) == 0) {
        Continue = FALSE
    }
    else {
        Continue = TRUE
    }
    ylims = c(GenesComps, Others)
    par(mar = margins, xpd = TRUE)
    plot(type = "n", x = 0, y = 0, xlim = c(0, ncol(GeneExpr)), 
        ylim = c(min(ylims), max(ylims)), ylab = "Gene Expression", 
        xlab = " ", xaxt = "n", ...)
    points(x = seq(1, length(GenesComps)), y = GenesComps, pch = 19, 
        col = "red")
    segments(x0 = 1, y0 = mean(ylims[1:length(GenesComps)]), 
        x1 = length(GenesComps), y1 = mean(ylims[1:length(GenesComps)]), 
        lwd = 3)
    Indices = c(colnames(GeneExpr)[which(colnames(GeneExpr) \%in\% 
        Comps)], colnames(GeneExpr)[which(!colnames(GeneExpr) \%in\% 
        Comps)])
    axis(1, at = seq(0.5, (ncol(GeneExpr) - 0.5)), labels = Indices, 
        las = 2, cex.axis = 0.7, xlab = " ")
    axis(2, ylim = c(min(ylims), max(ylims)), ylab = "Gene Expression")
    if (Continue == TRUE) {
        points(x = seq(length(GenesComps) + 1, ncol(GeneExpr)), 
            y = Others, pch = 19, col = "blue")
        segments(x0 = length(GenesComps) + 1, y0 = mean(ylims[(length(GenesComps) + 
            1):ncol(GeneExpr)]), x1 = ncol(GeneExpr), y1 = mean(ylims[(length(GenesComps) + 
            1):ncol(GeneExpr)]), lwd = 3)
    }
    if (!(is.null(Clusters))) {
        for (i in 1:length(Clusters)) {
            temp = Clusters[[i]][!(Clusters[[i]] \%in\% Comps)]
            index = which(Indices \%in\% temp)
            values = as.vector(GeneExpr[which(rownames(GeneExpr) == 
                Gene), colnames(GeneExpr) \%in\% temp])
            points(x = index, y = values, pch = 19, col = cols[i])
        }
    }
    if (AddLegend == TRUE) {
        if (is.null(names)) {
            stop("Provide the names of the used methods to draw a legend")
        }
        labels = c("Shared")
        colslegend = c("red")
        if (Continue == TRUE) {
            colslegend = c(colslegend, "blue")
            labels = c(labels, "Others")
        }
        if (!(is.null(Clusters))) {
            colslegend = c(colslegend, cols)
            for (j in 1:length(Clusters)) {
                labels = c(labels, paste("Extra in", names[j], 
                  sep = " "))
            }
        }
        par(xpd = T, mar = margins)
        legend(length(ylims) + extra, mean(c(min(ylims), max(ylims))), 
            legend = labels, col = colslegend, pch = 19)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
