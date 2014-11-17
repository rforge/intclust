\name{FeaturesPlot}
\alias{FeaturesPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot of a selection of features
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FeaturesPlot(LeadCpds, OrderedCpds, Features, Data, Color, nrclusters = 7, cols = my_palette2, name = c("FP"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LeadCpds}{
%%     ~~Describe \code{LeadCpds} here~~
}
  \item{OrderedCpds}{
%%     ~~Describe \code{OrderedCpds} here~~
}
  \item{Features}{
%%     ~~Describe \code{Features} here~~
}
  \item{Data}{
%%     ~~Describe \code{Data} here~~
}
  \item{Color}{
%%     ~~Describe \code{Color} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
}
  \item{name}{
%%     ~~Describe \code{name} here~~
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
function (LeadCpds, OrderedCpds, Features, Data, Color, nrclusters = 7, 
    cols = my_palette2, name = c("FP")) 
{
    par(mar = c(6, 3, 0, 6))
    temp = OrderedCpds[which(!(OrderedCpds \%in\% LeadCpds))]
    AllCpds = c(LeadCpds, temp)
    if (all(AllCpds \%in\% rownames(Data))) {
        Data = t(Data)
    }
    x <- c(1:length(AllCpds))
    y <- c(1:length(Features))
    PlotData <- t(Data[Features, AllCpds])
    image(x, y, PlotData, col = c("gray90", "blue"), axes = FALSE)
    image(x[1:length(LeadCpds)], y, PlotData[1:length(LeadCpds), 
        ], col = c("gray90", "green"), add = TRUE)
    if (!(is.null(Color))) {
        Data1 <- Color$Clust
        ClustData1 = cutree(Data1, nrclusters)
        ordercolors = ClustData1[Data1$order]
        names(ordercolors) = Data1$order.lab
        ClustData1 = ClustData1[Data1$order]
        order = seq(1, nrclusters)
        for (k in 1:length(unique(ClustData1))) {
            select = which(ClustData1 == unique(ClustData1)[k])
            ordercolors[select] = order[k]
        }
        colors <- cols[ordercolors]
        names(colors) <- names(ordercolors)
    }
    mtext(colnames(PlotData), side = 4, at = c(1:ncol(PlotData)), 
        line = 0, las = 2, cex = 0.8)
    mtext(name, side = 2, line = 1, las = 0, cex = 1)
    mtext(rownames(PlotData), side = 1, at = c(1:nrow(PlotData)), 
        line = 0, las = 2, cex = 0.8, col = colors[AllCpds])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
