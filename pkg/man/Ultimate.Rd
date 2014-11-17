\name{Ultimate}
\alias{Ultimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function that performs any aggregated data function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Ultimate(List, distmeasure, NN = 20, alpha = 0.5, T = 20, t = 10, r = NULL, nrclusters = 7, nrclusterssep = c(7, 7), nrclustersseq = NULL, weight = NULL, Clustweight = 0.5, clust = "agnes", linkage = "ward", gap = FALSE, maxK = 50, IntClust = c("ADC", "ADECa", "ADECb", "ADECc", "WonM", "CECa", "CECb", "CECc", "WeightedClust", "WeightedSim", "SNFa", "SNFb", "SNFc"), fusionsLog = TRUE, WeightClust = TRUE, PlotCompare = FALSE, cols = my_palette2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{distmeasure}{
%%     ~~Describe \code{distmeasure} here~~
}
  \item{NN}{
%%     ~~Describe \code{NN} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{t}{
%%     ~~Describe \code{t} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
}
  \item{nrclusterssep}{
%%     ~~Describe \code{nrclusterssep} here~~
}
  \item{nrclustersseq}{
%%     ~~Describe \code{nrclustersseq} here~~
}
  \item{weight}{
%%     ~~Describe \code{weight} here~~
}
  \item{Clustweight}{
%%     ~~Describe \code{Clustweight} here~~
}
  \item{clust}{
%%     ~~Describe \code{clust} here~~
}
  \item{linkage}{
%%     ~~Describe \code{linkage} here~~
}
  \item{gap}{
%%     ~~Describe \code{gap} here~~
}
  \item{maxK}{
%%     ~~Describe \code{maxK} here~~
}
  \item{IntClust}{
%%     ~~Describe \code{IntClust} here~~
}
  \item{fusionsLog}{
%%     ~~Describe \code{fusionsLog} here~~
}
  \item{WeightClust}{
%%     ~~Describe \code{WeightClust} here~~
}
  \item{PlotCompare}{
%%     ~~Describe \code{PlotCompare} here~~
}
  \item{cols}{
%%     ~~Describe \code{cols} here~~
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
function (List, distmeasure, NN = 20, alpha = 0.5, T = 20, t = 10, 
    r = NULL, nrclusters = 7, nrclusterssep = c(7, 7), nrclustersseq = NULL, 
    weight = NULL, Clustweight = 0.5, clust = "agnes", linkage = "ward", 
    gap = FALSE, maxK = 50, IntClust = c("ADC", "ADECa", "ADECb", 
        "ADECc", "WonM", "CECa", "CECb", "CECc", "WeightedClust", 
        "WeightedSim", "SNFa", "SNFb", "SNFc"), fusionsLog = TRUE, 
    WeightClust = TRUE, PlotCompare = FALSE, cols = my_palette2, 
    ...) 
{
    if (class(List) != "list") {
        stop("Data must be of type list")
    }
    if (clust != "agnes" | linkage != "ward") {
        print("Only hierarchical clustering with WARD link is implemented. Perform your choice of clustering on the resulting\n\t\t\t\t\t\tfused matrix.")
        clust = "agnes"
        linkage = "ward"
    }
    SeparateClustering = list()
    for (i in 1:length(List)) {
        SeparateClustering[[i]] = Cluster(List[[i]], distmeasure[i], 
            clust, linkage, gap, maxK)
        names(SeparateClustering)[i] = paste("Data", i, sep = " ")
    }
    if (is.null(nrclusters) & is.null(nrclusterssep) & gap == 
        FALSE) {
        stop("Please specify a number of clusters in nrclusters and optional in nrclusterssep")
    }
    if (!(is.null(nrclusters)) & is.null(nrclusterssep)) {
        nrclusterssep = rep(nrclusters, length(List))
    }
    if (is.null(nrclusters) & is.null(nrclusterssep) & gap == 
        TRUE) {
        nrclusterssep = c()
        for (i in 1:length(SeparateClustering)) {
            temp = SeparateClustering[[i]]
            nrclusterssep[i] = temp$k[which(colnames(temp$k) == 
                "Tibs2001SEmax")]
        }
        nrclusters = ceiling(mean(nrclusterssep))
    }
    if (is.null(nrclustersseq) & ("ADECb" \%in\% IntClust | "WonM" \%in\% 
        IntClust | "CECb" \%in\% IntClust)) {
        stop("Method ADECb, WonM and CECb need specification of nrclustersseq.")
    }
    IntClustering = list()
    FollowedSelect = list()
    for (i in 1:length(IntClust)) {
        if (IntClust[i] == "SNFa") {
            IntClustering[[i]] = SNFa(List, distmeasure, NN, 
                alpha, T, clust, linkage)
            names(IntClustering)[i] = "SNFa"
            print("SNFa")
        }
        if (IntClust[i] == "SNFb") {
            IntClustering[[i]] = SNFb(List, distmeasure, NN, 
                alpha, T, clust, linkage)
            names(IntClustering)[i] = "SNFb"
            print("SNFb")
        }
        if (IntClust[i] == "SNFc") {
            IntClustering[[i]] = SNFc(List, distmeasure, NN, 
                alpha, T, clust, linkage)
            names(IntClustering)[i] = "SNFc"
            print("SNFc")
        }
        if (IntClust[i] == "ADC") {
            if (length(unique(distmeasure)) != 1) {
                stop("When using ADClust, distance measures must be the same for all data matrices ")
            }
            IntClustering[[i]] = ADC(List, unique(distmeasure), 
                clust, linkage)
            names(IntClustering)[i] = "ADC"
            print("ADC")
        }
        if (IntClust[i] == "ADECa") {
            if (length(unique(distmeasure)) != 1) {
                stop("When using ADECa, distance measures must be the same for all data matrices ")
            }
            IntClustering[[i]] = ADECa(List, unique(distmeasure), 
                t, r, nrclusters, clust, linkage)
            names(IntClustering)[i] = "ADECa"
            print("ADECa")
        }
        if (IntClust[i] == "ADECb") {
            if (length(unique(distmeasure)) != 1) {
                stop("When using ADECb, distance measures must be the same for all data matrices ")
            }
            IntClustering[[i]] = ADECb(List, unique(distmeasure), 
                nrclusters = nrclustersseq, clust, linkage)
            names(IntClustering)[i] = "ADECb"
            print("ADECb")
        }
        if (IntClust[i] == "ADECc") {
            if (length(unique(distmeasure)) != 1) {
                stop("When using ADECc, distance measures must be the same for all data matrices ")
            }
            IntClustering[[i]] = ADECc(List, unique(distmeasure), 
                t, r, nrclusters = nrclustersseq, clust, linkage)
            names(IntClustering)[i] = "ADECc"
            print("ADECc")
        }
        if (IntClust[i] == "WonM") {
            IntClustering[[i]] = WonM(List, distmeasure, nrclusters = nrclustersseq, 
                clust, linkage)
            names(IntClustering)[i] = "WonM"
            print("WonM")
        }
        if (IntClust[i] == "CECa") {
            IntClustering[[i]] = CECa(List, distmeasure, t, r, 
                nrclusters = nrclusterssep, weight, clust, linkage, 
                Clustweight)
            names(IntClustering)[i] = "CECa"
            print("CECa")
        }
        if (IntClust[i] == "CECb") {
            IntClustering[[i]] = CECb(List, distmeasure, nrclusters = nrclustersseq, 
                weight, clust, linkage, Clustweight)
            names(IntClustering)[i] = "CECb"
            print("CECb")
        }
        if (IntClust[i] == "CECc") {
            IntClustering[[i]] = CECc(List, distmeasure, t, r, 
                nrclusters = nrclustersseq, weight, clust, linkage, 
                Clustweight)
            names(IntClustering)[i] = "CECc"
            print("CECb")
        }
        if (IntClust[i] == "WeightedClust") {
            IntClustering[[i]] = WeightedClust(List, distmeasure, 
                weight, Clustweight, clust, linkage)
            names(IntClustering)[i] = "WeightedClust"
            print("Weighted Clustering")
        }
        if (IntClust[i] == "WeightedSim") {
            IntClustering[[i]] = WeightedSimClust(List, type = "data", 
                w = weight, clust = clust, linkage = linkage, 
                distmeasure = distmeasure, gap = FALSE, maxK = 50, 
                nrclusters = nrclusters, names = c("Data1", "Data2"), 
                AllClusters = TRUE)
            names(IntClustering)[i] = "WeightedSim"
            print("WeightedSim")
        }
    }
    if (PlotCompare == TRUE) {
        if (!(is.null(IntClustering$CECa))) {
            dev.new()
            ComparePlot(list(IntClustering$CECa), nrclusters, 
                cols, fusionsLog, names = c(seq(1, 0, -0.1)), 
                main = "CECa: Weights", ...)
        }
        if (!(is.null(IntClustering$CECb))) {
            dev.new()
            ComparePlot(list(IntClustering$CECb), nrclusters, 
                cols, fusionsLog, names = c(seq(1, 0, -0.1)), 
                main = "CECb: Weights", ...)
        }
        if (!(is.null(IntClustering$CECc))) {
            dev.new()
            ComparePlot(list(IntClustering$CECc), nrclusters, 
                cols, fusionsLog, names = c(seq(1, 0, -0.1)), 
                main = "CECc: Weights", ...)
        }
        if (!(is.null(IntClustering$WeightedClust))) {
            dev.new()
            ComparePlot(list(IntClustering$WeightedClust), nrclusters, 
                cols, fusionsLog, names = c(seq(1, 0, -0.1)), 
                main = "WeightedClust: Weights", ...)
        }
        if (!(is.null(IntClustering$WeightedSim))) {
            dev.new()
            ComparePlot(list(IntClustering$WeightedSim), nrclusters, 
                cols, fusionsLog, names = c(seq(1, 0, -0.1)), 
                main = "WeightedSim: Weights", ...)
        }
        dev.new()
        ComparePlot(IntClustering, nrclusters, cols, fusionsLog, 
            WeightClust = TRUE, names = names(IntClustering), 
            main = "All Methods", ...)
        dev.new()
        L = list()
        L[[1]] = SeparateClustering[[1]]
        names(L)[[1]] = "Data 1"
        for (k in 1:length(IntClustering)) {
            L[[k + 1]] = IntClustering[[k]]
            names(L)[k + 1] = names(IntClustering)[k]
        }
        L[[length(L) + 1]] = SeparateClustering[[2]]
        names(L)[[length(L)]] = "Data 2"
        ComparePlot(L, nrclusters, cols, fusionsLog, WeightClust = TRUE, 
            names = names(L), main = "All Methods", ...)
    }
    out = c(SeparateClustering[1], SeparateClustering[2], IntClustering)
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
