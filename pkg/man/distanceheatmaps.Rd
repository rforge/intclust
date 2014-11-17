\name{distanceheatmaps}
\alias{distanceheatmaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal function for HeatMapCols
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
distanceheatmaps(Data1, Data2, names, nrclusters = 7)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Data1}{
%%     ~~Describe \code{Data1} here~~
}
  \item{Data2}{
%%     ~~Describe \code{Data2} here~~
}
  \item{names}{
%%     ~~Describe \code{names} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
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
function (Data1, Data2, names, nrclusters = 7) 
{
    ClustData1 = cutree(Data1, nrclusters)
    ClustData2 = cutree(Data2, nrclusters)
    ClustData1 = ClustData1[Data1$order]
    ClustData2 = ClustData2[Data2$order]
    trueorder1 = sort(Data1$order, index.return = TRUE)
    trueorder2 = sort(Data2$order, index.return = TRUE)
    ordercolors = ClustData1
    order = seq(1, nrclusters)
    for (k in 1:length(unique(ClustData1))) {
        select = which(ClustData1 == unique(ClustData1)[k])
        ordercolors[select] = order[k]
    }
    ordercolors2 = ClustData2
    for (k in 1:length(unique(ClustData2))) {
        select = which(ClustData2 == unique(ClustData2)[k])
        ordercolors2[select] = order[k]
    }
    ClustData1 = ordercolors[trueorder1$ix]
    ClustData2 = ordercolors2[trueorder2$ix]
    out = matrix(0, length(ClustData1), length(ClustData2))
    rownames(out) = names
    colnames(out) = names
    for (i in 1:length(names)) {
        focus = names[i]
        label = ClustData2[i]
        for (j in 1:length(names)) {
            other = names[j]
            found = FALSE
            k = 1
            while (found == FALSE & k <= nrclusters) {
                label2 = k
                if (other \%in\% names[ClustData1 == label2]) {
                  othercluster = names[ClustData1 == label2]
                  found = TRUE
                }
                k = k + 1
            }
            if (focus \%in\% othercluster) {
                out[j, i] = label
            }
        }
    }
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
