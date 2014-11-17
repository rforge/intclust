\name{MatrixFunction}
\alias{MatrixFunction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to order the outputs of many methods against a reference
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MatrixFunction(List, nrclusters = 7, fusionsLog = FALSE, WeightClust = FALSE, names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{List}{
%%     ~~Describe \code{List} here~~
}
  \item{nrclusters}{
%%     ~~Describe \code{nrclusters} here~~
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
function (List, nrclusters = 7, fusionsLog = FALSE, WeightClust = FALSE, 
    names = NULL) 
{
    matequal <- function(x, y) is.matrix(x) && is.matrix(y) && 
        dim(x) == dim(y) && all(x == y)
    ListNew = list()
    element = 0
    for (i in 1:length(List)) {
        if (class(List[[i]]) != "CEC" & class(List[[i]]) != "Weighted" & 
            class(List[[i]]) != "WeightedSim") {
            element = element + 1
            ListNew[[element]] = List[[i]]
        }
        else if (class(List[[i]]) == "CEC" | class(List[[i]]) == 
            "Weighted" | class(List[[i]]) == "WeightedSim") {
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
        names = seq(1, length(List), 1)
        for (i in 1:length(List)) {
            names[i] = paste("Method", i, sep = " ")
        }
    }
    ListClust = list()
    for (i in 1:length(List)) {
        ListClust[[i]] = List[[i]]$Clust
    }
    Clusters = list()
    Clusters = lapply(ListClust, cutree, nrclusters)
    xaxis = ListClust[[1]]$order
    ordercolors = Clusters[[1]][xaxis]
    order = seq(1, nrclusters)
    for (k in 1:length(unique(Clusters[[1]][xaxis]))) {
        select = which(Clusters[[1]][xaxis] == unique(Clusters[[1]][xaxis])[k])
        ordercolors[select] = order[k]
    }
    cols = unique(ordercolors)
    Ordered = list()
    autograph = list()
    for (i in cols) {
        autograph[[i]] = xaxis[which(ordercolors == i)]
    }
    for (j in 1:length(ListClust)) {
        print(j)
        ordercolorsj = Clusters[[j]][xaxis]
        order = seq(1, nrclusters)
        for (k in 1:length(unique(Clusters[[j]][xaxis]))) {
            select = which(Clusters[[j]][xaxis] == unique(Clusters[[j]][xaxis])[k])
            ordercolorsj[select] = order[k]
        }
        temp2 = ordercolorsj
        temp3 = xaxis
        fan = list()
        for (i in cols) {
            fan[[i]] = xaxis[which(temp2 == i)]
        }
        favors = matrix(0, length(autograph), length(fan))
        rownames(favors) = seq(1, length(autograph))
        colnames(favors) = seq(1, length(fan))
        for (a in 1:length(autograph)) {
            for (b in 1:length(fan)) {
                favorab = length(which(fan[[b]] \%in\% autograph[[a]]))/length(autograph[[a]])
                favors[a, b] = favorab
            }
        }
        tempfavors = favors
        matched = c(rep("Free", nrclusters))
        proposed = c(rep("No", nrclusters))
        Switches = c(rep("Open", nrclusters))
        proposals = matrix(0, length(autograph), length(fan))
        for (a in 1:dim(tempfavors)[1]) {
            for (b in 1:dim(tempfavors)[2]) {
                if (favors[a, b] == 1) {
                  matched[a] = b
                  proposed[b] = "Yes"
                  proposals[a, b] = 1
                  col = a
                  change = which(xaxis \%in\% fan[[b]])
                  temp3[change] = col
                  tempfavors[, b] = 0
                  tempfavors[a, ] = 0
                  Switches[a] = "Closed"
                }
            }
        }
        for (b in 1:dim(tempfavors)[2]) {
            if (length(which(tempfavors[, b] != 0)) == 1) {
                match = which(tempfavors[, b] != 0)
                test = which(tempfavors[match, ] == max(tempfavors[match, 
                  ]))[1]
                if (length(which(tempfavors[, test] != 0)) != 
                  1 | b \%in\% which(tempfavors[match, ] == max(tempfavors[match, 
                  ]))) {
                  matched[match] = b
                  proposed[b] = "Yes"
                  proposals[match, b] = 1
                  col = match
                  change = which(xaxis \%in\% fan[[b]])
                  temp3[change] = col
                  tempfavors[, b] = 0
                  tempfavors[match, ] = 0
                  Switches[match] = "Closed"
                }
            }
        }
        for (a in 1:dim(tempfavors)[1]) {
            if (length(which(tempfavors[a, ] != 0)) == 1) {
                propose = which(tempfavors[a, ] != 0)
                test = which(tempfavors[, propose] == max(tempfavors[, 
                  propose]))[1]
                if (length(which(tempfavors[test, ] != 0)) != 
                  1 | a \%in\% which(tempfavors[, propose] == max(tempfavors[, 
                  propose]))) {
                  matched[a] = propose
                  proposed[propose] = "Yes"
                  proposals[a, propose] = 1
                  col = a
                  change = which(xaxis \%in\% fan[[propose]])
                  temp3[change] = col
                  tempfavors[a, ] = 0
                  tempfavors[, propose] = 0
                  Switches[a] = "Closed"
                }
            }
        }
        Continue = TRUE
        if (length(which(matched == "Free")) == 0) {
            Continue = FALSE
        }
        while (length(which(matched == "Free")) != 0 | !(matequal(proposals[which(matched == 
            "Free"), ], matrix(1, length(which(matched == "Free")), 
            nrclusters))) | Continue != FALSE) {
            a = which(matched == "Free")[1]
            propose = which.max(tempfavors[a, ])
            if (tempfavors[a, propose] == 0) {
                if (length(which(matched == "Free")) == 1) {
                  Continue = FALSE
                }
                matched[a] = "Left"
            }
            else {
                if (proposed[propose] == "No") {
                  proposed[propose] = "Yes"
                  matched[a] = propose
                  proposals[a, propose] = 1
                  col = a
                  change = which(xaxis \%in\% fan[[propose]])
                  temp3[change] = col
                  tempfavors[a, propose] = 0
                  if (length(which(tempfavors[a, ] == 0)) == 
                    dim(tempfavors)[2]) {
                    Switches[a] = "Closed"
                    tempfavors[, propose] = 0
                    c = 1
                    while (c < a) {
                      if (Switches[c] != "Closed" & length(which(tempfavors[c, 
                        ] == 0)) == dim(tempfavors)[2]) {
                        Switches[c] = "Closed"
                        if (matched[c] == "Left") {
                          tempfavors[c, ] = 0
                        }
                        else {
                          tempfavors[, matched[c]] = 0
                        }
                        c = 1
                      }
                      else {
                        c = c + 1
                      }
                    }
                  }
                }
                else if (proposed[propose] == "Yes") {
                  if (favors[a, propose] > max(favors[which(matched == 
                    propose), propose]) & Switches[which(matched == 
                    propose)] == "Open") {
                    changeback = which(xaxis \%in\% fan[[propose]])
                    temp3[changeback] = xaxis[changeback]
                    matched[which(matched == propose)] = "Free"
                    matched[a] = propose
                    proposals[a, propose] = 1
                    col = a
                    change = which(xaxis \%in\% fan[[propose]])
                    temp3[change] = col
                    tempfavors[a, propose] = 0
                  }
                  else if (length(which(tempfavors[a, ] != 0)) == 
                    1) {
                    changeback = which(xaxis \%in\% fan[[propose]])
                    temp3[changeback] = xaxis[changeback]
                    matched[which(matched == propose)] = "Free"
                    matched[a] = propose
                    proposals[a, propose] = 1
                    col = a
                    change = which(xaxis \%in\% fan[[propose]])
                    temp3[change] = col
                    tempfavors[a, propose] = 0
                  }
                  else {
                    proposals[a, propose] = 1
                    tempfavors[a, propose] = 0
                  }
                }
            }
            if (length(which(matched == "Free")) == 0) {
                Continue = FALSE
            }
        }
        fusions = 0
        for (i in unique(matched)) {
            if (length(which(!(seq(1, nrclusters) \%in\% matched))) >= 
                1) {
                fusions = length(which(!(seq(1, nrclusters) \%in\% 
                  matched)))
            }
        }
        if (fusions != 0 & fusionsLog == FALSE) {
            print(paste("specify", fusions, "more color(s) and put fusionsLog equal to TRUE", 
                sep = " "))
        }
        premiumcol = c()
        for (i in 1:(fusions)) {
            premiumcol = c(premiumcol, length(matched) + i)
        }
        if ((length(which(matched == "Left")) != 0) | (length(which(proposed == 
            "No")) != 0)) {
            if (length(which(proposed == "No")) != 0) {
                for (i in 1:length(which(proposed == "No"))) {
                  Left = which(proposed == "No")[1]
                  maxLeft = which(favors[, Left] == max(favors[, 
                    Left]))
                  proposed[Left] = "Yes"
                  proposals[maxLeft, Left] = 1
                  col = premiumcol[i]
                  change = which(xaxis \%in\% fan[[Left]])
                  temp3[change] = col
                  tempfavors[, Left] = 0
                  tempfavors[maxLeft, ] = 0
                }
            }
            if (length(which(matched == "Left")) != 0) {
                for (i in 1:length(which(matched == "Left"))) Left = which(matched == 
                  "Left")[1]
                print(paste("Cluster", Left, "of the reference has found no suitable match.", 
                  sep = " "))
            }
        }
        Ordered[[j]] = temp3
    }
    Matrix = c()
    for (j in 1:length(Ordered)) {
        Matrix = rbind(Matrix, Ordered[[j]])
    }
    colnames(Matrix) = ListClust[[1]]$order.lab
    rownames(Matrix) = names
    return(Matrix)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
