\name{FindElement}
\alias{FindElement}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find an element in a data structure
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FindElement(What, Object, Element = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{What}{
%%     ~~Describe \code{What} here~~
}
  \item{Object}{
%%     ~~Describe \code{Object} here~~
}
  \item{Element}{
%%     ~~Describe \code{Element} here~~
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
function (What, Object, Element = list()) 
{
    if (class(Object) == "data.frame") {
        if (What \%in\% colnames(Object)) {
            Element <- Object[, What]
        }
        else if (What \%in\% rownames(Oject)) {
            Element <- Object[What, ]
        }
    }
    if (class(Object) == "list") {
        for (i in 0:length(Object)) {
            if (i == 0) {
                Names = names(Object)
                if (What \%in\% Names) {
                  for (j in which(What == Names)) {
                    Element[length(Element) + 1] = Object[j]
                    names(Element)[length(Element)] = paste(What, 
                      "_", length(Element), sep = "")
                    return(Element)
                  }
                }
            }
            else if (class(Object[[i]])[1] == "list") {
                Element = FindElement(What, Object[[i]], Element = Element)
            }
        }
    }
    return(Element)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
