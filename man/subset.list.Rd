\name{subset.list}
\alias{subset.list}
\title{ Subset mlist}
\description{ Return subsets of list elements which meet conditions. }
\usage{\method{subset}{list}(x,subset,...)}
\arguments{
\item{...}{}
\item{subset}{ logical expression.}
\item{x}{ object of class mlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 mvl<-subset(mvl,lengthmv>30)
}
\keyword{misc}
