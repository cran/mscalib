\name{subset.mlist}
\alias{subset.mlist}
\title{ Subset mlist}
\description{ Return subsets of list elements which meet conditions. }
\usage{\method{subset}{mlist}(x,subset,...)}
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
