\name{[[Assing.massvectorlist}
\alias{[[Assing.massvectorlist}
\alias{[[<-.massvectorlist}
\title{ Replace Parts of an Object}
\description{ Replace a massvector in the massvectorlist with a different one. }
\usage{}
\arguments{
\item{i}{ index or name (info) of massvector to replace}
\item{value}{ massvector}
\item{x}{ massvectorlist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 data(mv1)
 mvl[[10]]<-mv1
}
\keyword{misc}
