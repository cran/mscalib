\name{[[Assing.caliblist}
\alias{[[Assing.caliblist}
\alias{[[<-.caliblist}
\title{ Replace calibstat object in the caliblist.}
\description{ Replace a calibstat object in the caliblist with a different one. }
\usage{}
\arguments{
\item{i}{ index or name (info) of calibstat object to replace}
\item{value}{ object of class calibstat (e.g. calibrestat, calibextstat).}
\item{x}{ caliblist}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mvl)
 data(cal)
 res<-getintcalib(mvl,cal,error=300)
 res[[1]]<-res[[10]]
}
\keyword{misc}
