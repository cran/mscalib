\name{as.vector.calibrestat}
\alias{as.vector.calibrestat}
\title{ Coerces to vector}
\description{ coerces the calibrestat object into a vector. }
\usage{\method{as.vector}{calibrestat}(x,...)}
\arguments{
\item{...}{}
\item{x}{ calibrestat}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 res<-getrecalib(mv1)
 as.vector(res)
}
\keyword{misc}
