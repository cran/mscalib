\name{as.vector.calibintstat}
\alias{as.vector.calibintstat}
\title{ Coerces to vectors}
\description{ 'as.vector', a generic, attempts to coerce its argument into a   vector of mode 'mode' (the default is to coerce to whichever mode  is most convenient).  The attributes of 'x' are removed. }
\usage{\method{as.vector}{calibintstat}(x,mode="any")}
\arguments{
\item{mode}{ any}
\item{x}{ object of class calibintstat}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 data(cal)
 res<-getintcalib(mv1,cal,error=500,ppm=TRUE)
 as.vector(res)
}
\keyword{misc}
