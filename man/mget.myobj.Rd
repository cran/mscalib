\name{mget.myobj}
\alias{mget.myobj}
\title{ Field Access}
\description{ Acces fields in object of class myobj }
\usage{\method{mget}{myobj}(object,attrn,...)}
\arguments{
\item{...}{}
\item{attrn}{ name of field (Attribute)}
\item{object}{ object of class myobj}
}
\author{Witold Wolski \email{wolski@molgen.mpg.de}}
\examples{
 data(mv1)
 mget(mv1)
 mget(mv1,"info")
}
\keyword{misc}
